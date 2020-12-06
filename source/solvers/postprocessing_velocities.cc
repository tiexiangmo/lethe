
#include <deal.II/dofs/dof_tools.h>

#include <solvers/postprocessing_velocities.h>

template <int dim, typename VectorType, typename DofsType>
AverageVelocities<dim, VectorType, DofsType>::AverageVelocities()
  : average_calculation(false)
{}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                local_evaluation_point,
  const Parameters::PostProcessing &post_processing,
  const double &                    current_time,
  const double &                    time_step)
{
  const double epsilon      = 1e-6;
  const double initial_time = post_processing.initial_time;
  dt                        = time_step;

  // When averaging velocities begins
  if (current_time >= (initial_time - epsilon) && !average_calculation)
    {
      average_calculation = true;
      real_initial_time   = current_time;

      // Store the first dt value in case dt varies.
      dt_0 = dt;
    }

  // Calculate (u*dt) at each time step and accumulate the values
  velocity_dt.equ(dt, local_evaluation_point);
  sum_velocity_dt += velocity_dt;

  // Get the inverse of the total time with the first time step since
  // the weighted velocities are calculated with the first velocity when
  // total time = 0.
  inv_range_time = 1. / ((current_time - real_initial_time) + dt_0);

  // Calculate the average velocities.
  average_velocities.equ(inv_range_time, sum_velocity_dt);

  this->calculate_reynolds_stresses(local_evaluation_point);
}

// Since Trilinos vectors and block vectors data doesn't have the same
// structure, those vectors need to be processed in different ways.
// Same about the dimension, where 2D solutions are stored in (u,v,p) groups
// for vectors and [(u,v)][p] for block vectors and where 3D solutions
// are stored in (u,v,w,p) groups for vectors and [(u,v,w)][p] for
// block vectors.
template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stresses(
  const VectorType &local_evaluation_point)
{
  unsigned int begin_index, end_index;
  unsigned int (*k_index)(unsigned int) = [](unsigned int i) {
    return i + dim;
  };

  const TrilinosWrappers::MPI::Vector *local_solution, *local_average;
  TrilinosWrappers::MPI::Vector *      rns_dt, *rss_dt, *k_dt;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      begin_index    = local_evaluation_point.local_range().first;
      end_index      = local_evaluation_point.local_range().second;
      local_solution = &local_evaluation_point;
      local_average  = &average_velocities;
      rns_dt         = &reynolds_normal_stress_dt;
      rss_dt         = &reynolds_shear_stress_dt;
      k_dt           = &reynolds_normal_stress_dt;
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::BlockVector>)
    {
      n_dofs_per_vertex = dim;
      begin_index       = local_evaluation_point.block(0).local_range().first;
      end_index         = local_evaluation_point.block(0).local_range().second;
      local_solution    = &local_evaluation_point.block(0);
      local_average     = &average_velocities.block(0);
      rns_dt            = &reynolds_normal_stress_dt.block(0);
      rss_dt            = &reynolds_shear_stress_dt.block(0);
      k_dt              = &reynolds_normal_stress_dt.block(1);
      k_index           = [](unsigned int i) { return i / dim; };
    }

  for (unsigned int i = begin_index; i < end_index; i += n_dofs_per_vertex)
    {
      // u'u'*dt
      (*rns_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                     ((*local_solution)[i] - (*local_average)[i]) * dt;

      // v'v'*dt
      (*rns_dt)[i + 1] = ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                         ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                         dt;

      // u'v'*dt
      (*rss_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                     ((*local_solution)[i + 1] - (*local_average)[i + 1]) * dt;

      // k*dt = 1/2(u'u'+v'v')*dt (turbulence kinetic energy)
      // Note : k_dt and rns_dt are both pointers of reynolds_normal_stress_dt
      // for Trilinos vector (not block vectors)
      (*k_dt)[k_index(i)] = ((*rns_dt)[i] + (*rns_dt)[i + 1]) / 2;

      if (dim == 3)
        {
          // w'w'*dt
          (*rns_dt)[i + 2] =
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

          // v'w'*dt
          (*rss_dt)[i + 1] =
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

          // w'u'*dt
          (*rss_dt)[i + 2] =
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
            ((*local_solution)[i] - (*local_average)[i]) * dt;

          // k*dt = 1/2(u'u'+v'v'+w'w')*dt
          (*k_dt)[k_index(i)] = (*k_dt)[k_index(i)] + (*rns_dt)[i + 2] / 2;
        }
    }

  // Sum of all reynolds stresses during simulation.
  sum_reynolds_normal_stress_dt += reynolds_normal_stress_dt;
  sum_reynolds_shear_stress_dt += reynolds_shear_stress_dt;

  // Calculate the reynolds stresses.
  reynolds_normal_stresses.equ(inv_range_time, sum_reynolds_normal_stress_dt);
  reynolds_shear_stresses.equ(inv_range_time, sum_reynolds_shear_stress_dt);
}


template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::initialize_vectors(
  const DofsType &    locally_owned_dofs,
  const DofsType &    locally_relevant_dofs,
  const unsigned int &dofs_per_vertex,
  const MPI_Comm &    mpi_communicator)
{
  // Save the number of dofs per vertex. If solution vector is a block vector,
  // the variable will be later set to dim because only dofs related to the
  // velocity solution are desired in this particular case.
  n_dofs_per_vertex = dofs_per_vertex;

  // Reinitialisation of the average velocity and reynolds stress vectors
  // to get the right length.
  velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_velocities.reinit(locally_owned_dofs, mpi_communicator);
  get_av.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  // Reinitialize independent components of stress tensor vectors.
  reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_normal_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rns.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_shear_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rss.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
}

template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;
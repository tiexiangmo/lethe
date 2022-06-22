/* ---------------------------------------------------------------------
*
* Copyright (C) 2019 - 2019 by the Lethe authors
*
* This file is part of the Lethe library
*
* The Lethe library is free software; you can use it, redistribute
* it, and/or modify it under the terms of the GNU Lesser General
* Public License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
* The full text of the license can be found in the file LICENSE at
* the top level of the Lethe distribution.
*
* ---------------------------------------------------------------------

*
* Author: Lucka Barbeau, Shahab Golshan, Bruno Blais Polytechnique Montreal,
2021
*/

#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_linear_force.h>
#include <dem/particle_particle_nonlinear_force.h>
#include <dem/particle_wall_linear_force.h>
#include <dem/particle_wall_nonlinear_force.h>
#include <fem-dem/ib_particles_dem.h>

template <int dim>
void
IBParticlesDEM<dim>::initialize(
  const std::shared_ptr<Parameters::IBParticles<dim>> &p_nsparam,
  const MPI_Comm &                                     mpi_communicator_input,
  const std::vector<IBParticle<dim>> &                 particles)
{
  parameters       = p_nsparam;
  mpi_communicator = mpi_communicator_input;
  dem_particles    = particles;
  boundary_cells.resize(dem_particles.size());

  dem_parameters.model_parameters.particle_particle_contact_force_method =
    Parameters::Lagrangian::ModelParameters::ParticleParticleContactForceModel::
      hertz;
  dem_parameters.model_parameters.particle_wall_contact_force_method =
    Parameters::Lagrangian::ModelParameters::ParticleWallContactForceModel::
      nonlinear;
  particle_particle_contact_force_object =
    std::make_shared<ParticleParticleHertz<dim>>(dem_parameters);
  std::vector<types::boundary_id> boundary_index(0);
  particle_wall_contact_force_object =
    std::make_shared<ParticleWallNonLinearForce<dim>>(
      dem_parameters.boundary_conditions.boundary_translational_velocity,
      dem_parameters.boundary_conditions.boundary_rotational_speed,
      dem_parameters.boundary_conditions.boundary_rotational_vector,
      triangulation_cell_diameter,
      dem_parameters,
      boundary_index);
}
template <int dim>
void
IBParticlesDEM<dim>::update_particles(
  const std::vector<IBParticle<dim>> &particles,
  double                              time)
{
  dem_particles = particles;
  cfd_time      = time;
}

template <int dim>
void
IBParticlesDEM<dim>::update_contact_candidates()
{
  particles_contact_candidates.resize(dem_particles.size());

  double radius_factor = parameters->contact_search_radius_factor;

  for (auto &particle_one : dem_particles)
    {
      for (auto &particle_two : dem_particles)
        {
          if (particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;

              if ((particle_one_location - particle_two_location).norm() <
                  (particle_one.radius + particle_two.radius) * radius_factor)
                {
                  particles_contact_candidates[particle_one.particle_id].insert(
                    particle_two.particle_id);
                }
            }
        }
    }
}



template <int dim>
void
IBParticlesDEM<dim>::calculate_pp_contact_force(
  const double               dt_dem,
  std::vector<Tensor<1, 3>> &contact_force,
  std::vector<Tensor<1, 3>> &contact_torque)
{
  for (auto &particle_one : dem_particles)
    {
      for (auto particle_contact_candidates_id =
             particles_contact_candidates[particle_one.id].begin();
           particle_contact_candidates_id !=
           particles_contact_candidates[particle_one.id].end();
           ++particle_contact_candidates_id)
        {
          const auto &particle_contact_id = *particle_contact_candidates_id;
          auto &      particle_two        = dem_particles[particle_contact_id];
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              particle_particle_contact_info_struct<dim> contact_info;

              // Check if there is already information on the contact of these
              // to particles. If not initialize it in the contact map with 0
              // values.
              try
                {
                  contact_info = pp_contact_map[particle_one.particle_id]
                                               [particle_two.particle_id];
                }
              catch (...)
                {
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d]           = 0;
                      contact_info.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_info;
                }

              // Calculation of normal overlap
              double normal_overlap =
                (particle_one.radius + particle_two.radius) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                // This means that the adjacent particles are in contact
                {
                  Tensor<1, 3> normal_force;
                  Tensor<1, 3> tangential_force;
                  Tensor<1, 3> particle_one_tangential_torque;
                  Tensor<1, 3> particle_two_tangential_torque;
                  Tensor<1, 3> rolling_resistance_torque;

                  particle_particle_contact_force_object
                    ->calculate_IB_particle_particle_contact_force(
                      normal_overlap,
                      pp_contact_map[particle_one.particle_id]
                                    [particle_two.particle_id],
                      normal_force,
                      tangential_force,
                      particle_one_tangential_torque,
                      particle_two_tangential_torque,
                      rolling_resistance_torque,
                      particle_one,
                      particle_two,
                      particle_one_location,
                      particle_two_location,
                      dt_dem,
                      particle_one.radius,
                      particle_two.radius,
                      particle_one.mass,
                      particle_two.mass);



                  contact_force[particle_one.particle_id] -=
                    (normal_force + tangential_force);
                  contact_force[particle_two.particle_id] +=
                    (normal_force + tangential_force);

                  contact_torque[particle_one.particle_id] =
                    contact_torque[particle_one.particle_id] -
                    particle_one_tangential_torque + rolling_resistance_torque;
                  contact_torque[particle_two.particle_id] =
                    contact_torque[particle_two.particle_id] -
                    particle_two_tangential_torque - rolling_resistance_torque;
                }

              else
                {
                  // if the adjacent pair is not in contact anymore
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d]           = 0;
                      contact_info.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id].erase(
                    particle_two.particle_id);
                }
            }
        }
    }
}

template <int dim>
void
IBParticlesDEM<dim>::calculate_pp_lubrication_force(
  const double /*dt_dem*/,
  const double               h_max,
  const double               h_min,
  const double               mu,
  std::vector<Tensor<1, 3>> &lubrication_force,
  std::vector<Tensor<1, 3>> & /*lubrication_torque*/)
{
  using numbers::PI;
  // loop over all particles to find pair of close partilces
  for (auto &particle_one : dem_particles)
    {
      for (auto particle_contact_candidates_id =
             particles_contact_candidates[particle_one.id].begin();
           particle_contact_candidates_id !=
           particles_contact_candidates[particle_one.id].end();
           ++particle_contact_candidates_id)
        {
          const auto &particle_contact_id = *particle_contact_candidates_id;
          auto &      particle_two        = dem_particles[particle_contact_id];
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              Tensor<1, 3>     radial_vector;
              double           radial_velocity;
              Tensor<1, 3>     f_lub;

              // Calculation of normal overlap
              double gap =
                particle_one_location.distance(particle_two_location) -
                (particle_one.radius + particle_two.radius);
              radial_vector =
                tensor_nd_to_3d(particle_one.position - particle_two.position);
              if (gap > 0 and gap < h_max)
                // This means that the adjacent particles are very close but not
                // in contact
                {
                  // Limit the smallest gap calculated
                  if (gap < h_min)
                    {
                      gap = h_min;
                    }

                  // Evaluate the force
                  radial_velocity =
                    scalar_product(-radial_vector, particle_one.velocity) +
                    scalar_product(radial_vector, particle_two.velocity);
                  f_lub =
                    3 / 2 * PI * mu *
                      (particle_one.radius * 2 * particle_two.radius * 2 /
                       (particle_one.radius * 2 + particle_two.radius * 2)) *
                      (particle_one.radius * 2 * particle_two.radius * 2 /
                       (particle_one.radius * 2 + particle_two.radius * 2)) /
                      gap * radial_velocity * radial_vector /
                      radial_vector.norm() -
                    3 / 2 * PI * mu *
                      (particle_one.radius * 2 * particle_two.radius * 2 /
                       (particle_one.radius * 2 + particle_two.radius * 2)) *
                      (particle_one.radius * 2 * particle_two.radius * 2 /
                       (particle_one.radius * 2 + particle_two.radius * 2)) /
                      h_max * radial_velocity * radial_vector /
                      radial_vector.norm();
                  ;
                }

              lubrication_force[particle_one.particle_id] = f_lub;
              lubrication_force[particle_two.particle_id] = -f_lub;
            }
        }
    }
}

template <int dim>
void
IBParticlesDEM<dim>::update_particles_boundary_contact(
  const std::vector<IBParticle<dim>> &particles,
  const DoFHandler<dim> &             dof_handler,
  const Quadrature<dim - 1> &         face_quadrature_formula,
  const Mapping<dim> &                mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      // Clear the last boundary cell candidates.
      boundary_cells[p_i].clear();

      // Find the new cells that are at a boundary and in proximity of the
      // particle.
      auto cells_at_boundary = LetheGridTools::find_boundary_cells_in_sphere(
        dof_handler, particles[p_i].position, particles[p_i].radius * 1.5);

      // Loop over the cells at the boundary.
      for (unsigned int i = 0; i < cells_at_boundary.size(); ++i)
        {
          unsigned int      n_face_q_points = face_quadrature_formula.size();
          FEFaceValues<dim> fe_face_values(mapping,
                                           fe,
                                           face_quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_normal_vectors);
          // Loop over the faces of the cell at the boundary.
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              // Find the face at the boundary
              if (cells_at_boundary[i]->face(face_id)->at_boundary())
                {
                  fe_face_values.reinit(cells_at_boundary[i], face_id);
                  // Loop over the quadrature point of the face at the boundary
                  // to store information about the location and normals of the
                  // boundary.
                  for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                       ++f_q_point)
                    {
                      Tensor<1, dim> normal_vector =
                        -fe_face_values.normal_vector(f_q_point);
                      BoundaryCellsInfo boundary_information;
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_boundary =
                        fe_face_values.quadrature_point(f_q_point);
                      boundary_cells[p_i].push_back(boundary_information);
                    }
                }
            }
        }

      // Regroup the information of all the processor
      auto global_boundary_cell =
        Utilities::MPI::all_gather(this->mpi_communicator, boundary_cells[p_i]);
      boundary_cells[p_i].clear();
      for (unsigned int i = 0; i < global_boundary_cell.size(); ++i)
        {
          boundary_cells[p_i].insert(boundary_cells[p_i].end(),
                                     global_boundary_cell[i].begin(),
                                     global_boundary_cell[i].end());
        }
    }
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pw_contact_force(
  const double               dt_dem,
  std::vector<Tensor<1, 3>> &contact_force,
  std::vector<Tensor<1, 3>> &contact_torque)
{
  double wall_youngs_modulus = parameters->wall_youngs_modulus;
  double wall_poisson_ratio  = parameters->wall_poisson_ratio;
  double wall_rolling_friction_coefficient =
    parameters->wall_rolling_friction_coefficient;
  double wall_friction_coefficient = parameters->wall_friction_coefficient;
  double wall_restitution_coefficient =
    parameters->wall_restitution_coefficient;

  // Loop over the particles
  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      unsigned int best_index     = UINT_MAX;
      // For each particle loop over the point and normal identified as
      // potential contact candidate.
      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          // find the best candidate (the closest point).
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      // Do the particle wall contact calculation with the best candidate.
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];

          auto boundary_cell_information = boundary_cell;
          particle_wall_contact_info_struct<dim> contact_info;

          // Check if there is already information on the contact between this
          // particle and this boundary contact point. If not initialize the
          // contact history with 0 values.
          try
            {
              contact_info =
                pw_contact_map[particle.particle_id][boundary_index];
            }
          catch (...)
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d]           = 0;
                  contact_info.tangential_relative_velocity[d] = 0;
                }

              // BB temporary fix for unused variables
              contact_info.global_face_id           = 0;
              contact_info.boundary_id              = 0;
              contact_info.normal_overlap           = 0;
              contact_info.normal_relative_velocity = 0;
              // End BB

              pw_contact_map[particle.particle_id][boundary_index] =
                contact_info;
            }

          Tensor<1, 3> normal =
            tensor_nd_to_3d(boundary_cell_information.normal_vector);
          auto point_on_boundary = boundary_cell_information.point_on_boundary;


          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance

          Point<3> particle_position_3d = point_nd_to_3d(particle.position);
          Point<3> point_on_boundary_3d = point_nd_to_3d(point_on_boundary);

          Tensor<1, 3> point_to_particle_vector =
            particle_position_3d - point_on_boundary_3d;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            particle_wall_contact_force_object->find_projection(
              point_to_particle_vector, normal);

          // Find the normal overlap
          double normal_overlap = particle.radius - (projected_vector.norm());

          if (normal_overlap > 0)
            {
              // Do the calculation to evaluate the particle wall contact force.
              contact_info.normal_overlap = normal_overlap;

              Tensor<1, 3> normal_force;
              Tensor<1, 3> tangential_force;
              Tensor<1, 3> tangential_torque;
              Tensor<1, 3> rolling_resistance_torque;

              pw_contact_map[particle.particle_id][boundary_index]
                .normal_overlap = normal_overlap;
              pw_contact_map[particle.particle_id][boundary_index]
                .normal_vector = normal;
              pw_contact_map[particle.particle_id][boundary_index].boundary_id =
                0;

              particle_wall_contact_force_object
                ->calculate_IB_particle_wall_contact_force(
                  pw_contact_map[particle.particle_id][boundary_index],
                  normal_force,
                  tangential_force,
                  tangential_torque,
                  rolling_resistance_torque,
                  particle,
                  wall_youngs_modulus,
                  wall_poisson_ratio,
                  wall_restitution_coefficient,
                  wall_friction_coefficient,
                  wall_rolling_friction_coefficient,
                  dt_dem,
                  particle.mass,
                  particle.radius);
              // Updating the force of particles in the particle handler
              contact_force[particle.particle_id] +=
                normal_force + tangential_force;
              // Updating the torque acting on particles
              contact_torque[particle.particle_id] +=
                tangential_torque + rolling_resistance_torque;
            }
          else
            {
              // Set to 0 the tangential overlap if the particle is not in
              // contact with the wall anymore.
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d] = 0;
                }
            }
        }
    }
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pw_lubrication_force(
  const double /*dt_dem*/,
  const double               h_max,
  const double               h_min,
  const double               mu,
  std::vector<Tensor<1, 3>> &lubrication_force,
  std::vector<Tensor<1, 3>> & /*lubrication_torque*/)
{
  using numbers::PI;

  // Loop over the particles
  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      unsigned int best_index     = UINT_MAX;
      // For each particle loop over the point and normal identified as
      // potential contact candidate.
      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          // find the best candidate (the closest point).
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      // Do the particle wall lubrication calculation with the best candidate.
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];

          auto         boundary_cell_information = boundary_cell;
          Tensor<1, 3> normal =
            tensor_nd_to_3d(boundary_cell_information.normal_vector);
          auto point_on_boundary = boundary_cell_information.point_on_boundary;

          // Calculation of gap
          Point<3> particle_position_3d = point_nd_to_3d(particle.position);
          Point<3> point_on_boundary_3d = point_nd_to_3d(point_on_boundary);

          Tensor<1, 3> point_to_particle_vector =
            particle_position_3d - point_on_boundary_3d;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            particle_wall_contact_force_object->find_projection(
              point_to_particle_vector, normal);


          double gap = (projected_vector.norm()) - particle.radius;

          Tensor<1, 3> radial_vector;
          double       radial_velocity;
          Tensor<1, 3> f_lub;


          radial_vector = projected_vector;
          if (gap > 0 and gap < h_max)
            // This means that the particles is very close to the wall but not
            // in contact
            {
              // Limit the smallest gap calculated
              if (gap < h_min)
                {
                  gap = h_min;
                }
              // Evaluate the force
              radial_velocity =
                scalar_product(-radial_vector, particle.velocity);
              f_lub = 3 / 2 * PI * mu * (particle.radius) * (particle.radius) /
                        gap * radial_velocity * radial_vector /
                        radial_vector.norm() -
                      3 / 2 * PI * mu * (particle.radius) * (particle.radius) /
                        h_max * radial_velocity * radial_vector /
                        radial_vector.norm();
            }
          lubrication_force[particle.particle_id] = f_lub;
        }
    }
}

template <int dim>
void
IBParticlesDEM<dim>::integrate_particles_motion(const double dt,
                                                const double h_max,
                                                const double h_min,
                                                const double rho,
                                                const double mu)
{
  // Initialize local containers and physical variables
  using numbers::PI;
  double dt_dem = dt / parameters->coupling_frequency;

  std::vector<Tensor<1, 3>> contact_force(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_wall_force(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_wall_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> current_fluid_force(dem_particles.size());
  std::vector<Tensor<1, 3>> current_fluid_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrication_force(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrication_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrication_wall_force(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrication_wall_torque(dem_particles.size());

  std::vector<Tensor<1, 3>> velocity(dem_particles.size());
  std::vector<Point<dim>>   position(dem_particles.size());

  // Local time for the dem step.
  double t = 0;
  // The gravitational acceleration.
  Tensor<1, 3> g;
  this->parameters->f_gravity->set_time(cfd_time);
  // The gravitational force on the particle.
  Tensor<1, 3> gravity;

  // Initialize the particles
  for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
    {
      dem_particles[p_i].position  = dem_particles[p_i].previous_positions[0];
      dem_particles[p_i].velocity  = dem_particles[p_i].previous_velocity[0];
      dem_particles[p_i].omega     = dem_particles[p_i].previous_omega[0];
      dem_particles[p_i].impulsion = 0;
      dem_particles[p_i].omega_impulsion         = 0;
      dem_particles[p_i].contact_impulsion       = 0;
      dem_particles[p_i].omega_contact_impulsion = 0;
      // Initialized the gravity at the particle position.
      g[0] = this->parameters->f_gravity->value(dem_particles[p_i].position, 0);
      g[1] = this->parameters->f_gravity->value(dem_particles[p_i].position, 1);
      if (dim == 3)
        g[2] =
          this->parameters->f_gravity->value(dem_particles[p_i].position, 2);
    }


  // Integrate with the sub_time_step
  while (t + 0.5 * dt_dem < dt)
    {
      // Initialize vector to contain the RK 4 step
      std::vector<Tensor<1, 3>> last_velocity;
      std::vector<Point<dim>>   last_position;
      std::vector<Tensor<1, 3>> last_omega;
      last_velocity = std::vector<Tensor<1, 3>>(dem_particles.size());
      last_position = std::vector<Point<dim>>(dem_particles.size());
      last_omega    = std::vector<Tensor<1, 3>>(dem_particles.size());

      std::vector<std::vector<Tensor<1, 3>>> k_velocity =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));
      std::vector<std::vector<Tensor<1, dim>>> k_position =
        std::vector<std::vector<Tensor<1, dim>>>(
          dem_particles.size(), std::vector<Tensor<1, dim>>(4));
      std::vector<std::vector<Tensor<1, 3>>> k_omega =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));
      std::vector<std::vector<Tensor<1, 3>>> k_impulsion =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));
      std::vector<std::vector<Tensor<1, 3>>> k_omega_impulsion =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));
      std::vector<std::vector<Tensor<1, 3>>> k_contact_impulsion =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));
      std::vector<std::vector<Tensor<1, 3>>> k_omega_contact_impulsion =
        std::vector<std::vector<Tensor<1, 3>>>(dem_particles.size(),
                                               std::vector<Tensor<1, 3>>(4));

      // Solve each of the 4 step of the RK method
      for (unsigned int step = 0; step < 4; ++step)
        {
          std::fill(current_fluid_force.begin(), current_fluid_force.end(), 0);
          std::fill(current_fluid_torque.begin(),
                    current_fluid_torque.end(),
                    0);
          std::fill(contact_torque.begin(), contact_torque.end(), 0);
          std::fill(contact_force.begin(), contact_force.end(), 0);
          std::fill(contact_wall_force.begin(), contact_wall_force.end(), 0);
          std::fill(contact_wall_torque.begin(), contact_wall_torque.end(), 0);
          std::fill(lubrication_force.begin(), lubrication_force.end(), 0);
          std::fill(lubrication_torque.begin(), lubrication_torque.end(), 0);
          std::fill(lubrication_wall_force.begin(),
                    lubrication_wall_force.end(),
                    0);
          std::fill(lubrication_wall_torque.begin(),
                    lubrication_wall_torque.end(),
                    0);

          // Calculate particle-particle and particle-wall contact force
          calculate_pp_contact_force(dt_dem, contact_force, contact_torque);
          calculate_pw_contact_force(dt_dem,
                                     contact_wall_force,
                                     contact_wall_torque);
          calculate_pp_lubrication_force(
            dt_dem, h_max, h_min, mu, lubrication_force, lubrication_torque);
          calculate_pw_lubrication_force(dt_dem,
                                         h_max,
                                         h_min,
                                         mu,
                                         lubrication_wall_force,
                                         lubrication_wall_torque);
          // define local time of the rk step
          double local_dt = dt_dem * 0.5;
          if (step == 3)
            {
              local_dt = local_dt * 2;
            }

          for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
            {
              if (step == 0)
                {
                  last_velocity[p_i] = dem_particles[p_i].velocity;
                  last_position[p_i] = dem_particles[p_i].position;
                  last_omega[p_i]    = dem_particles[p_i].omega;
                }
              auto inv_inertia = invert(dem_particles[p_i].inertia);
              if (dim == 2)
                {
                  gravity = g * (dem_particles[p_i].mass -
                                 dem_particles[p_i].radius *
                                   dem_particles[p_i].radius * PI * rho);
                }
              else
                {
                  gravity = g * (dem_particles[p_i].mass -
                                 4.0 / 3.0 * dem_particles[p_i].radius *
                                   dem_particles[p_i].radius *
                                   dem_particles[p_i].radius * PI * rho);
                }

              // We consider only the force at t+dt so the scheme is consistent
              // to a BDFn scheme on the fluid side. If there is no contact.

              current_fluid_force[p_i]  = dem_particles[p_i].fluid_forces;
              current_fluid_torque[p_i] = dem_particles[p_i].fluid_torque;

              // Store each rk step of the variable we integrate in its
              // respective vector.

              k_velocity[p_i][step] =
                (current_fluid_force[p_i] + contact_force[p_i] +
                 contact_wall_force[p_i] + gravity + lubrication_force[p_i] +
                 lubrication_wall_force[p_i]) /
                dem_particles[p_i].mass;

              for (unsigned int d = 0; d < dim; ++d)
                {
                  k_position[p_i][step][d] = dem_particles[p_i].velocity[d];
                }

              k_omega[p_i][step] =
                inv_inertia * (current_fluid_torque[p_i] + contact_torque[p_i] +
                               contact_wall_torque[p_i]);

              k_impulsion[p_i][step] =
                (current_fluid_force[p_i] + gravity + contact_wall_force[p_i] +
                 contact_force[p_i] + lubrication_force[p_i] +
                 lubrication_wall_force[p_i]);

              k_omega_impulsion[p_i][step] =
                (current_fluid_torque[p_i] + contact_torque[p_i] +
                 contact_wall_torque[p_i]);

              k_contact_impulsion[p_i][step] =
                (contact_wall_force[p_i] + contact_force[p_i]);

              k_omega_contact_impulsion[p_i][step] =
                (contact_torque[p_i] + contact_wall_torque[p_i]);


              // Integrate the relevant state Variable for the next RK step.
              dem_particles[p_i].velocity =
                last_velocity[p_i] + k_velocity[p_i][step] * local_dt;

              for (unsigned int d = 0; d < dim; ++d)
                {
                  dem_particles[p_i].position[d] =
                    last_position[p_i][d] + k_position[p_i][step][d] * local_dt;
                }

              dem_particles[p_i].omega =
                last_omega[p_i] + k_omega[p_i][step] * local_dt;
            }
        }


      for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
        {
          // Define the integral by combining each of the RK step.

          dem_particles[p_i].velocity =
            last_velocity[p_i] +
            dt_dem *
              (k_velocity[p_i][0] + 2 * k_velocity[p_i][1] +
               2 * k_velocity[p_i][2] + k_velocity[p_i][3]) /
              6;

          for (unsigned int d = 0; d < dim; ++d)
            {
              dem_particles[p_i].position[d] =
                last_position[p_i][d] +
                dt_dem *
                  (k_position[p_i][0][d] + 2 * k_position[p_i][1][d] +
                   2 * k_position[p_i][2][d] + k_position[p_i][3][d]) /
                  6;
            }

          dem_particles[p_i].omega =
            last_velocity[p_i] + dt_dem *
                                   (k_omega[p_i][0] + 2 * k_omega[p_i][1] +
                                    2 * k_omega[p_i][2] + k_omega[p_i][3]) /
                                   6;

          // Integration of the impulsion applied to the particle.
          // This is what will be transferred to the CFD to integrate the
          // particle.
          dem_particles[p_i].impulsion +=
            dt_dem *
            (k_impulsion[p_i][0] + 2 * k_impulsion[p_i][1] +
             2 * k_impulsion[p_i][2] + k_impulsion[p_i][3]) /
            6;

          dem_particles[p_i].contact_impulsion +=
            dt_dem *
            (k_contact_impulsion[p_i][0] + 2 * k_contact_impulsion[p_i][1] +
             2 * k_contact_impulsion[p_i][2] + k_contact_impulsion[p_i][3]) /
            6;

          dem_particles[p_i].omega_impulsion +=
            dt_dem *
            (k_omega_impulsion[p_i][0] + 2 * k_omega_impulsion[p_i][1] +
             2 * k_omega_impulsion[p_i][2] + k_omega_impulsion[p_i][3]) /
            6;

          dem_particles[p_i].omega_contact_impulsion +=
            dt_dem *
            (k_omega_contact_impulsion[p_i][0] +
             2 * k_omega_contact_impulsion[p_i][1] +
             2 * k_omega_contact_impulsion[p_i][2] +
             k_omega_contact_impulsion[p_i][3]) /
            6;
        }

      t += dt_dem;
    }
}


template class IBParticlesDEM<2>;
template class IBParticlesDEM<3>;
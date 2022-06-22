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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef LETHE_GLSSHARPNS_H
#define LETHE_GLSSHARPNS_H

#include <core/ib_particle.h>
#include <core/ib_stencil.h>

#include <solvers/gls_navier_stokes.h>

#include <fem-dem/ib_particles_dem.h>

#include <deal.II/dofs/dof_tools.h>

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization and
 * Sharp-Edge immersed boundaries
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Lucka Barbeau, Bruno Blais, 2020
 */

template <int dim>
class GLSSharpNavierStokesSolver : public GLSNavierStokesSolver<dim>
{
public:
  GLSSharpNavierStokesSolver(SimulationParameters<dim> &nsparam);

  ~GLSSharpNavierStokesSolver();

  void
  solve();

  /**
   * @brief Call for the assembly of the matrix
   */
  void
  assemble_system_matrix()
  {
    assemble_matrix_and_rhs();
  }

  /**
   * @brief Call for the assembly of the right-hand side
   */
  void
  assemble_system_rhs()
  {
    assemble_rhs();
    sharp_edge();
  }


private:
  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data) override;

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data) override;

  /**
   * @brief sets up the vector of assembler functions
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle. Refer to 2 different assembler
   * depending on what type of equations is assembled inside the particles.
   */
  void
  setup_assemblers() override;


  /**
   * @brief Copy local cell information to global matrix
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copy local cell rhs information to global rhs
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Call for the assembly of the matrix and the right hand side
   *
   * @deprecated This function is to be deprecated when the non-linear solvers
   * have been refactored to call for rhs and matrix assembly seperately.
   */

  /*
   Modified version of assemble_matrix_and_rhs to include the presence of
   extra steps. For more detail see the same function in the
   gls_navier_stokes.h solver.
   */

  virtual void
  assemble_matrix_and_rhs()
  {
    if (this->simulation_parameters.particlesParameters->integrate_motion)
      {
        force_on_ib();
        integrate_particles();
        generate_cut_cells_map();
      }
    // this->simulation_control->set_assembly_method(this->time_stepping_method);
    {
      this->GLSNavierStokesSolver<
        dim>::assemble_system_matrix_without_preconditioner();
    }

    sharp_edge();

    // Assemble the preconditioner
    this->setup_preconditioner();
  }


  /**
   * @brief Call for the assembly of the right hand side
   *
   * @deprecated This function is to be deprecated when the non-linear solvers
   * have been refactored to call for rhs and matrix assembly seperately.
   *
   * Modified version of assemble_matrix_and_rhs to include the presence of
   * extra steps. For more detail see the same function in the
   * gls_navier_stokes.h solver.
   */
  virtual void
  assemble_rhs()
  {
    this->GLSNavierStokesSolver<dim>::assemble_system_rhs();
  }

  /**
   * @brief
   * Map the vertex index to the cell that includes that vertex.
   * This map is used to find all the cell close to a specific vertex.
   */
  void
  vertices_cell_mapping();

  /**
   * @brief
   * Defines the particle structure and value based on the parameter file.
   * This structure gives access to the position, velocity, force and other
   * proprieties of each IB particle. All the variables defined for each of the
   * particle are described in the class: IBparticle. see:
   * \include\core\ib_particle.h
   */
  void
  define_particles();

  /**
   * @brief
   * Evaluate the forces applied on each of the IB particle.
   */
  void
  force_on_ib();


  /**
   * @brief
   * Modify the system matrix to impose IB condition using the sharp_edge
   * approach. The detail of this approach are presented in : L. Barbeau, S.
   * Étienne, C. Béguin & B. Blais, «Development of a high-order continuous
   * Galerkin sharp-interface immersed boundary method and its application to
   * incompressible flow problems,» Computers & Fluids, 2020, in press, ref.
   * CAF-D-20-00773
   */
  void
  sharp_edge();

  /**
   * @brief
   * Write in a specifique file for each of the paticles its forces, velocity,
   * position at each time step.
   */
  void
  write_force_ib();

  /**
   * @brief
   * Integrate the particle velocity and position based on the forces and
   * torques and applies the next value to the particle.
   */
  void
  integrate_particles();

  /**
   * @brief
   * Store the solution of the particles dynamics parameters for integration.
   * Defines the table to store the history of each of the particles.
   */
  void
  finish_time_step_particles();

  /**
   * @brief
   * Evaluate the L2 error on the computational domain if an analytical solution
   * is given. This function is slightly different from its standard GLS
   * counterpart as the cells that are cut by an immersed boundary should not be
   * taken into account in the error evaluation. See "computation domain"
   * definition in: L. Barbeau, S. Étienne, C. Béguin & B. Blais, «Development
   * of a high-order continuous Galerkin sharp-interface immersed boundary
   * method and its application to incompressible flow problems,» Computers &
   * Fluids, 2020, in press, ref. CAF-D-20-0077
   */
  std::pair<double, double>
  calculate_L2_error_particles();


  /**
   * @brief
   * Same function as its standard GLS counterpart but it used the error
   * evaluation that takes into account the particle’s position.
   */
  virtual void
  postprocess_fd(bool firstIter) override;

  /**
   * @brief
   * Allow a refinement around each of the particles.
   * The zone where the cells will be refined is defined by a ring in 2D and a
   *shell in 3D. The outside and inside radius of the ring\shell is defined in
   *relation to the diameter of the particle by the immersed boundaries
   *parameter: "refine mesh inside radius factor" and "refine mesh outside
   *radius factor". These factors multiply the radius of the particle to define
   *the outside and inside radius of the ring\shell.
   */
  void
  refine_ib();

  /**
   * @brief
   *This function create a map (cut_cells_map) that indicates if a cell is cut,
   *and the particle id of the particle that cut it.
   */
  void
  generate_cut_cells_map();


  /**
   * @brief
   * Return a bool to define if a cell is cut by an IB particle and the local
   * DOFs of the cell for later us. If the cell is cut, the function will return
   * the id of the particle that cut it, else it returns 0.
   *
   * @param cell , the cell that we verify whether it is cut or not.
   *
   * @param local_dof_indices, a container for the local dof indices used during the function call.
   *
   * @param support_points, a mapping of support points for the DOFs.
   *
   */
  std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
  cell_cut(const typename DoFHandler<dim>::active_cell_iterator &cell,
           std::vector<types::global_dof_index> &         local_dof_indices,
           std::map<types::global_dof_index, Point<dim>> &support_points);
  bool
  cell_cut_by_p(std::vector<types::global_dof_index> &local_dof_indices,
                std::map<types::global_dof_index, Point<dim>> &support_points,
                unsigned int                                   p);
  /**
   * @brief
   * Return a bool to define if a cell is contains inside an IB particle and the
   * local DOFs of the cell for later us. If the cell is cut, the function will
   * return the id of the particle, else it returns 0.
   *
   * @param cell , the cell that we verify whether it is cut or not.
   *
   * @param local_dof_indices, a container for the local dof indices used during the function call.
   *
   * @param support_points, a mapping of support points for the DOFs.
   *
   */
  std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
  cell_inside(const typename DoFHandler<dim>::active_cell_iterator &cell,
              std::vector<types::global_dof_index> &         local_dof_indices,
              std::map<types::global_dof_index, Point<dim>> &support_points);



  /**
   * @brief Write a gls_sharp simulation checkpointing to allow for gls_sharp simulation restart
   * This function stores all the previous' particles states in one file. Each
   * row corresponds to one particle state. The file is structured as follows:
   *
   *
   * P0 state at time t
   *
   * P0 state at time t-dt
   *
   * P0 state at time t-2dt
   *
   * P1 state at time t
   *
   * P1 state at time t-dt
   *
   * P1 state at time t-2dt
   *
   * etc
   *
   */

  virtual void
  write_checkpoint() override;
  /*
   * @brief Read a gls_sharp simulation checkpoint and initiate simulation restart.
   * See the description of the function write_checkpoint for more details about
   * the file structure.
   * */
  virtual void
  read_checkpoint() override;

  /**
* @brief
Return a bool that describes  if a cell contains a specific point
*
* @param cell , The initial cell for which we want to check if the point is inside.
*
* @param point, The point that we wish to check
*/
  bool
  point_inside_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                    Point<dim> point);

  /**
   * @brief
   *  A override of the get_current_residual to take into account the particles
   * coupling residual.
   */
  double
  get_current_residual() override
  {
    double scalling = this->simulation_parameters.non_linear_solver.tolerance /
                      this->simulation_parameters.particlesParameters
                        ->particle_nonlinear_tolerance;
    return std::max(this->system_rhs.l2_norm(), particle_residual * scalling);
  }



  /**
   * @brief Defines a struct with methods that allow the generation of a visualisation of the IB_particles. This is equivalent to the corresponding class in the DEM solver.
   */
  struct Visualization_IB : public dealii::DataOutInterface<0, dim>
  {
  public:
    /**
     * Carries out building the patches of properties of particles for
     * visualization
     *
     * @param particles The vector fo IB_particles
     */
    void
    build_patches(std::vector<IBParticle<dim>> particles);


    ~Visualization_IB();

  private:
    /**
     * Implementation of the corresponding function of the base class.
     */
    virtual const std::vector<DataOutBase::Patch<0, dim>> &
    get_patches() const;

    /**
     * Implementation of the corresponding function of the base class.
     */
    virtual std::vector<std::string>
    get_dataset_names() const;

    virtual std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
    get_nonscalar_data_ranges() const;


    /**
     * Output information that is filled by build_patches() and
     * written by the write function of the base class.
     */
    std::vector<DataOutBase::Patch<0, dim>> patches;

    /**
     * A list of field names for all data components stored in patches.
     */
    std::vector<std::string> dataset_names;

    /**
     * Store which of the data fields are vectors.
     */

    std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      vector_datasets;


    /**
     * Particle properties that are written in output files
     */
    std::vector<std::pair<std::string, int>> properties_to_write;
  };



  /**
   * Members
   */

private:
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;
  /*
   * This map uses the cell as the key, and store the following information:
   * if that cell is cut (bool) and what particle cut this cell (unsigned int).
   */
  std::map<typename DoFHandler<dim>::active_cell_iterator,
           std::tuple<bool, unsigned int>>
    cut_cells_map;

  std::map<typename DoFHandler<dim>::active_cell_iterator,
           std::tuple<bool, unsigned int>>
    cells_inside_map;
  /*
   * This map is used to keep in memory which DOFs already have an IB equation
   * imposed on them in order to avoid writing multiple time the same equation.
   */
  std::map<unsigned int,
           std::pair<bool, typename DoFHandler<dim>::active_cell_iterator>>
    ib_done;

  // Special assembler of the cells inside an IB particle
  std::vector<std::shared_ptr<NavierStokesAssemblerBase<dim>>>
    assemblers_inside_ib;

  PVDHandler ib_particles_pvdhandler;


  std::vector<IBParticle<dim>> particles;
  double                       particle_residual;

  std::vector<TableHandler> table_p;

  // Object used to sub-time step the particle dynamics to allow contact between
  // particles.
  IBParticlesDEM<dim> ib_dem;
};


#endif
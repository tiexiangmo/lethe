
// Lethe includes
#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/periodic_hills_grid.h>

// Deal.II includes
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

// Std
#include <fstream>
#include <iostream>


template <int dim, int spacedim>
void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
  const Parameters::Mesh &                               mesh_parameters)

{
  // GMSH input
  if (mesh_parameters.type == Parameters::Mesh::Type::gmsh)
    {
      if (mesh_parameters.simplex)
        {
          auto        comm      = triangulation.get_communicator();
          std::string file_name = mesh_parameters.file_name;

          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation_in_groups<dim, spacedim>(
              [file_name](dealii::Triangulation<dim, spacedim> &basetria) {
                GridIn<dim, spacedim> grid_in;
                grid_in.attach_triangulation(basetria);
                std::ifstream input_file(file_name);

                grid_in.read_msh(input_file);
              },
              [](dealii::Triangulation<dim, spacedim> &basetria,
                 const MPI_Comm                        comm,
                 const unsigned int /*group_size*/) {
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm), basetria);
              },
              comm,
              Utilities::MPI::n_mpi_processes(comm) /* group size */,
              dealii::Triangulation<dim, spacedim>::none);

          triangulation.create_triangulation(construction_data);
        }
      else
        {
          GridIn<dim, spacedim> grid_in;
          grid_in.attach_triangulation(triangulation);
          std::ifstream input_file(mesh_parameters.file_name);
          grid_in.read_msh(input_file);
        }
    }
  // Dealii grids
  else if (mesh_parameters.type == Parameters::Mesh::Type::dealii)
    {
      if (mesh_parameters.simplex)
        {
          Triangulation<dim, spacedim> temporary_quad_triangulation;
          GridGenerator::generate_from_name_and_arguments(
            temporary_quad_triangulation,
            mesh_parameters.grid_type,
            mesh_parameters.grid_arguments);

          // initial refinement
          const unsigned int initial_refinement =
            mesh_parameters.initial_refinement;
          temporary_quad_triangulation.refine_global(initial_refinement);

          // flatten the triangulation
          Triangulation<dim, spacedim> flat_temp_quad_triangulation;
          GridGenerator::flatten_triangulation(temporary_quad_triangulation,
                                               flat_temp_quad_triangulation);

          Triangulation<dim, spacedim> temporary_tri_triangulation(
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
          GridGenerator::convert_hypercube_to_simplex_mesh(
            flat_temp_quad_triangulation, temporary_tri_triangulation);

          GridTools::partition_triangulation_zorder(
            Utilities::MPI::n_mpi_processes(triangulation.get_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              triangulation.get_communicator(),
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          triangulation.create_triangulation(construction_data);
        }
      else
        {
          GridGenerator::generate_from_name_and_arguments(
            triangulation,
            mesh_parameters.grid_type,
            mesh_parameters.grid_arguments);
        }
    }
  // Customizable cylinder mesh
  else if (mesh_parameters.type == Parameters::Mesh::Type::cylinder)
    {
      if (mesh_parameters.simplex)
        {
          throw std::runtime_error(
            "Unsupported mesh type - custom cylinder mesh with simplex is not supported. Use a dealii cylinder to use simplex mesh.");
        }
      else if (dim != 3)
        {
          throw std::runtime_error(
            "Unsupported mesh type - custom cylinder mesh is only supported in 3d.");
        }

      if constexpr (dim == 3)
        {
          // Separate arguments of the string
          std::vector<std::string> arguments;
          std::stringstream        s_stream(mesh_parameters.grid_arguments);
          while (s_stream.good())
            {
              std::string substr;
              getline(s_stream, substr, ':');
              arguments.push_back(substr);
            }

          // Arguments declaration
          unsigned int subdivisions;
          double       radius, half_height;
          if (arguments.size() != 3)
            {
              throw std::runtime_error(
                "Mandatory cylinder parameters are (x subdivisions: radius : half height)");
            }
          else
            {
              std::vector<double> arguments_double =
                dealii::Utilities::string_to_double(arguments);
              subdivisions = static_cast<int>(arguments_double[0]);
              radius       = arguments_double[1];
              half_height  = arguments_double[2];
            }

          if (mesh_parameters.grid_type == "classic")
            {
              // Create a subdivided cylinder from deal.ii
              GridGenerator::subdivided_cylinder(triangulation,
                                                 subdivisions,
                                                 radius,
                                                 half_height);
            }
          else
            {
              // Create a temporary 2d mesh
              Triangulation<2, spacedim - 1> temporary_triangulation;

              // Create a spherical manifold for 2d mesh
              Point<2>                                 center(0.0, 0.0);
              const SphericalManifold<2, spacedim - 1> m0(center);

              if (mesh_parameters.grid_type == "regularized" ||
                  mesh_parameters.grid_type == "squared")
                {
                  // Create a square mesh
                  double real_radius = radius * std::sin(M_PI_4);
                  GridGenerator::hyper_cube(temporary_triangulation,
                                            -real_radius,
                                            real_radius,
                                            true);

                  // Assign boundary 0 to perimeter as for cylinder
                  for (const auto &cell :
                       temporary_triangulation.active_cell_iterators())
                    {
                      if (cell->is_locally_owned())
                        {
                          // Looping through all the faces of the cell
                          for (const auto &face : cell->face_iterators())
                            {
                              // Check to see if the face is located at boundary
                              if (face->at_boundary())
                                {
                                  face->set_boundary_id(0);
                                }
                            }
                        }
                    }
                }
              else if (mesh_parameters.grid_type == "balanced")
                {
                  GridGenerator::hyper_ball_balanced(temporary_triangulation,
                                                     center,
                                                     radius);
                }
              else
                {
                  throw std::runtime_error(
                    "Unknown grid type. Choices are <classic|balanced|squared|regularized>.");
                }

              temporary_triangulation.reset_all_manifolds();
              temporary_triangulation.set_all_manifold_ids_on_boundary(0);
              temporary_triangulation.set_manifold(0, m0);

              if (mesh_parameters.grid_type == "regularized")
                {
                  // Pre-refinement to reduce mesh size at corners before
                  // regularization
                  temporary_triangulation.refine_global(2);
                  GridTools::regularize_corner_cells(temporary_triangulation);

                  // Flatten the triangulation
                  Triangulation<2, spacedim - 1> flat_temporary_triangulation;
                  flat_temporary_triangulation.copy_triangulation(
                    temporary_triangulation);
                  temporary_triangulation.clear();
                  GridGenerator::flatten_triangulation(
                    flat_temporary_triangulation, temporary_triangulation);
                }

              // Extrude the 2d temporary mesh to 3d cylinder
              GridGenerator::extrude_triangulation(temporary_triangulation,
                                                   subdivisions + 1,
                                                   2.0 * half_height,
                                                   triangulation,
                                                   true);

              // Rotate mesh in x-axis and set the (0,0,0) at the barycenter
              // to be comparable to dealii cylinder meshes
              Tensor<1, spacedim> axis_vector({0.0, 1.0, 0.0});
              GridTools::rotate(axis_vector, M_PI_2, triangulation);
              Tensor<1, spacedim> shift_vector({-half_height, 0.0, 0.0});
              GridTools::shift(shift_vector, triangulation);

              // Add a cylindrical manifold on the final unrefined mesh
              const CylindricalManifold<3, spacedim> m1(0);
              triangulation.set_manifold(0, m1);
            }
        }
    }
  // Periodic Hills grid
  else if (mesh_parameters.type == Parameters::Mesh::Type::periodic_hills &&
           !mesh_parameters.simplex)
    {
      if (mesh_parameters.simplex)
        {
          throw std::runtime_error(
            "Unsupported mesh type - periodic hills mesh with simplex is not supported");
        }
      else
        {
          PeriodicHillsGrid<dim, spacedim> grid(mesh_parameters.grid_arguments);
          grid.make_grid(triangulation);
        }
    }

  // Radial injector grid
  else if (mesh_parameters.type == Parameters::Mesh::Type::radial_injector)
    {
      const double c1_length = 0.05172;
      const double c2_length = 0.09+0.0337;
      const double c1_radius = 0.131525;
      const double s1_radius = 0.00635 + c1_radius;
      const double s2_radius = 0.00635 + s1_radius;
      const double s3_length = 0.0822;
      const double s4_length = 0.0337;

      const unsigned int c1_division = 5;
      const unsigned int c2_division = 5;
      const unsigned int s3_division = 5;
      const unsigned int s4_division = 5;



      // Create cylinders c1 and c2 then merge them
      Triangulation<3> cylinder_1;
      GridGenerator::subdivided_cylinder(cylinder_1,
                                         c1_division,
                                         c1_radius,
                                         c1_length / 2.);
      GridTools::shift(Tensor<1, 3>({c1_length / 2, 0., 0}), cylinder_1);

      // Shift cylinder 1
      Triangulation<3> cylinder_2;
      GridGenerator::subdivided_cylinder(cylinder_2,
                                         c2_division,
                                         c1_radius,
                                         c2_length / 2.);
      GridTools::shift(Tensor<1, 3>({c1_length + c2_length / 2., 0., 0}),
                       cylinder_2);

      Triangulation<3> merged_cylinder;
      GridGenerator::merge_triangulations(
        cylinder_1, cylinder_2, merged_cylinder, 1e-12, true, true);
      GridTools::rotate(Tensor<1, 3>({0., 1., 0.}), -M_PI / 2, merged_cylinder);
      GridTools::rotate(Tensor<1, 3>({0., 0, 1.}), -M_PI / 4, merged_cylinder);



      // Create cylinder shells and merge them
      Triangulation<3> shell_1;
      GridGenerator::cylinder_shell(
        shell_1, c1_length, c1_radius, s1_radius, 4, c1_division);

      Triangulation<3> shell_2;
      GridGenerator::cylinder_shell(
        shell_2, c1_length, s1_radius, s2_radius, 4, c1_division);

      Triangulation<3> shell_3;
      GridGenerator::cylinder_shell(
        shell_3, s3_length, s1_radius, s2_radius, 4, s3_division);
      GridTools::shift(Tensor<1, 3>({0, 0., c1_length}), shell_3);

      Triangulation<3> shell_4;
      GridGenerator::cylinder_shell(
        shell_4, s4_length, s1_radius, s2_radius, 4, s4_division);
      GridTools::shift(Tensor<1, 3>({0, 0., c1_length + s3_length}), shell_4);

      Triangulation<3> s1_s2;
      GridGenerator::merge_triangulations(
        shell_1, shell_2, s1_s2, 1e-12, true, false);

      Triangulation<3> s1_s2_s3;
      GridGenerator::merge_triangulations(
        shell_3, s1_s2, s1_s2_s3, 1e-12, true, false);

      Triangulation<3> s1_s2_s3_s4;
      GridGenerator::merge_triangulations(
        shell_4, s1_s2_s3, s1_s2_s3_s4, 1e-12, true, false);

      Triangulation<3> final_merge;
      GridGenerator::merge_triangulations(
        s1_s2_s3_s4, merged_cylinder, final_merge, 1e-12, true, false);


      // Copy triangulation
      if constexpr (dim == 3)
        {
          final_merge.set_manifold(0, CylindricalManifold<3>(2)); // axis z


          Triangulation<3>::cell_iterator cell = final_merge.begin();

          for (; cell != final_merge.end(); ++cell)
            {
              // if (cell->at_boundary())
              {
                for (const unsigned int f : GeometryInfo<3>::face_indices())
                  {
                    if (cell->face(f)->at_boundary())
                      {
                        Point<3>     center = cell->face(f)->center();
                        const double radius = std::sqrt(center[0] * center[0] +
                                                        center[1] * center[1]);
                        const double z      = center[2];
                        if (z > (c2_length + c1_length - 1e-3))
                          {
                            if (radius < (c1_radius-1e-3))
                              {
                                std::cout << "radius: " << radius << " center radius: " << c1_radius << std::endl;
                                cell->face(f)->set_boundary_id(1);
                              }
                              else
                              {
                                cell->face(f)->set_boundary_id(0);
                                std::cout << "radius: " << radius << " center radius: " << c1_radius << std::endl;
                              }
                          }
                        for (unsigned int j = 0;
                             j < GeometryInfo<3>::lines_per_face;
                             ++j)
                          {
                            Point<3> center = cell->face(f)->line(j)->center();
                            const double radius = std::sqrt(
                              center[0] * center[0] + center[1] * center[1]);
                            const double z = center[2];
                            if (z <
                                ((s3_length + s4_length + c1_length) - 1e-10))
                              {
                                if (z > ((s3_length + c1_length) - 1e-10))
                                  {
                                    if (radius > (s1_radius + 1e-3))
                                      {
                                        cell->face(f)->set_boundary_id(2);
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
            }
          triangulation.copy_triangulation(final_merge);
        }

      // Manually set boundary conditions



      // GridGenerator::subdivided_cylinder(cylinder_2,
      //                                    c2_division,
      //                                    c1_radius,
      //                                    c2_length / 2.);
      // Triangulation<3> shell_2;
      // Triangulation<3> shell_3;



      // Create bottom cylinder shell
      //      Triangulation<3> cylinder_shell_bottom;
      //      GridGenerator::cylinder_shell(cylinder_shell_bottom,
      //                                    length * (1 - depth_ratio),
      //                                    inner_radius,
      //                                    outer_radius,
      //                                    4,
      //                                    half_n_bottom_cells);

      // Create top cylinder shell
      // Triangulation<3> cylinder_shell_top;
      //      GridGenerator::cylinder_shell(cylinder_shell_top,
      //                                    length * depth_ratio,
      //                                    inner_radius,
      //                                    outer_radius,
      //                                    4,
      //                                    half_n_top_cells);
      // Shift top cylinder shell to fit on top of the
      // first
      //      GridTools::shift(Tensor<1, 3>({0., 0., length * (1 -
      //      depth_ratio)}),
      //                       cylinder_shell_top);
      // Create full cylinder to fill bottom cylinder shell
      //      double           full_cyl_half_length = (length * (1 -
      //      depth_ratio)) / 2.;
      //      GridGenerator::subdivided_cylinder(full_cylinder,
      //                                         half_n_bottom_cells,
      //                                         inner_radius,
      //                                         full_cyl_half_length);
      // Rotate cylinder such that x->z and z->-x
      // GridTools::rotate(Tensor<1, 3>({0., 1., 0.}), -M_PI / 2,
      // full_cylinder);
      // Shift full cylinder so that the bottom part is at 0
      //      GridTools::shift(Tensor<1, 3>({0., 0., full_cyl_half_length}),
      //                       full_cylinder);
      // Create temp triangulation to merge bottom cylinder
      // shell and full cylinder

      // GridTools::rotate(Tensor<1, 3>({0., 0., 1.}), -M_PI / 4,
      // full_cylinder);


      // Triangulation<3> temp;
      // GridGenerator::merge_triangulations(
      //   full_cylinder, cylinder_shell_bottom, temp, 1e-12, true, false);
      //  Merge bottom part to top cylinder shell
      //       GridGenerator::merge_triangulations(
      //         temp, cylinder_shell_top, triangulation, 1e-12,
      //         true, false);
      //  Keep previously copied manifolds and set Cylindrical
      //  manifold along z-axis
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");
}


template <int dim, int spacedim>
void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<dim, spacedim> & triangulation,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)

{
  // Setup periodic boundary conditions
  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      if (boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          auto triangulation_ptr = &triangulation;

          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim, spacedim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(
            *dynamic_cast<Triangulation<dim, spacedim> *>(triangulation_ptr),
            boundary_conditions.id[i_bc],
            boundary_conditions.periodic_id[i_bc],
            boundary_conditions.periodic_direction[i_bc],
            periodicity_vector);
          triangulation.add_periodicity(periodicity_vector);
        }
    }
}

template <int dim, int spacedim>
void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<dim, spacedim> & triangulation,
  const Parameters::Mesh &                                mesh_parameters,
  const Parameters::Manifolds &                           manifolds_parameters,
  const bool &                                            restart,
  const BoundaryConditions::BoundaryConditions<spacedim> &boundary_conditions)
{
  attach_grid_to_triangulation(triangulation, mesh_parameters);

  setup_periodic_boundary_conditions(triangulation, boundary_conditions);

  attach_manifolds_to_triangulation(triangulation, manifolds_parameters);

  if (mesh_parameters.simplex)
    {
      // Refinement isn't possible yet
    }
  else
    {
      if (mesh_parameters.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size = mesh_parameters.target_size;
          unsigned int number_refinement =
            floor(std::log(minimal_cell_size / target_size) / std::log(2));
          triangulation.refine_global(number_refinement);
        }
      else if (!restart)
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }
}

template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<2> &triangulation,
  const Parameters::Mesh &                   mesh_parameters);
template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<3> &triangulation,
  const Parameters::Mesh &                   mesh_parameters);
template void
attach_grid_to_triangulation(
  parallel::DistributedTriangulationBase<2, 3> &triangulation,
  const Parameters::Mesh &                      mesh_parameters);


template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 2> &   triangulation,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<2, 3> &   triangulation,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void
setup_periodic_boundary_conditions(
  parallel::DistributedTriangulationBase<3, 3> &   triangulation,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);

template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2> &      triangulation,
  const Parameters::Mesh &                         mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<2> &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<3> &      triangulation,
  const Parameters::Mesh &                         mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);
template void
read_mesh_and_manifolds(
  parallel::DistributedTriangulationBase<2, 3> &   triangulation,
  const Parameters::Mesh &                         mesh_parameters,
  const Parameters::Manifolds &                    manifolds_parameters,
  const bool &                                     restart,
  const BoundaryConditions::BoundaryConditions<3> &boundary_conditions);

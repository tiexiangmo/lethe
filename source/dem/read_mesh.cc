#include <dem/read_mesh.h>

template <int dim, int spacedim>
void
read_mesh(const Parameters::Mesh &      mesh_params,
          const bool &                  restart,
          const ConditionalOStream &    pcout,
          Triangulation<dim, spacedim> &triangulation,
          double &                      triangulation_cell_diameter)
{
  pcout << "Reading triangulation " << std::endl;
  // GMSH input
  if (mesh_params.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(mesh_params.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (mesh_params.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        triangulation, mesh_params.grid_type, mesh_params.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  triangulation_cell_diameter = 0.5 * GridTools::diameter(triangulation);

  if (!restart)
    {
      if (mesh_params.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size = mesh_params.target_size;
          unsigned int number_refinement =
            floor(std::log(minimal_cell_size / target_size) / std::log(2));
          pcout << "Automatically refining grid until target size : "
                << target_size << std::endl;
          triangulation.refine_global(number_refinement);
          pcout << "Mesh was automatically refined : " << number_refinement
                << " times" << std::endl;
        }
      else
        {
          const int initial_refinement = mesh_params.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }

  pcout << std::endl << "Finished reading triangulation " << std::endl;
}

template void
read_mesh<1, 2>(const Parameters::Mesh &  mesh_params,
                const bool &              restart,
                const ConditionalOStream &pcout,
                Triangulation<1, 2> &     triangulation,
                double &                  triangulation_cell_diameter);

template void
read_mesh<2, 2>(const Parameters::Mesh &  mesh_params,
                const bool &              restart,
                const ConditionalOStream &pcout,
                Triangulation<2, 2> &     triangulation,
                double &                  triangulation_cell_diameter);

template void
read_mesh<2, 3>(const Parameters::Mesh &  mesh_params,
                const bool &              restart,
                const ConditionalOStream &pcout,
                Triangulation<2, 3> &     triangulation,
                double &                  triangulation_cell_diameter);

template void
read_mesh<3, 3>(const Parameters::Mesh &  mesh_params,
                const bool &              restart,
                const ConditionalOStream &pcout,
                Triangulation<3, 3> &     triangulation,
                double &                  triangulation_cell_diameter);
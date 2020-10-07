#include <dem/update_ghost_iterator_pp_contact_container.h>

using namespace dealii;

template <int dim>
void
update_ghost_iterator_pp_contact_container(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &ghost_particle_container)
{
  for (auto adjacent_particles_iterator = ghost_adjacent_particles.begin();
       adjacent_particles_iterator != ghost_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();
           ++pp_map_iterator)
        {
          int particle_two_id = pp_map_iterator->first;
          pp_map_iterator->second.particle_two =
            ghost_particle_container.at(particle_two_id);
        }
    }
}

template void
update_ghost_iterator_pp_contact_container(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &ghost_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<2>>
    &ghost_particle_container);

template void
update_ghost_iterator_pp_contact_container(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<3>>>
    &ghost_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<3>>
    &ghost_particle_container);
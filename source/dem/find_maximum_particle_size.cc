#include <dem/find_maximum_particle_size.h>

double
find_maximum_particle_size_for_normal_distribution(
  const Parameters::Lagrangian::LagrangianPhysicalProperties
              &physical_properties,
  const double standard_deviation_multiplier)
{
  double maximum_particle_size = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      maximum_particle_size =
        std::max(maximum_particle_size,
                 physical_properties.particle_average_diameter.at(i) +
                   standard_deviation_multiplier *
                     physical_properties.particle_size_std.at(i));
    }
  return maximum_particle_size;
}


double
find_maximum_particle_size_for_lognormal_distribution(
  const Parameters::Lagrangian::LagrangianPhysicalProperties
              &physical_properties,
  const double standard_deviation_multiplier)
{
  double maximum_particle_size = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      maximum_particle_size =
        std::max(maximum_particle_size,
                 physical_properties.particle_average_diameter.at(i) +
                   standard_deviation_multiplier *
                     physical_properties.particle_size_std.at(i));
    }
  return maximum_particle_size;
}

double
find_maximum_particle_size_for_histogram_distribution(
  const Parameters::Lagrangian::LagrangianPhysicalProperties
    &physical_properties)
{
  double maximum_particle_size = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      for (double j = 0;
           j < physical_properties.particle_list_diameter.at(i).size();
           j++)
        {
          maximum_particle_size =
            std::max(maximum_particle_size,
                     physical_properties.particle_list_diameter.at(i).at(j));
        }
    }

  return maximum_particle_size;
}

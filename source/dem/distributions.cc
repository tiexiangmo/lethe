#include <dem/distributions.h>

// Normal
NormalDistribution::NormalDistribution(
  const std::unordered_map<unsigned int, double> &d_averages,
  const std::unordered_map<unsigned int, double> &d_standard_deviations)
{
  diameter_averages   = d_averages;
  standard_deviations = d_standard_deviations;
}

void
NormalDistribution::particle_size_sampling(const unsigned int &particle_number,
                                           const unsigned int &particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{
    diameter_averages.at(particle_type), standard_deviations.at(particle_type)};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(distribution(gen));
}

// Uniform
UniformDistribution::UniformDistribution(
  const std::unordered_map<unsigned int, double> &d_values)
{
  diameter_values = d_values;
}

void
UniformDistribution::particle_size_sampling(const unsigned int &particle_number,
                                            const unsigned int &particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back(this->diameter_values.at(particle_type));
}

LogNormalDistribution::LogNormalDistribution(
  const std::unordered_map<unsigned int, double> &d_averages,
  const std::unordered_map<unsigned int, double> &d_standard_deviations)
{
  diameter_averages   = d_averages;
  standard_deviations = d_standard_deviations;
}

// LogNormal
void
LogNormalDistribution::particle_size_sampling(
  const unsigned int &particle_number,
  const unsigned int &particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  double mu    = diameter_averages.at(particle_type);
  double sigma = standard_deviations.at(particle_type);
  double log_standard_deviation =
    std::sqrt(std::log(1 + (sigma * sigma) / (mu * mu)));

  std::random_device            rd{};
  std::mt19937                  gen{rd()};
  std::lognormal_distribution<> distribution{std::log(mu),
                                             log_standard_deviation};

  for (unsigned int n = 0; n < particle_number; ++n)
    this->particle_sizes.push_back((distribution(gen)));
}

// Histogram
HistogramDistribution::HistogramDistribution(
  const std::unordered_map<unsigned int, std::vector<double>> &d_list,
  const std::unordered_map<unsigned int, std::vector<double>> &d_probabilities)
{
  diameter_list = d_list;

  for (const auto &it : d_probabilities)
    {
      std::vector<double> cummulative_probability_vector;
      cummulative_probability_vector.reserve(it.second.size());

      double cumulative_value;
      for (const auto &prob : it.second)
        {
          cumulative_value += prob;
          cummulative_probability_vector.push_back(cumulative_value);
        }

      diameter_cumulative_probability.insert(
        {it.first, cummulative_probability_vector});
    }
}

void
HistogramDistribution::particle_size_sampling(
  const unsigned int &particle_number,
  const unsigned int &particle_type)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(particle_number);

  std::random_device               rd{};
  std::mt19937                     gen{rd()};
  std::uniform_real_distribution<> dis(0.0, 1.0);


  for (unsigned int i = 0; i < particle_number; ++i)
    {
      std::cout << __LINE__ << std::endl;

      // Search to find the appropriate diameter index
      auto it = std::upper_bound(
        diameter_cumulative_probability.at(particle_type).begin(),
        diameter_cumulative_probability.at(particle_type).end(),
        dis(gen));

      unsigned int index =
        std::distance(diameter_cumulative_probability.at(particle_type).begin(),
                      it);
      std::cout << diameter_cumulative_probability.at(particle_type).at(index)
                << std::endl;
      this->particle_sizes.push_back(
        diameter_cumulative_probability.at(particle_type).at(index));
    }
}

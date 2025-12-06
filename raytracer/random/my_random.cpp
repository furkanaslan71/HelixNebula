#include "my_random.h"
#include "../external/gsl/gsl"

float generateRandomFloat(float start, float end)
{
  static thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<float> distribution(start, end);
  return distribution(generator);
}

std::vector<float> generateNRandomFloats(float start, float end, float N)
{
  std::vector<float> result(N);
  for (int i = 0; i < N; i++)
  {
    result.push_back(generateRandomFloat(start, end));
  }
  return result;
}

std::vector<std::pair<float, float>> generateJitteredSamples(int num_samples)
{
  Expects(num_samples > 0);
  int n = std::sqrt(num_samples);
  Expects(n * n == num_samples);
  std::vector<std::pair<float, float>> samples(num_samples);
  int i = 0;
  for (int y = 0; y < n; y ++)
  {
    for (int x = 0; x < n; x++)
    {
      float psi1 = generateRandomFloat(0, 1);
      float psi2 = generateRandomFloat(0, 1);
      samples[i].first = (x + psi1) / n;
      samples[i].second = (y + psi2) / n;
      i++;
    }
  }
  return samples;
}
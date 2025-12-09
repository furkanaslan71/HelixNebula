#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>

inline float generateRandomFloat(float start, float end)
{
  static thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<float> distribution(start, end);
  return distribution(generator);
}

std::vector<float> generateNRandomFloats(float start, float end, float N);

std::vector<std::pair<float, float>> generateJitteredSamples(int num_samples);
#endif // !RANDOM_H

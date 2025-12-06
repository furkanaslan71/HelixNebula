#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>

float generateRandomFloat(float start, float end);

std::vector<float> generateNRandomFloats(float start, float end, float N);

std::vector<std::pair<float, float>> generateJitteredSamples(int num_samples);
#endif // !RANDOM_H

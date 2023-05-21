#pragma once

#include "bloomRF/bloomRF.h"
#include <random>

namespace filters {

namespace test {

inline BloomFilterRFParameters genParams(size_t filterSizeBytes, int maxDelta, size_t maxDeltaSum) {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> layer(1, maxDelta);

  std::vector<size_t> layers((rand() % 8) + 2);
  std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  while (std::accumulate(layers.begin(), layers.end(), 0) > maxDeltaSum) {
    std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  }
  return BloomFilterRFParameters{filterSizeBytes, 0, layers};
}

}

}
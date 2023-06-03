#pragma once

#include <algorithm>
#include <random>
#include <sstream>
#include <vector>
#include "bloomRF/bloomRF.h"

namespace filters {

namespace test {

inline BloomFilterRFParameters genParams(size_t filterSizeBytes,
                                         int deltaSize,
                                         int maxDelta,
                                         size_t maxDeltaSum) {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> layer(1, maxDelta);

  std::vector<size_t> layers((rand() % deltaSize) + 2);
  std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  while (std::accumulate(layers.begin(), layers.end(), 0) > maxDeltaSum) {
    std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  }
  std::vector<int> r(layers.size());
  std::generate(r.begin(), r.end(), []() { return rand() % 2 + 1; });
  return BloomFilterRFParameters{filterSizeBytes, 0, layers, r};
}

template<typename T>
std::ostringstream genParamString(const BloomRF<T>& bf) {
  std::ostringstream os;
  os << "Delta vector: ";
  for (auto&& v : bf.getDelta()) {
    os << v << " ";
  }
  os << std::endl << "R vector: ";
  for (auto&& r : bf.getR()) {
    os << r << " ";
  }
  os << std::endl;
  return os;
}

}  // namespace test

}  // namespace filters
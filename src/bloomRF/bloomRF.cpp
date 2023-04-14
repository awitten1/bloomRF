#include "bloomRF.h"
#include <_types/_uint64_t.h>
#include <_types/_uint8_t.h>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <exception>

namespace filters {

namespace {

constexpr uint64_t SEED_GEN_A = 845897321;
constexpr uint64_t SEED_GEN_B = 217728422;

constexpr uint64_t MAX_BLOOM_FILTER_SIZE = 1 << 30;




}  // namespace

BloomFilterRFParameters::BloomFilterRFParameters(size_t filter_size_,
                                                 size_t filter_hashes_,
                                                 size_t seed_,
                                                 uint16_t delta_)
    : filter_size(filter_size_),
      filter_hashes(filter_hashes_),
      seed(seed_),
      delta(delta_) {
  if (filter_size == 0)
    throw "The size of bloom filter cannot be zero";
  if (filter_hashes == 0)
    throw "The number of hash functions for bloom filter cannot be zero";
}

template <typename T>
BloomRF<T>::BloomRF(const BloomFilterRFParameters& params)
    : BloomRF<T>(params.filter_size,
                 params.filter_hashes,
                 params.seed,
                 params.delta) {}

template <typename T>
size_t BloomRF<T>::bloomRFHash(T data, size_t i) {
  auto hash = this->hash(data >> (delta - 1), i);
  return ((hash % (numBits() >> (delta - 1))) << (delta - 1)) +
         (data & (1 << (delta - 1)) - 1);
}

template <typename T>
size_t BloomRF<T>::hash(T data, size_t i) {
  size_t hash1 = CityHash64WithSeed(reinterpret_cast<const char*>(&data),
                                    sizeof(data), seed);
  size_t hash2 =
      CityHash64WithSeed(reinterpret_cast<const char*>(&data), sizeof(data),
                         SEED_GEN_A * seed + SEED_GEN_B);
  return hash1 + i * hash2 + i * i;
}

template <typename T>
void BloomRF<T>::add(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHash(data, i) % numBits();
    filter[pos / (8 * sizeof(UnderType))] |=
        (1ULL << (pos % (8 * sizeof(UnderType))));
    data >>= delta;
  }
}

template <typename T>
bool BloomRF<T>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHash(data, i) % numBits();
    if (!(filter[pos / (8 * sizeof(UnderType))] &
          (1ULL << (pos % (8 * sizeof(UnderType)))))) {
      return false;
    }
    data >>= delta;
  }
  return true;
}

template<typename T>
std::vector<std::vector<T>> decomposeIntoDyadicIntervals(T lkey, T hkey) {
  if (lkey > hkey) {
    throw std::logic_error{"lkey must be less than hkey."};
  }
  uint16_t domain_width = sizeof(T) * 8;

  std::vector<std::vector<T>> ret(domain_width);

  uint16_t iterations = 0;

  T mid;
  T low = 0;
  T high = ~low;

  // Phase 1: Descend through DIs that fully cover the target range.
  for (; ;) {
    ret[iterations].push_back(low);

    mid = high - ((high - low) >> 1);
    if (mid <= lkey) {
      low = mid;
    } else if (mid >= hkey) {
      high = mid - 1;
    } else {
      break;
    }
    ++iterations;
  }
  assert(mid > lkey && mid < hkey);

  // Phase 2: Descend left and then descend right.
  uint8_t iterations_copy = iterations;
  T left_low = low, left_high = mid;
  for (; ;) {
    T left_mid = left_high - ((left_high - left_low) >> 1);
    if (lkey < left_mid) {
      ret[iterations_copy].push_back(left_mid);
      left_high = left_mid;
      ++iterations_copy;
    } else if (lkey > left_mid) {
      left_low = left_mid;
      ++iterations_copy;
    } else {
      ret[iterations_copy].push_back(left_mid);
      break;
    }
  }

  T right_low = mid, right_high = high;
  for (; ;) {
    T right_mid = right_high - ((right_high - right_low) >> 1);
    if (hkey > right_mid) {
      ret[iterations].push_back(right_low);
      right_low = right_mid;
      ++iterations;
    } else if (hkey < right_mid) {
      right_high = right_mid - 1;
      ++iterations;
    } else {
      ret[iterations].push_back(right_low);
      break;
    }
  }


  return ret;
}



template <typename T>
BloomRF<T>::BloomRF(size_t size_, size_t hashes_, size_t seed_, uint16_t delta_)
    : hashes(hashes_),
      seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)),
      delta(delta_),
      filter(words, 0) {}

template class BloomRF<uint16_t>;
template class BloomRF<uint32_t>;
template class BloomRF<uint64_t>;

template std::vector<std::vector<uint64_t>> decomposeIntoDyadicIntervals<uint64_t>(uint64_t, uint64_t);

}  // namespace filters

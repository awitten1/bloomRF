#include "bloomRF.h"
#include <_types/_uint16_t.h>
#include <_types/_uint64_t.h>
#include <_types/_uint8_t.h>
#include <cassert>
#include <cstddef>
#include <exception>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace filters {

namespace {

constexpr uint64_t SEED_GEN_A = 845897321;
constexpr uint64_t SEED_GEN_B = 217728422;

constexpr uint64_t MAX_BLOOM_FILTER_SIZE = 1 << 30;

}  // namespace

BloomFilterRFParameters::BloomFilterRFParameters(size_t filter_size_,
                                                 size_t filter_hashes_,
                                                 size_t seed_)
    : filter_size(filter_size_),
      filter_hashes(filter_hashes_),
      seed(seed_) {
  if (filter_size == 0)
    throw "The size of bloom filter cannot be zero";
  if (filter_hashes == 0)
    throw "The number of hash functions for bloom filter cannot be zero";
}

template <typename T, typename UnderType, size_t Delta>
BloomRF<T, UnderType, Delta>::BloomRF(const BloomFilterRFParameters& params)
    : BloomRF<T, UnderType>(params.filter_size,
                 params.filter_hashes,
                 params.seed) {}

template <typename T, typename UnderType, size_t Delta>
size_t BloomRF<T, UnderType, Delta>::bloomRFHash(T data, size_t i) {
  auto hash = this->hash(data >> (delta - 1), i);
  return ((hash % (numBits() >> (delta - 1))) << (delta - 1)) +
         (data & (1 << (delta - 1)) - 1);
}

template <typename T, typename UnderType, size_t Delta>
size_t BloomRF<T, UnderType, Delta>::hash(T data, size_t i) {
  size_t hash1 = CityHash64WithSeed(reinterpret_cast<const char*>(&data),
                                    sizeof(data), seed);
  size_t hash2 =
      CityHash64WithSeed(reinterpret_cast<const char*>(&data), sizeof(data),
                         SEED_GEN_A * seed + SEED_GEN_B);
  return hash1 + i * hash2 + i * i;
}

template <typename T, typename UnderType, size_t Delta>
void BloomRF<T, UnderType, Delta>::add(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHash(data, i) % numBits();
    filter[pos / (8 * sizeof(UnderType))] |=
        (1ULL << (pos % (8 * sizeof(UnderType))));
    data >>= delta;
  }
}

template <typename T, typename UnderType, size_t Delta>
bool BloomRF<T, UnderType, Delta>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHash(data, i) % numBits();
    if (!checkBit(pos)) {
      return false;
    }
    data >>= delta;
  }
  return true;
}

template <typename T, typename UnderType, size_t Delta>
bool BloomRF<T, UnderType, Delta>::findRange(T low, T high) {
  Checks checks(low, high);
  checks.initChecks(hashes, delta);

  for (size_t layer = hashes - 1; layer >= 0; --layer) {
    Checks new_checks(low, high);
    for (const auto& check : checks.getChecks()) {
      if (check.is_covering) {
        size_t pos = bloomRFHash(check.low, layer);
        if (checkBit(pos)) {
          Checks check_for_interval{low, high};
          check_for_interval.advanceChecks(delta);
          new_checks.concatenateChecks(check_for_interval);
        }
      } else {
        size_t pos = bloomRFHash(check.low, layer);
        UnderType word = getWord(pos);
        return (pos & word) != 0;
      }
    }
    checks = std::move(new_checks);
  }

  return false;
}

template<typename T, typename UnderType, size_t Delta>
void BloomRF<T, UnderType, Delta>::Checks::advanceChecks(size_t times) {
  for (int i = 0; i < times; ++i) {
    decltype(checks) new_checks;
    bool split = checks.size() > 1;  // Possibly a superfluous check.
    for (const auto& check : checks) {
      if (check.is_covering) {
        T low = check.low;
        T high = check.high;
        T mid = low + ((high - low) >> 1);
        if (mid <= lkey && !split) {
          new_checks.push_back({mid, high, true, false});
        } else if (mid >= hkey && !split) {
          new_checks.push_back({low, static_cast<T>(mid - 1), true, false});
        } else {
          if (check.on_left) {
            bool is_left_covering = mid > lkey;
            new_checks.push_back({mid, high, !is_left_covering, true});
            if (is_left_covering) {
              new_checks.push_back({low, static_cast<T>(mid - 1), true, true});
            }
          } else {
            bool is_right_covering = mid < hkey;
            new_checks.push_back(
                {low, static_cast<T>(mid - 1), !is_right_covering, false});
            if (is_right_covering) {
              new_checks.push_back({mid, high, true, false});
            }
          }
        }
      }
    }
    checks = std::move(new_checks);
  }
}

template <typename T, typename UnderType, size_t Delta>
void BloomRF<T, UnderType, Delta>::Checks::initChecks(size_t num_hashes, uint16_t delta) {
  static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>,
                "Must be integral and unsigned. Must provide alternate "
                "implementation for other types.");

  T low = 0;
  T high = ~low;

  checks.push_back({low, high, true, false});
  advanceChecks(domain_width - (num_hashes * delta));

}

template <typename T, typename UnderType, size_t Delta>
BloomRF<T, UnderType, Delta>::BloomRF(size_t size_, size_t hashes_, size_t seed_)
    : hashes(hashes_),
      seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)),
      delta(Delta),
      filter(words, 0) {}

template class BloomRF<uint16_t>;
template class BloomRF<uint32_t>;
template class BloomRF<uint64_t>;

}  // namespace filters

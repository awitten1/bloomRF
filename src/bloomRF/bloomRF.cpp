#include "bloomRF.h"
#include "city/city.h"
#include <_types/_uint16_t.h>
#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <_types/_uint8_t.h>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iostream>
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
    : filter_size(filter_size_), filter_hashes(filter_hashes_), seed(seed_) {
  if (filter_size == 0)
    throw "The size of bloom filter cannot be zero";
  if (filter_hashes == 0)
    throw "The number of hash functions for bloom filter cannot be zero";
}

template <typename T, typename UnderType, size_t Delta>
BloomRF<T, UnderType, Delta>::BloomRF(const BloomFilterRFParameters& params)
    : BloomRF<T, UnderType, Delta>(params.filter_size,
                            params.filter_hashes,
                            params.seed) {}

template <typename T, typename UnderType, size_t Delta>
size_t BloomRF<T, UnderType, Delta>::bloomRFHashToWord(T data, size_t i) {
  auto hash = this->hash(data >> ((i * delta) + delta - 1), i);
  return hash % (numBits() >> (Delta - 1));
}

template <typename T, typename UnderType, size_t Delta>
UnderType BloomRF<T, UnderType, Delta>::bloomRFRemainder(T data, size_t i) {
  T offset = (data >> (i * Delta)) & remainderMask;
  return 1 << offset;
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
    size_t pos = bloomRFHashToWord(data, i);
    filter[pos] |= bloomRFRemainder(data, i);
  }
}

template <typename T, typename UnderType, size_t Delta>
bool BloomRF<T, UnderType, Delta>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHashToWord(data, i);
    if (!(filter[pos] & bloomRFRemainder(data, i))) {
      return false;
    }
  }
  return true;
}

template <typename T, typename UnderType, size_t Delta>
bool BloomRF<T, UnderType, Delta>::findRange(T low, T high) {
  Checks checks(low, high, {});

  checks.initChecks(hashes, delta);

  for (int layer = hashes - 1; layer >= 0; --layer) {
    Checks new_checks(low, high, {});
    for (const auto& check : checks.getChecks()) {
      if (check.is_covering) {
        size_t pos = bloomRFHashToWord(check.low, layer);
        if (filter[pos] & bloomRFRemainder(check.low, layer)) {
          Checks check_for_interval{
              low,
              high,
              {typename Checks::Check{check.low, check.high, true, check.loc}}};
          check_for_interval.advanceChecks(delta);
          new_checks.concatenateChecks(check_for_interval);
        }
      } else {
        UnderType bitmask = buildBitMaskForRange(check.low, check.high, layer);
        UnderType word = filter[bloomRFHashToWord(check.low, layer)];
        if ((bitmask & word) != 0) {
          return true;
        }
      }
    }
    checks = std::move(new_checks);
  }

  return false;
}

template <typename T, typename UnderType, size_t Delta>
void BloomRF<T, UnderType, Delta>::Checks::advanceChecks(size_t times) {
  for (int i = 0; i < times; ++i) {
    decltype(checks) new_checks;
    bool split = checks.size() > 1;
    for (const auto& check : checks) {
      if (check.is_covering) {
        T low = check.low;
        T high = check.high;
        T mid = high - ((high - low) >> 1);
        if (check.loc == IntervalLocation::NotYetSplit) {
          if (mid <= lkey) {
            bool covering = mid < lkey || high > hkey;
            new_checks.push_back(
                {mid, high, covering, IntervalLocation::NotYetSplit});
          } else if (mid - 1 >= hkey) {
            bool covering = low < lkey || mid - 1 > hkey;
            new_checks.push_back({low, static_cast<T>(mid - 1), covering,
                                  IntervalLocation::NotYetSplit});
          } else {
            new_checks.push_back({low, static_cast<T>(mid - 1), lkey < low,
                                  IntervalLocation::Left});
            new_checks.push_back(
                {mid, high, hkey > high, IntervalLocation::Right});
          }
        } else if (check.loc == IntervalLocation::Left) {
          bool is_left_covering = mid > lkey;
          new_checks.push_back(
              {mid, high, !is_left_covering, IntervalLocation::Left});
          if (is_left_covering) {
            new_checks.push_back(
                {low, static_cast<T>(mid - 1), true, IntervalLocation::Left});
          }
        } else {
          bool is_right_covering = mid < hkey;
          new_checks.push_back({low, static_cast<T>(mid - 1),
                                !is_right_covering, IntervalLocation::Right});
          if (is_right_covering) {
            new_checks.push_back({mid, high, true, IntervalLocation::Right});
          }
        }
      } else {
        new_checks.push_back(check);
      }
    }
    checks = std::move(new_checks);
  }
}

template <typename T, typename UnderType, size_t Delta>
void BloomRF<T, UnderType, Delta>::Checks::initChecks(size_t num_hashes,
                                                      uint16_t delta) {
  T low = 0;
  T high = ~low;

  if (checks.size() > 0) {
    throw std::logic_error{
        "Cannot init checks on a non-empty checks instance."};
  }

  checks.push_back({low, high, true, IntervalLocation::NotYetSplit});
  advanceChecks(domain_width - ((num_hashes - 1) * delta));
}

template <typename T, typename UnderType, size_t Delta>
BloomRF<T, UnderType, Delta>::BloomRF(size_t size_,
                                      size_t hashes_,
                                      size_t seed_)
    : hashes(hashes_),
      seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)),
      delta(Delta),
      filter(words, 0) {}

template class BloomRF<uint16_t>;
template class BloomRF<uint32_t>;
template class BloomRF<uint64_t>;
template class BloomRF<uint64_t, uint32_t, 6>;



#ifdef __SIZEOF_INT128__
template class BloomRF<uint64_t, uint128_t, 8>;
#endif


}  // namespace filters

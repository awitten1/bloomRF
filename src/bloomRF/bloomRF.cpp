#include "bloomRF.h"
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
#include <cstdlib>
#include "city/city.h"

namespace filters {

namespace {

constexpr uint64_t SEED_GEN_A = 845897321;
constexpr uint64_t SEED_GEN_B = 217728422;

constexpr uint64_t MAX_BLOOM_FILTER_SIZE = 1 << 30;

}  // namespace

BloomFilterRFParameters::BloomFilterRFParameters(size_t filter_size_,
                                                 size_t seed_,
                                                 std::vector<size_t> delta_)
    : filter_size(filter_size_), seed(seed_), delta(std::move(delta_)) {
  if (filter_size == 0)
    throw "The size of bloom filter cannot be zero";
}

template <typename T, typename UnderType>
BloomRF<T, UnderType>::BloomRF(const BloomFilterRFParameters& params)
    : BloomRF<T, UnderType>(params.filter_size,
                                   params.seed, params.delta) {}

template <typename T, typename UnderType>
size_t BloomRF<T, UnderType>::bloomRFHashToWord(T data, size_t i) {
  auto hash = this->hash(data >> (shifts[i] + delta[i] - 1), i);
  return hash % (numBits() >> (delta[i] - 1));
}

template <typename T, typename UnderType>
UnderType BloomRF<T, UnderType>::bloomRFRemainder(T data, size_t i, int wordPos) {
  T offset = (data >> shifts[i]) & ((1 << (delta[i] - 1)) - 1);
  return (1 << offset) << (wordPos * (delta[i] - 1));
}


template <typename T, typename UnderType>
size_t BloomRF<T, UnderType>::hash(T data, size_t i) {
  size_t hash1 = CityHash64WithSeed(reinterpret_cast<const char*>(&data),
                                    sizeof(data), seed);
  size_t hash2 =
      CityHash64WithSeed(reinterpret_cast<const char*>(&data), sizeof(data),
                         SEED_GEN_A * seed + SEED_GEN_B);
  return hash1 + i * hash2 + i * i;
}

template <typename T, typename UnderType>
void BloomRF<T, UnderType>::add(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHashToWord(data, i);
    std::div_t div = std::div(8 * sizeof(UnderType), (1 << (delta[i] - 1)));
    filter[pos / div.quot] |= bloomRFRemainder(data, i, div.rem);
  }
}

template <typename T, typename UnderType>
bool BloomRF<T, UnderType>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHashToWord(data, i);
    std::div_t div = std::div(8 * sizeof(UnderType), (1 << (delta[i] - 1)));
    if (!(filter[pos / div.quot] & bloomRFRemainder(data, i, div.rem))) {
      return false;
    }
  }
  return true;
}

template <typename T, typename UnderType>
bool BloomRF<T, UnderType>::findRange(T low, T high) {
  Checks checks(low, high, {});

  checks.initChecks(shifts.back());

  for (int layer = hashes - 1; layer >= 0; --layer) {
    Checks new_checks(low, high, {});
    for (const auto& check : checks.getChecks()) {
      if (check.is_covering) {
        size_t pos = bloomRFHashToWord(check.low, layer);
        std::div_t div = std::div(8 * sizeof(UnderType), (1 << (delta[layer] - 1)));
        if (filter[pos / div.quot] & bloomRFRemainder(check.low, layer, div.rem)) {
          Checks check_for_interval{
              low,
              high,
              {typename Checks::Check{check.low, check.high, true, check.loc}}};
          check_for_interval.advanceChecks(delta[layer]);
          new_checks.concatenateChecks(check_for_interval);
        }
      } else {
        size_t pos = bloomRFHashToWord(check.low, layer);
        std::div_t div = std::div(8 * sizeof(UnderType), (1 << (delta[layer] - 1)));
        UnderType bitmask = buildBitMaskForRange(check.low, check.high, layer, div.rem);
        UnderType word = filter[pos / div.quot];
        if ((bitmask & word) != 0) {
          return true;
        }
      }
    }
    checks = std::move(new_checks);
  }

  return false;
}

template <typename T, typename UnderType>
void BloomRF<T, UnderType>::Checks::advanceChecks(size_t times) {
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

template <typename T, typename UnderType>
void BloomRF<T, UnderType>::Checks::initChecks(size_t delta_sum) {
  T low = 0;
  T high = ~low;

  if (checks.size() > 0) {
    throw std::logic_error{
        "Cannot init checks on a non-empty checks instance."};
  }

  checks.push_back({low, high, true, IntervalLocation::NotYetSplit});
  advanceChecks(domain_width - delta_sum);
}

template <typename T, typename UnderType>
BloomRF<T, UnderType>::BloomRF(size_t size_,
                                      size_t seed_,
                                      std::vector<size_t> delta_)
    : hashes(delta_.size()),
      seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)),
      filter(words, 0), delta(delta_), shifts(delta.size()) {

  for (auto d : delta) {
    if ((d - 1) > 8 * sizeof(UnderType)) {
      throw std::logic_error{"The width of a PMHF word cannot exceed the width of the UnderType; \
        For all layer widths d, (d - 1) cannot exceed 8 * sizeof(UnderType)."};
    }
  }

  /// Compute prefix sums.
  for (int i = 1; i < delta.size(); ++i) {
    shifts[i] = shifts[i - 1] + delta[i - 1];
  }

}

template class BloomRF<uint16_t>;
template class BloomRF<uint32_t>;
template class BloomRF<uint64_t>;
template class BloomRF<uint64_t, uint32_t>;

#ifdef __SIZEOF_INT128__
template class BloomRF<uint64_t, uint128_t>;
#endif

}  // namespace filters

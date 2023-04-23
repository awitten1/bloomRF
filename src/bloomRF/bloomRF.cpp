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
#include <sstream>
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

template<typename T>
void printBinary(T t) {
  std::vector<int> print;
  auto x = 8 * sizeof(T);
  while(x--) {
    print.push_back(t % 2);
    t >>= 1;
  }
  std::reverse(print.begin(), print.end());
  std::copy(print.begin(), print.end(), std::ostream_iterator<int>(std::cerr));
  std::cerr << std::endl;
}

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
  UnderType offset = (data >> shifts[i]) & ((1 << (delta[i] - 1)) - 1);
  UnderType ret = (UnderType{1} << offset);
  ret = ret << (wordPos * (1 << (delta[i] - 1)));
  return ret;
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
    //printBinary(data >> shifts[i]);
    auto div = getFilterPosAndOffset(pos, i);
    //std::cout << "pos: " << div.quot << ", offset: " << div.rem << std::endl;
    //std::cout << "layer: " << i << std::endl;
    filter[div.quot] |= bloomRFRemainder(data, i, div.rem);
  }
  //std::cout << "---- done adding ------" << std::endl;
}

template <typename T, typename UnderType>
bool BloomRF<T, UnderType>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = bloomRFHashToWord(data, i);
    auto div = getFilterPosAndOffset(pos, i);
    if (!(filter[div.quot] & bloomRFRemainder(data, i, div.rem))) {
      return false;
    }
  }
  return true;
}

template<typename T, typename UnderType>
std::ldiv_t BloomRF<T, UnderType>::getFilterPosAndOffset(size_t pos, size_t i) {
  size_t wordsPerUnderType = 8 * sizeof(UnderType) / (1 << (delta[i] - 1));
  return std::ldiv(pos, wordsPerUnderType);
}

template <typename T, typename UnderType>
bool BloomRF<T, UnderType>::findRange(T lkey, T hkey) {
  Checks checks(lkey, hkey, {});

  checks.initChecks(shifts.back());

  for (int layer = hashes - 1; layer >= 0; --layer) {
    Checks new_checks(lkey, hkey, {});
    //std::cout << "layer: " << layer << std::endl;
    for (const auto& check : checks.getChecks()) {
      //std::cout << "low: " << check.low << ", high: " << check.high << std::endl;
      //printBinary(check.low);
      if (check.low < lkey || check.high > hkey) {
        assert((check.high - check.low + 1 == (T{1} << shifts[layer])));
        //std::cout << "is covering" << std::endl;
        size_t pos = bloomRFHashToWord(check.low, layer);
        auto div = getFilterPosAndOffset(pos, layer);
        if (filter[div.quot] & bloomRFRemainder(check.low, layer, div.rem)) {
          Checks check_for_interval{
              lkey,
              hkey,
              {typename Checks::Check{check.low, check.high, check.loc}}};
          check_for_interval.advanceChecks(delta[layer - 1]);
          new_checks.concatenateChecks(check_for_interval);
        }
      } else {
        //std::cout << "not covering" << std::endl;
        //std::cout << check.high - check.low << std::endl;
        size_t pos = bloomRFHashToWord(check.low, layer);
        auto div = getFilterPosAndOffset(pos, layer);
        //printBinary(check.low);
        //std::cout << "pos: " << div.quot << ", offset: " << div.rem << std::endl;
        UnderType bitmask = buildBitMaskForRange(check.low, check.high, layer, div.rem);
        //printBinary(bitmask);
        UnderType word = filter[div.quot];
        //printBinary(word);
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
    //std::cout << "iter: " << i << std::endl;
    for (const auto& check : checks) {
      if (check.low < lkey || check.high > hkey) {
        T low = check.low;
        T high = check.high;
        T mid = high - ((high - low) >> 1);
        if (check.loc == IntervalLocation::NotYetSplit) {
          if (mid <= lkey) {
            new_checks.push_back(
                {mid, high, IntervalLocation::NotYetSplit});
          } else if (mid - 1 >= hkey) {
            bool covering = low < lkey || mid - 1 > hkey;
            new_checks.push_back({low, static_cast<T>(mid - 1),
                                  IntervalLocation::NotYetSplit});
          } else {
            new_checks.push_back({low, static_cast<T>(mid - 1),
                                  IntervalLocation::Left});
            new_checks.push_back(
                {mid, high, IntervalLocation::Right});
          }
        } else if (check.loc == IntervalLocation::Left) {
          bool is_left_covering = mid > lkey;
          new_checks.push_back(
              {mid, high, IntervalLocation::Left});
          if (is_left_covering) {
            new_checks.push_back(
                {low, static_cast<T>(mid - 1), IntervalLocation::Left});
          }
        } else {
          bool is_right_covering = mid <= hkey;
          new_checks.push_back({low, static_cast<T>(mid - 1),
                                 IntervalLocation::Right});
          if (is_right_covering) {
            new_checks.push_back({mid, high, IntervalLocation::Right});
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

  checks.push_back({low, high, IntervalLocation::NotYetSplit});
  advanceChecks(domain_width - delta_sum);
}

template<typename T, typename UnderType>
UnderType BloomRF<T, UnderType>::buildBitMaskForRange(T low, T high, size_t i, int wordPos) {
    UnderType lowOffset = ((low >> shifts[i]) & ((1 << (delta[i] - 1)) - 1));
    UnderType highOffset = ((high >> shifts[i]) & ((1 << (delta[i] - 1)) - 1));
    UnderType ret = ~0;
    ret ^= (UnderType{1} << lowOffset) - 1;

    if (highOffset < 8 * sizeof(UnderType) - 1) {
      ret &= (UnderType{1} << (highOffset + 1)) - 1;
    }

    ret = ret << (wordPos * (1 << (delta[i] - 1)));
    return ret;
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
    if (8 * sizeof(UnderType) % (1 << (d - 1))) {
      std::ostringstream os;
      os << "The width of a PMHF word, must divide the width of the UnderType; ";
      os << "this does not hold for distance layer: " << d << " and UnderType width: "
         << 8 * sizeof(UnderType)
         << std::endl;
      throw std::logic_error{os.str()};
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

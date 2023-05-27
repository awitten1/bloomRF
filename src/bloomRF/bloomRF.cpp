#include "bloomRF.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
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
    throw std::logic_error{"The size of bloom filter cannot be zero"};
}

namespace detail {

template <typename T, typename UnderType>
BloomRfImpl<T, UnderType>::BloomRfImpl(const BloomFilterRFParameters& params)
    : BloomRfImpl<T, UnderType>(params.filter_size, params.seed, params.delta) {
}

template <typename T, typename UnderType>
size_t BloomRfImpl<T, UnderType>::bloomRFHashToWord(T data, size_t i) const {
  auto hash = this->hash(data >> (shifts[i] + delta[i] - 1), i);
  return hash % (numBits() >> (delta[i] - 1));
}

template <typename T, typename UnderType>
UnderType BloomRfImpl<T, UnderType>::bloomRFRemainder(T data,
                                                      size_t i,
                                                      int wordPos) const {
  UnderType offset =
      (data >> shifts[i]) & ((UnderType{1} << (delta[i] - 1)) - 1);
  UnderType ret = (UnderType{1} << offset);
  assert(offset < (1 << (delta[i] - 1)));

  ret = ret << (wordPos * (1 << (delta[i] - 1)));

  return ret;
}

template <typename T, typename UnderType>
size_t BloomRfImpl<T, UnderType>::hash(T data, size_t i) const {
  size_t hash1 = CityHash64WithSeed(reinterpret_cast<const char*>(&data),
                                    sizeof(data), seed);
  size_t hash2 =
      CityHash64WithSeed(reinterpret_cast<const char*>(&data), sizeof(data),
                         SEED_GEN_A * seed + SEED_GEN_B);
  return hash1 + i * hash2 + i * i;
}

template <typename T, typename UnderType>
void BloomRfImpl<T, UnderType>::add(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    auto hash = hashToIndexAndBitMask(data, i);
    filter[hash.first] |= hash.second;
  }
}

template <typename T, typename UnderType>
bool BloomRfImpl<T, UnderType>::find(T data) const {
  for (size_t i = 0; i < hashes; ++i) {
    const auto& [filterPos, bitmask] = hashToIndexAndBitMask(data, i);
    if (!(filter[filterPos] & bitmask)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename UnderType>
std::pair<size_t, UnderType> BloomRfImpl<T, UnderType>::hashToIndexAndBitMask(
    T data,
    size_t i) const {
  size_t pos = bloomRFHashToWord(data, i);

  if (1 << (delta[i] - 1) <= 8 * sizeof(UnderType)) {
    // Case 1: Size of PMHF word is less than or equal to the size of the
    // UnderType.
    size_t wordsPerUnderType = 8 * sizeof(UnderType) / (1 << (delta[i] - 1));
    std::ldiv_t div = std::ldiv(pos, wordsPerUnderType);
    return {div.quot, bloomRFRemainder(data, i, div.rem)};

  } else {
    // Case 2: Size of PMHF word is greater than that of the UnderType.
    int pmhfWordsPerUT = (1 << (delta[i] - 1)) / (8 * sizeof(UnderType));
    auto filterPos = pos * pmhfWordsPerUT;
    UnderType offset =
        ((data >> shifts[i]) & ((UnderType{1} << (delta[i] - 1)) - 1));
    assert(offset < (1 << (delta[i] - 1)));
    std::ldiv_t div = std::ldiv(offset, 8 * sizeof(UnderType));
    assert(div.rem < 8 * sizeof(UnderType));

    filterPos += div.quot;
    return {filterPos, (UnderType{1} << div.rem)};
  }
}

template <typename T, typename UnderType>
bool BloomRfImpl<T, UnderType>::checkDIOfDecomposition(T low,
                                                       T high,
                                                       int layer) const {
  size_t pos = bloomRFHashToWord(low, layer);

  if (1 << (delta[layer] - 1) <= 8 * sizeof(UnderType)) {
    // Case 1: Size of PMHF word is less than or equal to the size of the
    // UnderType.
    size_t wordsPerUnderType =
        8 * sizeof(UnderType) / (1 << (delta[layer] - 1));
    std::ldiv_t div = std::ldiv(pos, wordsPerUnderType);
    UnderType bitmask = buildBitMaskForRange(low, high, layer, div.rem);
    UnderType word = filter[div.quot];
    if ((bitmask & word) != 0) {
      return true;
    }
  } else {
    // Case 2: Size of PMHF word is greater than that of the UnderType.
    // In this case we need to iterate over the UnderTypes that comprise the
    // PMHF word.
    int pmhfWordsPerUT = (1 << (delta[layer] - 1)) / (8 * sizeof(UnderType));
    int filterPos = pos * pmhfWordsPerUT;
    UnderType lowOffset =
        ((low >> shifts[layer]) & ((UnderType{1} << (delta[layer] - 1)) - 1));
    filterPos += (lowOffset / (8 * sizeof(UnderType)));
    UnderType highOffset =
        ((high >> shifts[layer]) & ((UnderType{1} << (delta[layer] - 1)) - 1));
    size_t iters = (highOffset / (8 * sizeof(UnderType))) -
                   (lowOffset / (8 * sizeof(UnderType))) + 1;
    for (int i = 0; i < iters; ++i) {
      UnderType bitmask = ~UnderType{0};
      if (i == 0) {
        bitmask ^= (UnderType{1} << (lowOffset % (8 * sizeof(UnderType)))) - 1;
      }
      if (i == iters - 1 && (highOffset % (8 * sizeof(UnderType))) <
                                (8 * sizeof(UnderType) - 1)) {
        bitmask &=
            (UnderType{1} << ((highOffset % (8 * sizeof(UnderType))) + 1)) - 1;
      }
      if ((bitmask & filter[filterPos]) != 0) {
        return true;
      }
      ++filterPos;
    }
  }
  return false;
}

template <typename T, typename UnderType>
bool BloomRfImpl<T, UnderType>::findRange(T lkey, T hkey) const {
  Checks checks(lkey, hkey, {});

  checks.initChecks(shifts.back(), delta.back());

  for (int layer = hashes - 1; layer >= 0; --layer) {
    Checks new_checks(lkey, hkey, {});
    //checks.compressChecks(shifts[layer] + delta[layer] - 1);

    for (const auto& check : checks.getChecks()) {
            std::cout << "[" << check.low << "," << check.high << "], ";
      if (check.low < lkey || check.high > hkey) {
        auto hash = hashToIndexAndBitMask(check.low, layer);
        if (filter[hash.first] & hash.second) {
          Checks check_for_interval{
              lkey,
              hkey,
              {typename Checks::Check{check.low, check.high}}};
          check_for_interval.advanceChecks(shifts[layer - 1], delta[layer - 1]);
          new_checks.concatenateChecks(check_for_interval);
        }
      } else {
        if (checkDIOfDecomposition(check.low, check.high, layer)) {
          return true;
        }
      }
    }
    std::cout << std::endl;

    checks = std::move(new_checks);
  }

  return false;
}

template <typename T, typename UnderType>
void BloomRfImpl<T, UnderType>::Checks::advanceChecks(size_t shifts, size_t delta) {
  T target_width = 1 << shifts;
  std::vector<Check> new_checks;
  for (const auto& check : checks) {
    T low = (lkey / target_width) * target_width;
    for (T counter = low; counter <= check.high; counter += target_width) {
      T curr_high = counter + target_width;
      if (!new_checks.empty() && new_checks.back().low >= lkey && curr_high <= hkey &&
        (counter >> (shifts + delta - 1) == new_checks.back().low >> (shifts + delta - 1))) {
          new_checks.back().high = curr_high;
      } else {
        new_checks.push_back({counter, curr_high});
      }
    }
  }
  checks = std::move(new_checks);
}

template <typename T, typename UnderType>
void BloomRfImpl<T, UnderType>::Checks::initChecks(size_t delta_sum,
                                                   size_t delta_back) {
  T low = 0;
  T high = ~low;

  if (checks.size() > 0) {
    throw std::logic_error{
        "Cannot init checks on a non-empty checks instance."};
  }

  checks.push_back({lkey, hkey});
  advanceChecks(delta_sum, delta_back);
}

template <typename T, typename UnderType>
UnderType BloomRfImpl<T, UnderType>::buildBitMaskForRange(T low,
                                                          T high,
                                                          size_t i,
                                                          int wordPos) const {
  UnderType lowOffset = ((low >> shifts[i]) & ((1 << (delta[i] - 1)) - 1));
  UnderType highOffset = ((high >> shifts[i]) & ((1 << (delta[i] - 1)) - 1));
  UnderType bitmask = ~UnderType{0};
  bitmask ^= (UnderType{1} << lowOffset) - 1;
  if (highOffset < 8 * sizeof(UnderType) - 1) {
    bitmask &= (UnderType{1} << (highOffset + 1)) - 1;
  }
  bitmask = bitmask << (wordPos * (1 << (delta[i] - 1)));
  return bitmask;
}

template <typename T, typename UnderType>
BloomRfImpl<T, UnderType>::BloomRfImpl(size_t size_,
                                       size_t seed_,
                                       std::vector<size_t> delta_)
    : hashes(delta_.size()),
      seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)),
      filter(words, 0),
      delta(delta_),
      shifts(delta.size()) {
  if (delta.empty()) {
    throw std::logic_error{"Delta vector cannot be empty."};
  }

  if (std::accumulate(delta.begin(), delta.end(), 0) > 8 * sizeof(T)) {
    throw std::logic_error{
        "Sum of delta vector should not exceed width of key."};
  }

  for (const auto& d : delta) {
    if (d == 0) {
      throw std::logic_error{"Delta vector cannot have zero values."};
    }
  }

  /// Compute prefix sums.
  for (int i = 1; i < delta.size(); ++i) {
    shifts[i] = shifts[i - 1] + delta[i - 1];
  }
}

template class BloomRfImpl<uint16_t>;
template class BloomRfImpl<uint32_t>;
template class BloomRfImpl<uint64_t>;
template class BloomRfImpl<uint64_t, uint32_t>;

} // namespace detail

}  // namespace filters

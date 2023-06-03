

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>
#include <cassert>

#include "bloomRF/bloomRF.h"

using filters::BloomFilterRFParameters;
using filters::BloomRF;

template <typename T, typename Generator, typename QueryGenerator>
class ExperimentDriver {

  static_assert(std::is_same_v<T, decltype(std::declval<Generator>()())>);
  static_assert(std::is_same_v<T, decltype(std::declval<QueryGenerator>()())>);

 public:
  ExperimentDriver(const BloomFilterRFParameters& params, Generator g, QueryGenerator qg)
      : gen(rd()), bf{params}, keyGenerator(g), queryKeyGenerator(qg) { }

  std::pair<bool, bool> find(T data) {
    return {bf.find(data), std::binary_search(s.begin(), s.end(), data)};
  }

  std::pair<bool, bool> findRange(T low, T high) {
    return {bf.findRange(low, high),
        [=, this]() {
          auto lower_bound = std::lower_bound(s.begin(), s.end(), low);
          if (lower_bound == s.end()) {
            return false;
          }
          return *lower_bound <= high;
        }()
    };
  }

  void insert(T data) {
    bf.add(data);
    s.push_back(data);
  }

  void doInserts(int n) {
    for (int i = 0; i < n; ++i) {
      T x = keyGenerator();
      insert(x);
    }
  }

  double randomQuerys(int denominator) {
    int false_positives = 0;
    int true_negative = 0;

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < denominator; ++i) {
      T q = keyGenerator();
      const auto& [inFilter, actuallyIn] = find(q);
      if (inFilter) {
        if (!actuallyIn) {
          ++false_positives;
        }
      }
      if (!inFilter) {
        ++true_negative;
      }

      // Extra sanity check.
      if (actuallyIn) {
        assert(inFilter);
      }
    }
    auto fp = static_cast<double>(false_positives) /
              static_cast<double>(false_positives + true_negative);
    return fp;
  }

  /// Assumes genQuery always returns a range that doesn't contain a key.
  double randomRangeQuerys(int denominator,
                           T interval_size) {
    int false_positives = 0;
    int true_negative = 0;

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < denominator; ++i) {
      T low = keyGenerator();
      T high = low + interval_size;
      if (high < low) high = std::numeric_limits<T>::max();
      const auto& [inFilter, actuallyIn] = findRange(low, high);
      if (inFilter) {
        if (!actuallyIn) {
          ++false_positives;
        }
      }
      if (!inFilter) {
        ++true_negative;
      }

      // Extra sanity check.
      if (actuallyIn) {
        assert(inFilter);
      }
    }

    auto fp =
        static_cast<double>(false_positives) / static_cast<double>(false_positives + true_negative);
    return fp;
  }

 private:

  std::random_device rd;
  std::mt19937 gen;
  BloomRF<T> bf;
  std::vector<T> s;
  Generator keyGenerator;
  QueryGenerator queryKeyGenerator;
};

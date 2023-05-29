

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bloomRF/bloomRF.h"

using filters::BloomFilterRFParameters;
using filters::BloomRF;

template <typename T, typename Generator>
class ExperimentDriver {

 static_assert(std::is_same_v<T, decltype(std::declval<Generator>()())>);

 public:
  ExperimentDriver(const BloomFilterRFParameters& params, Generator g)
      : gen(rd()), bf{params}, keyGenerator(g) { }

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
      insert(static_cast<T>(x));
    }
    std::sort(s.begin(), s.end());
  }

  double randomQuerys(int denominator) {
    int false_positives = 0;
    int true_positive = 0;

    for (int i = 0; i < denominator; ++i) {
      T q = keyGenerator();
      const auto& [inFilter, in] = find(q);
      if (inFilter) {
        if (!in) {
          ++false_positives;
        } else {
          ++true_positive;
        }
      }
    }
    auto fp = static_cast<double>(false_positives) /
              static_cast<double>(denominator - true_positive);
    return fp;
  }

  /// Assumes genQuery always returns a range that doesn't contain a key.
  double randomRangeQuerys(int denominator,
                           T interval_size) {
    int false_positives = 0;
    int true_positive = 0;

    for (int i = 0; i < denominator; ++i) {
      T low = keyGenerator();
      T high = low + interval_size;
      if (high < low) high = std::numeric_limits<T>::max();
      const auto& [inFilter, in] = findRange(low, high);
      if (inFilter) {
        if (!in) {
          ++false_positives;
        } else {
          ++true_positive;
        }
      }
    }
    auto fp =
        static_cast<double>(false_positives) / static_cast<double>(denominator - true_positive);
    return fp;
  }

 private:

  std::random_device rd;
  std::mt19937 gen;
  BloomRF<T> bf;
  std::vector<T> s;
  Generator keyGenerator;
};

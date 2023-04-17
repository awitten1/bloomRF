

#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>

#include "bloomRF/bloomRF.h"

using filters::BloomFilterRFParameters;
using filters::BloomRF;

template <typename T, typename RandomNumberGenerator, typename UnderType, size_t Delta>
class ExperimentDriver {
 public:
  ExperimentDriver(const BloomFilterRFParameters& params,
                   RandomNumberGenerator dt)
      : bf{params}, dist{std::move(dt)} {}

  std::pair<bool, bool> find(T data) {
    return {bf.find(data), s.find(data) != s.end()};
  }

  void insert(T data) {
    bf.add(data);
    s.insert(data);
  }

  void doInserts(int n) {
    for (int i = 0; i < n; ++i) {
      insert(static_cast<T>(std::round(dist())));
    }
  }

  double randomQuerys(int denominator) {
    int false_positives = 0;
    int true_positive = 0;

    for (int i = 0; i < denominator; ++i) {
      T q = dist();
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

 private:
  BloomRF<T, UnderType, Delta> bf;
  std::unordered_set<T> s;
  RandomNumberGenerator dist;
};

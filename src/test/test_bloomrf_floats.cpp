#include <gtest/gtest.h>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include "gtest/gtest.h"

#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <ostream>
#include <random>
#include <unordered_set>

#include "city/city.h"
#include "test_helpers.h"

namespace filters {
namespace test {

template <typename FloatingPoint>
FloatingPoint randomUniformFloat(FloatingPoint low, FloatingPoint high) {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_real_distribution<FloatingPoint> dist(low, high);
  return dist(rd);
}

class BloomFilterUniformFloat : public ::testing::Test,
                                public testing::WithParamInterface<
                                    std::pair<int, BloomFilterRFParameters>> {
 protected:
  void SetUp() override {
    for (int i = 0; i < GetParam().first; ++i) {
      auto x = randomUniformFloat(std::numeric_limits<float>::min(),
                                  std::numeric_limits<float>::max());
      s.push_back(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::vector<float> s;

  BloomRF<float> bf{GetParam().second};
};

TEST_P(BloomFilterUniformFloat, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_P(BloomFilterUniformFloat, NoFalseNegativesRangeQuerySmallRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - randomUniformFloat(.001, .01);
    auto high = *it + randomUniformFloat(.001, .01);
    ;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

TEST_P(BloomFilterUniformFloat, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - 1;
    auto high = *it + 1;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

TEST_P(BloomFilterUniformFloat, NoFalseNegativesRangeQueryExtraLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10;
    auto high = *it + rand() % 10;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniformFloat);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniformFloat,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        size_t numKeys = 10000;
        return std::pair<int, BloomFilterRFParameters>{
            numKeys, genParams((rand() % numKeys) + numKeys, 5, 7, 32)};
      });
      return ret;
    }()));

TEST(OneOffFloat, RangeQuery) {
  BloomRF<float, uint64_t> bf{
      BloomFilterRFParameters{16000, 0, {7, 6, 6, 4, 3}}};
  float key = 0;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 1, key + 1));
  ASSERT_TRUE(bf.findRange(key - .0001, key + .0001));
}

}  // namespace test
}  // namespace filters
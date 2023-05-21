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

#include "test_helpers.h"
#include "city/city.h"

namespace filters {
namespace test {

int64_t randomUniformInt64() {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<int64_t> dist(
      std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max());
  return dist(rd);
}

class BloomFilterUniformSigned64Test : public ::testing::Test,
                                  public testing::WithParamInterface<
                                      std::pair<int, BloomFilterRFParameters>> {
 protected:
  void SetUp() override {
    for (int i = 0; i < GetParam().first; ++i) {
      auto x = randomUniformInt64();
      s.push_back(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::vector<int64_t> s;

  BloomRF<int64_t> bf{GetParam().second};
};

TEST_P(BloomFilterUniformSigned64Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_P(BloomFilterUniformSigned64Test, NoFalseNegativesRangeQuerySmallRange) {
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

TEST_P(BloomFilterUniformSigned64Test, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10000;
    auto high = *it + rand() % 10000;
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

TEST_P(BloomFilterUniformSigned64Test, NoFalseNegativesRangeQueryExtraLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 100000;
    auto high = *it + rand() % 100000;
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

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniformSigned64Test);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniformSigned64Test,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        size_t numKeys = 10000;
        return std::pair<int, BloomFilterRFParameters>{numKeys, genParams((rand() % numKeys) + numKeys, 9, 64)};
      });
      return ret;
    }()));

}

}

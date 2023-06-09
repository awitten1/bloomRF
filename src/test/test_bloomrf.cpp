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

uint64_t randomUniformUint64() {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  return dist(rd);
}

class BloomFilterUniform32Test : public ::testing::Test,
                                 public testing::WithParamInterface<
                                     std::pair<int, BloomFilterRFParameters>> {
 protected:
  void SetUp() override {
    for (int i = 0; i < GetParam().first; ++i) {
      auto x = randomUniformUint64();
      s.push_back(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::vector<uint64_t> s;

  BloomRF<uint64_t, uint32_t> bf{GetParam().second};
};

class BloomFilterUniform64Test : public ::testing::Test,
                                 public testing::WithParamInterface<
                                     std::pair<int, BloomFilterRFParameters>> {
 protected:
  void SetUp() override {
    for (int i = 0; i < GetParam().first; ++i) {
      auto x = randomUniformUint64();
      s.push_back(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::vector<uint64_t> s;

  BloomRF<uint64_t, uint64_t> bf{GetParam().second};
};

TEST_P(BloomFilterUniform64Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_P(BloomFilterUniform64Test, NoFalseNegativesRangeQuerySmallRange) {
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

TEST_P(BloomFilterUniform64Test, NoFalseNegativesRangeQueryLargeRange) {
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

TEST_P(BloomFilterUniform64Test, NoFalseNegativesRangeQueryExtraLargeRange) {
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

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniform64Test);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniform64Test,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        size_t numKeys = 10000;
        return std::pair<int, BloomFilterRFParameters>{
            numKeys, genParams((rand() % numKeys) + numKeys, 8, 9, 64)};
      });
      return ret;
    }()));

TEST(OneOff, RangeQuery2) {
  BloomRF<uint64_t, uint64_t> bf{
      BloomFilterRFParameters{16000, 0, {8, 3, 3, 4}}};
  uint64_t key = 3068990209559152388;

  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 9, key + 2));
}

TEST(OneOff, RangeQuery) {
  BloomRF<uint64_t, uint64_t> bf{BloomFilterRFParameters{16000, 0, {9, 8, 6}}};
  uint64_t key = 17183560791176864955ULL;
  // bf.add(key);
  // ASSERT_TRUE(bf.findRange(key - 100, key + 100));

  key = 13539885930325430328ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 9, key + 9));

  key = 13482642926757329959ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 8, key + 9));

  key = 4944684668419138897ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 5, key + 8));

  key = 12836727673998169215ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 4, key + 2));

  key = 6734315744289451875ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key, key + 1));

  key = 16343179362131379382ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key, key));

  key = 1894361899248432030ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key, key + 1));

  key = 994988673032400334ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 3, key));

  key = 6005695518738970761ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 9, key + 7));

  key = 9910494239719928678ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 7, key + 9));

  key = 7947621528143548327;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 9, key + 8));
}

TEST_P(BloomFilterUniform32Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_P(BloomFilterUniform32Test, NoFalseNegativesRangeQuerySmallRange) {
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

TEST_P(BloomFilterUniform32Test, NoFalseNegativesRangeQueryLargeRange) {
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

TEST_P(BloomFilterUniform32Test, NoFalseNegativesRangeQueryExtraLargeRange) {
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

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniform32Test);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniform32Test,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        size_t numKeys = 10000;
        return std::pair<int, BloomFilterRFParameters>{
            numKeys, genParams((rand() % numKeys) + numKeys, 8, 11, 64)};
      });
      return ret;
    }()));
}  // namespace test
}  // namespace filters
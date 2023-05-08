#include <_types/_uint16_t.h>
#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <gtest/gtest.h>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include "gtest/gtest.h"

#include <iomanip>
#include <iterator>
#include <limits>
#include <ostream>
#include <random>
#include <unordered_set>

#include "bloomRF/bloomRF.h"
#include "city/city.h"

using filters::BloomFilterRFParameters;
using filters::BloomRF;

uint64_t randomUniformUint64() {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  return dist(rd);
}

class BloomFilterUniform128Test : public ::testing::Test,
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

  BloomRF<uint64_t, filters::uint128_t> bf{GetParam().second};
};

class BloomFilterUniform32Test : public ::testing::Test {
 protected:
  void SetUp() override {
    for (int i = 0; i < 10000; ++i) {
      auto x = randomUniformUint64();
      s.insert(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::unordered_set<uint64_t> s;

  BloomRF<uint64_t, uint32_t> bf{
      BloomFilterRFParameters{16000, 0, {6, 6, 6, 6, 6, 6}}};
};

TEST_P(BloomFilterUniform128Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_P(BloomFilterUniform128Test, NoFalseNegativesRangeQuerySmallRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10;
    auto high = *it + rand() % 10;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

TEST_P(BloomFilterUniform128Test, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10000;
    auto high = *it + rand() % 10000;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

TEST_P(BloomFilterUniform128Test, NoFalseNegativesRangeQueryExtraLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 100000;
    auto high = *it + rand() % 100000;
    bool found = bf.findRange(low, high);
    if (!found) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

BloomFilterRFParameters genParams() {
  std::random_device rd;
  std::mt19937_64 e2(rd());

  std::uniform_int_distribution<uint64_t> layer(2, 8);

  std::vector<size_t> layers(7);
  std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  return BloomFilterRFParameters{16000, 0, layers};
}

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniform128Test);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniform128Test,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        return std::pair<int, BloomFilterRFParameters>{10000, genParams()};
      });
      return ret;
    }()));

TEST(OneOff, RangeQuery) {
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 0, {2, 2}}};
  uint64_t key = 17183560791176864955ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key + 2));

  key = 6734315744289451875ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key + 1));

  key = 16343179362131379382ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key));

  key = 1894361899248432030ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key + 1));

  key = 994988673032400334ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key - 3, key));

  key = 7947621528143548327;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key - 9, key + 8));
}

#include <_types/_uint16_t.h>
#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <gtest/gtest.h>
#include <algorithm>
#include <cstddef>
#include <cstdlib>

#include <iomanip>
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

class BloomFilterUniform64Test : public ::testing::Test {
 protected:
  void SetUp() override {
    for (int i = 0; i < 10000; ++i) {
      auto x = randomUniformUint64();
      s.push_back(x);
      bf.add(x);
    }
  }
  // Keep track of what has actually been inserted.
  std::vector<uint64_t> s;

  BloomRF<uint64_t, uint64_t> bf{
      BloomFilterRFParameters{16000, 0, {7, 7, 7, 7, 7, 7}}};
};

class BloomFilterUniform128Test : public ::testing::Test {
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

  BloomRF<uint64_t, filters::uint128_t> bf{
      BloomFilterRFParameters{16000, 0, {8, 8, 8, 8, 8, 8}}};
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

class BloomFilterVariableLayers : public ::testing::Test {
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

  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 0, {7, 7, 4, 4, 2, 2}}};
};

TEST_F(BloomFilterUniform64Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST(OneOff, RangeQuery) {
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 0, {7, 7, 7, 7, 7, 7}}};
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
}

TEST_F(BloomFilterUniform64Test, NoFalseNegativesRangeQuerySmallRange) {
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

TEST_F(BloomFilterUniform64Test, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10000;
    auto high = *it + rand() % 10000;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

TEST_F(BloomFilterUniform128Test, NoFalseNegativesRangeQuerySmallRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10;
    auto high = *it + rand() % 10;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

TEST_F(BloomFilterUniform128Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    EXPECT_TRUE(bf.find(*it));
  }
}

TEST(BloomFilterVariableLayersOneOff, BasicRange) {
  auto key = 2978291708368540195ULL;
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 0, {7, 7, 4, 4, 2, 2}}};
  bf.add(key);
  ASSERT_TRUE(bf.findRange(2978291708368540122ULL, 2978291708368543853ULL));

  key = 3641009636154320258ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(3641009636154315818ULL, 3641009636154328423ULL));

  key = 11949888074462238032ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(11949888074462174450ULL, 11949888074462324763ULL));
}

TEST_F(BloomFilterVariableLayers, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_F(BloomFilterVariableLayers, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 100000;
    auto high = *it + rand() % 100000;
    if (!bf.findRange(low, high)) {
      std::ostringstream os;
      os << "Failed lookup. Query range: [" << low << "," << high << "]. ";
      os << "Have key: " << *it << std::endl;
      std::cerr << os.str();
    }
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

TEST_F(BloomFilterUniform128Test, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10000;
    auto high = *it + rand() % 10000;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

TEST_F(BloomFilterUniform32Test, NoFalseNegativesRangeQuerySmallRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10;
    auto high = *it + rand() % 10;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

TEST_F(BloomFilterUniform32Test, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    EXPECT_TRUE(bf.find(*it));
  }
}

TEST_F(BloomFilterUniform32Test, NoFalseNegativesRangeQueryLargeRange) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10000;
    auto high = *it + rand() % 10000;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

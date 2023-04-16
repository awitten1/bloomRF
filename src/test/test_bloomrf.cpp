#include <_types/_uint16_t.h>
#include <_types/_uint64_t.h>
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdlib>

#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <unordered_set>

#include "bloomRF/bloomRF.h"

using filters::BloomFilterRFParameters;
using filters::BloomRF;

uint64_t randomUniformUint64() {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  return dist(rd);
}

class BloomFilterUniformTest : public ::testing::Test {
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

  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 6, 0}};
};

TEST_F(BloomFilterUniformTest, NoFalseNegativesPointQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    EXPECT_TRUE(bf.find(*it));
  }
}

TEST(Basic, OneOffRangeQuery) {
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 6, 0}};
  uint64_t key = 17183560791176864955ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key + 2));

  key = 6734315744289451875ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key + 1));

  key = 16343179362131379382ULL;
  bf.add(key);
  EXPECT_TRUE(bf.findRange(key, key));
}

TEST_F(BloomFilterUniformTest, NoFalseNegativesRangeQuery) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    auto low = *it - rand() % 10;
    auto high = *it + rand() % 10;
    ASSERT_TRUE(bf.findRange(low, high));
  }
}

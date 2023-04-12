#include <gtest/gtest.h>

#include <iomanip>
#include <limits>
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
  std::unordered_set<uint64_t> s;
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 6, 0, 3}};
};

TEST_F(BloomFilterUniformTest, NoFalseNegatives) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

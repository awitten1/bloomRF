#include "bloomRF/bloomRF.h"

#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <random>
#include <unordered_set>

using filters::BloomFilterRFParameters;
using filters::BloomRF;

uint64_t randomUint64() {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());

  return dist(rd);
}

class BloomFilterTest : public ::testing::Test {
 protected:
  void SetUp() override {
    for (int i = 0; i < 10000; ++i) {
      auto x = randomUint64();
      s.insert(x);
      bf.add(x);
    }
  }
  std::unordered_set<uint64_t> s;
  BloomRF<uint64_t> bf{BloomFilterRFParameters{18000, 6, 0, 4}};
};

TEST_F(BloomFilterTest, NoFalseNegatives) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

TEST_F(BloomFilterTest, FalsePositiveRateBound) {
  int false_positives = 0;
  int denominator = 10000;
  int true_positive = 0;
  for (int i = 0; i < denominator; ++i) {
    if (s.find(i) == s.end()) {
      if (bf.find(i)) {
        ++false_positives;
      }
    } else {
      ++true_positive;
    }
  }
  auto fp = static_cast<double>(false_positives) / static_cast<double>(denominator - true_positive);

  ASSERT_LE(fp, .01) << std::setprecision(4) << "False Positive Rate: " << fp;

}

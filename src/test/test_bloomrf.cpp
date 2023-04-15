#include <_types/_uint16_t.h>
#include <_types/_uint64_t.h>
#include <gtest/gtest.h>
#include <algorithm>

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
  std::unordered_set<uint64_t> s;
  BloomRF<uint64_t> bf{BloomFilterRFParameters{16000, 6, 0, 3}};
};

TEST_F(BloomFilterUniformTest, NoFalseNegatives) {
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    ASSERT_TRUE(bf.find(*it));
  }
}

template <typename T>
void printBinary(T t) {
  std::cerr << t << ": ";
  std::vector<int> print;
  while (t > 0) {
    print.push_back(t % 2);
    t >>= 1;
  }
  std::reverse(print.begin(), print.end());
  std::copy(print.begin(), print.end(), std::ostream_iterator<int>(std::cerr));
}

TEST(decomposeDyadicIntervals, Basic) {
  auto x = filters::decomposeIntoDyadicIntervals<uint64_t>(45, 60);
  for (int i = 0; i < x.size(); ++i) {
    std::cerr << i << ": " << std::endl;

    for (auto prefix : x[i]) {
      printBinary(prefix);
      std::cerr << ", ";
    }
    std::cerr << std::endl;
  }
}

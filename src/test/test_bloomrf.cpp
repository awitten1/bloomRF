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

BloomFilterRFParameters genParams(size_t filterSizeBytes, int maxDelta, size_t maxDeltaSum) {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> layer(2, maxDelta);

  std::vector<size_t> layers((rand() % 8) + 2);
  std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  while (std::accumulate(layers.begin(), layers.end(), 0) > maxDeltaSum) {
    std::generate(layers.begin(), layers.end(), [&]() { return layer(rd); });
  }
  return BloomFilterRFParameters{filterSizeBytes, 0, layers};
}

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
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
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
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
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
      os << "Delta vector: ";
      for (auto&& v : bf.getDelta()) {
        os << v << " ";
      }
      std::cerr << os.str();
    }
    ASSERT_TRUE(found);
  }
}

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(BloomFilterUniform128Test);

INSTANTIATE_TEST_SUITE_P(
    NoFalseNegatives,
    BloomFilterUniform128Test,
    testing::ValuesIn([]() {
      std::vector<std::pair<int, BloomFilterRFParameters>> ret;
      std::generate_n(std::back_inserter(ret), 15, []() {
        size_t numKeys = 10000;
        return std::pair<int, BloomFilterRFParameters>{numKeys, genParams((rand() % numKeys) + numKeys, 9, 64)};
      });
      return ret;
    }()));

TEST(OneOff, RangeQuery) {
  BloomRF<uint64_t, uint64_t> bf{
      BloomFilterRFParameters{16000, 0, {9, 8, 6}}};
  uint64_t key = 17183560791176864955ULL;
  bf.add(key);
  ASSERT_TRUE(bf.findRange(key - 100, key + 100));

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
        return std::pair<int, BloomFilterRFParameters>{numKeys, genParams((rand() % numKeys) + numKeys, 9, 64)};
      });
      return ret;
    }()));

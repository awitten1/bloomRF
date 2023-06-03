
#include "experiments.h"
#include "city/city.h"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

namespace {

constexpr uint64_t uint64_max = std::numeric_limits<uint64_t>::max();



template <typename T, typename Generator, typename QueryGenerator>
void runRangeExperiments(T interval_size, Generator d, QueryGenerator qg) {
  static_assert(std::is_same_v<decltype(d()), T>);

  ExperimentDriver<T, Generator, QueryGenerator> ed64U{
    BloomFilterRFParameters{6000000, 0, {10, 8, 4, 4, 4, 4, 4, 3}, {1, 1, 1, 1, 1, 1, 2, 2}}, d, qg};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomRangeQuerys(100000, interval_size);

  std::cout << fp << std::endl;
}

template <typename T, typename Generator, typename QueryGenerator>
void runPointExperiments(Generator d, QueryGenerator qg) {
  static_assert(std::is_same_v<decltype(d()), T>);
  ExperimentDriver<T, Generator, QueryGenerator> ed64U{
    BloomFilterRFParameters{3200000, 0, {8, 8, 6, 6, 5, 5, 4, 3}, {1, 1, 1, 1, 1, 1, 1, 2}}, d, qg};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(100000);

  std::cout << fp << std::endl;
}

}

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());

  // Unsigned integer, normal distribution.
  runPointExperiments<uint64_t>([&]() mutable {
    uint64_t val = std::round(std::normal_distribution<double>(1e8, 1e7)(gen));
    return val;
  },
  [&]() mutable {
    return std::uniform_int_distribution<uint64_t>(0,
      std::numeric_limits<uint64_t>::max())(gen);
  });

  // Unsigned integer, uniform distribution.
  runPointExperiments<uint64_t>([&]() mutable {
    return std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(gen);
  },
  [&]() mutable {
      return std::uniform_int_distribution<uint64_t>(std::numeric_limits<uint64_t>::lowest(),
        std::numeric_limits<uint64_t>::max())(gen);
    }
  );

  std::cout << "-----------RUNNING RANGE EXPERIMENTS-----------" << std::endl;

  // Unsigned integers, normal distribution.
  runRangeExperiments<uint64_t>(1e5, [&]() mutable {
    uint64_t val = std::round(std::normal_distribution<double>(1 << 30, 1 << 28)(gen));
    return val;
  },
  [&]() mutable {
      return std::uniform_int_distribution<uint64_t>(0,
        std::numeric_limits<uint64_t>::max())(gen);
    }
  );

  // Unsigned integers, uniform distribution.
  runRangeExperiments<uint64_t>(1e8, [&]() mutable {
    return std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(gen);
  },
  [&]() mutable {
      return std::uniform_int_distribution<uint64_t>(std::numeric_limits<uint64_t>::lowest(),
        std::numeric_limits<uint64_t>::max())(gen);
    });

  // Floats, uniform distribution.
  runRangeExperiments<double>(100, [&]() mutable {
    return std::uniform_real_distribution<double>(0, std::numeric_limits<double>::max())(gen);
  },
  [&]() mutable {
      return std::uniform_real_distribution<double>(0,
        std::numeric_limits<double>::max())(gen);
    });
}

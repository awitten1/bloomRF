
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



template <typename T, typename Generator>
void runRangeExperiments(T interval_size, Generator d) {
  static_assert(std::is_same_v<decltype(d()), T>);

  ExperimentDriver<T, Generator> ed64U{
    BloomFilterRFParameters{2600000, 0, {10, 8, 6, 6, 5, 5, 4, 3}}, d};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomRangeQuerys(100000, interval_size);

  std::cout << fp << std::endl;
}

template <typename T, typename Generator>
void runPointExperiments(Generator d) {
  static_assert(std::is_same_v<decltype(d()), T>);
  ExperimentDriver<T, Generator> ed64U{
    BloomFilterRFParameters{3200000, 0, {8, 8, 6, 6, 5, 5, 4, 3}}, d};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(100000);

  std::cout << fp << std::endl;
}

}

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());
  runPointExperiments<uint64_t>([=]() mutable {
    auto val = std::normal_distribution<double>(1 << 30, 1 << 30)(gen);
    return static_cast<uint64_t>(val);
  });
  runPointExperiments<uint64_t>([=]() mutable {
    return std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(gen);
  });
  std::cout << "-----------RUNNING RANGE EXPERIMENTS-----------" << std::endl;
  runRangeExperiments<uint64_t>(1e7, [=]() mutable {
    auto val = std::normal_distribution<double>(1 << 30, 1 << 10)(gen);
    return static_cast<uint64_t>(val);
  });
  runRangeExperiments<uint64_t>(1e7, [=]() mutable {
    return std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(gen);
  });

}

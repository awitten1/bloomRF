
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
void runRangeExperiments(T interval_size, Generator d, std::string msg) {
  static_assert(std::is_same_v<decltype(d()), T>);

  ExperimentDriver<T, Generator> ed64U{
      BloomFilterRFParameters{6000000, 0, {7, 7, 7, 4, 4, 2, 2}}, d};

  std::cout << "Running experiment: " << msg << std::endl;

  ed64U.doInserts(2000000);
  double fp = ed64U.randomRangeQuerys(100000, interval_size);

  std::cout << fp << std::endl;
}

template <typename T, typename Generator>
void runPointExperiments(Generator d, std::string msg) {
  static_assert(std::is_same_v<decltype(d()), T>);
  ExperimentDriver<T, Generator> ed64U{
      BloomFilterRFParameters{3200000, 0, {8, 8, 6, 6, 5, 5, 4, 3}}, d};

  std::cout << "Running experiment: " << msg << std::endl;

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(100000);

  std::cout << fp << std::endl;
}

}  // namespace

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());
  runPointExperiments<uint64_t>(
      [&]() mutable {
        auto val =
            std::normal_distribution<double>(1ULL << 30, 1ULL << 29)(gen);
        return static_cast<uint64_t>(val);
      },
      "point query, unsigned integer, normal distribution");
  runPointExperiments<uint64_t>(
      [&]() mutable {
        return std::uniform_int_distribution<uint64_t>(
            0, std::numeric_limits<uint64_t>::max())(gen);
      },
      "point query, unsigned integer, uniform distribution");
  std::cout << "-----------RUNNING RANGE EXPERIMENTS-----------" << std::endl;
  runRangeExperiments<uint64_t>(
      1e9,
      [&]() mutable {
        uint64_t val = std::round(
            std::normal_distribution<double>(1ULL << 30, 1ULL << 29)(gen));
        return val;
      },
      "unsigned integers, normal distriubtion");
  runRangeExperiments<uint64_t>(
      1e8,
      [&]() mutable {
        return std::uniform_int_distribution<uint64_t>(
            0, std::numeric_limits<uint64_t>::max())(gen);
      },
      "unsigned integers, uniform distribution");
  runRangeExperiments<double>(
      10,
      [&]() mutable {
        return std::uniform_real_distribution<double>(
            0, std::numeric_limits<double>::max())(gen);
      },
      "floats, uniform distribution");
  runRangeExperiments<double>(
      1,
      [&]() mutable {
        return std::normal_distribution<double>(1ULL << 30, 1ULL << 29)(gen);
      },
      "floats, normal distribution");
}


#include "experiments.h"
#include "city/city.h"

#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <random>

template <typename T, typename UnderType>
void runExperimentsForUniform() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  auto genUniform = [=, uniform = std::uniform_int_distribution<T>{
                            0}]() mutable { return uniform(mt); };
  ExperimentDriver<T, decltype(genUniform), UnderType> ed64U{
      BloomFilterRFParameters{3200000, 0, {6,6,6,6,6}}, genUniform};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(30000);

  std::cout << fp << std::endl;
}

template <typename T, typename UnderType>
void runRangeExperimentsForUniform() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  auto genUniform = [=, uniform = std::uniform_int_distribution<T>{
                            0, std::numeric_limits<T>::max() >> 1}]() mutable {
    return uniform(mt);
  };
  ExperimentDriver<T, decltype(genUniform), UnderType> ed64U{
      BloomFilterRFParameters{3200000, 0, {6,6,6,6,6,6}}, genUniform};

  ed64U.doInserts(2000000);
  double fp =
      ed64U.randomRangeQuerys(30000,
                              [=, uniform = std::uniform_int_distribution<T>{
                                      (std::numeric_limits<T>::max() >> 1) + 1,
                                      std::numeric_limits<T>::max() -
                                          10000}]() mutable -> std::pair<T, T> {
                                auto start = uniform(mt);
                                return {start, start + 5000};
                              });

  std::cout << fp << std::endl;
}

template <typename T, typename UnderType>
void runExperimentsForNormal() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  std::normal_distribution<> normal{1ULL << 31, 1ULL << 29};
  auto genNormal = [&]() mutable {
    auto ret = normal(mt);
    return ret;
  };

  ExperimentDriver<T, decltype(genNormal), UnderType> ed64N{
      BloomFilterRFParameters{3200000, 0, {6,6,6,6,6,6}}, genNormal};

  ed64N.doInserts(2000000);

  auto fp = ed64N.randomQuerys(30000);

  std::cout << fp << std::endl;
}

int main() {
  runExperimentsForUniform<uint64_t, uint32_t>();
  runExperimentsForNormal<uint64_t, uint32_t>();
  runExperimentsForUniform<uint64_t, uint64_t>();
  runExperimentsForNormal<uint64_t, uint64_t>();
  //runExperimentsForUniform<uint64_t, filters::uint128_t>();
  //runExperimentsForNormal<uint64_t, filters::uint128_t>();
  std::cout << "-----------RUNNING RANGE EXPERIMENTS-----------" << std::endl;
  runRangeExperimentsForUniform<uint64_t, uint32_t>();
  runRangeExperimentsForUniform<uint64_t, uint64_t>();
  //runRangeExperimentsForUniform<uint64_t, filters::uint128_t>();
}

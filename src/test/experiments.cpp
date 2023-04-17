
#include "experiments.h"

#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <cstddef>
#include <iomanip>
#include <random>

template<typename T, typename UnderType, size_t Delta>
void runExperimentsForUniform() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  auto genUniform = [=, uniform = std::uniform_int_distribution<>{
                            0}]() mutable { return uniform(mt); };
  ExperimentDriver<T, decltype(genUniform), UnderType, Delta> ed64U{
      BloomFilterRFParameters{3200000, 6, 0}, genUniform};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(30000);

  std::cout << fp << std::endl;
}


template<typename T, typename UnderType, size_t Delta>
void runExperimentsForNormal() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  std::normal_distribution<> normal{1ULL << 31, 1ULL << 29};
  auto genNormal = [&]() mutable {
    auto ret = normal(mt);
    return ret;
  };

  ExperimentDriver<T, decltype(genNormal), UnderType, Delta> ed64N{
      BloomFilterRFParameters{3200000, 6, 0}, genNormal};

  ed64N.doInserts(2000000);

  auto fp = ed64N.randomQuerys(30000);

  std::cout << fp << std::endl;
}


int main() {
  runExperimentsForUniform<uint64_t, uint32_t, 6>();
  runExperimentsForNormal<uint64_t, uint32_t, 6>();
  runExperimentsForUniform<uint64_t, uint64_t, 7>();
  runExperimentsForNormal<uint64_t, uint64_t, 7>();
  runExperimentsForUniform<uint64_t, filters::uint128_t, 8>();
  runExperimentsForNormal<uint64_t, filters::uint128_t, 8>();
}

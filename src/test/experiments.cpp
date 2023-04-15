
#include "experiments.h"

#include <iomanip>
#include <random>

void runExperimentsForUniform() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  auto genUniform = [=, uniform = std::uniform_int_distribution<>{
                            0}]() mutable { return uniform(mt); };
  ExperimentDriver<uint64_t, decltype(genUniform)> ed64U{
      BloomFilterRFParameters{4000000, 11, 0}, genUniform};

  ed64U.doInserts(2000000);
  double fp = ed64U.randomQuerys(30000);

  std::cout << fp << std::endl;
}

void runExperimentsForNormal() {
  std::random_device rd;
  std::mt19937 mt{rd()};
  std::normal_distribution<> normal{1ULL << 31, 1ULL << 29};
  auto genNormal = [&]() mutable {
    auto ret = normal(mt);
    return ret;
  };

  ExperimentDriver<uint64_t, decltype(genNormal)> ed64N{
      BloomFilterRFParameters{3000000, 11, 0}, genNormal};

  ed64N.doInserts(2000000);

  auto fp = ed64N.randomQuerys(30000);

  std::cout << fp << std::endl;
}

int main() {
  runExperimentsForUniform();
  runExperimentsForNormal();
}

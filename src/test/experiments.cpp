
#include "experiments.h"
#include "city/city.h"

#include <chrono>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

namespace {

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

template <typename T>
void runRangeExperiments(T interval_size,
                         std::function<T()> d,
                         std::function<T()> qkg,
                         const std::string& msg) {
  ExperimentDriver<T> ed64U{
      BloomFilterRFParameters{4000000, 0, {7, 7, 7, 4, 4, 2, 2, 2}}, d, qkg};
  std::cout << "------------------------" << std::endl;

  std::cout << "Running experiment: " << msg << std::endl;

  ed64U.doInserts(2000000);

  auto t1 = high_resolution_clock::now();
  double fp = ed64U.randomRangeQuerys(100000, interval_size);
  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;

  std::cout << "time for 100000 range queries: " << ms_double.count() << "ms\n";
  std::cout << fp << std::endl;
  std::cout << "------------------------" << std::endl;
}

template <typename T>
void runPointExperiments(std::function<T()> d,
                         std::function<T()> qkg,
                         const std::string& msg) {
  std::cout << "------------------------" << std::endl;
  static_assert(std::is_same_v<decltype(d()), T>);
  ExperimentDriver<T> ed64U{
      BloomFilterRFParameters{3200000, 0, {8, 8, 6, 6, 5, 5, 4, 3}}, d, qkg};

  std::cout << "Running experiment: " << msg << std::endl;

  ed64U.doInserts(2000000);

  auto t1 = high_resolution_clock::now();
  double fp = ed64U.randomQuerys(100000);
  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;

  std::cout << "time for 100000 point queries: " << ms_double.count() << "ms\n";

  std::cout << fp << std::endl;
  std::cout << "------------------------" << std::endl;
}

}  // namespace

namespace {

std::random_device rd;
std::mt19937 gen(rd());

uint64_t genNormalUInt(uint64_t mean, uint64_t stddev) {
  uint64_t val =
      std::round(std::normal_distribution<double>(mean, stddev)(gen));
  // std::cout << val << std::endl;
  return val;
}

uint64_t genUniformUInt(uint64_t low, uint64_t high) {
  uint64_t val = std::round(  // NOLINT
      std::uniform_int_distribution<uint64_t>(low, high)(gen));
  // std::cout << val << std::endl;
  return val;
}

double genUniformDouble(double low, double high) {
  return std::uniform_real_distribution<double>(low, high)(gen);
}

double genNormalDouble(double mean, double stddev) {
  return std::normal_distribution<double>(mean, stddev)(gen);
}

}  // namespace

int main() {
  runPointExperiments<uint64_t>(
      []() { return genNormalUInt(1ULL << 33, 1ULL << 31); },
      []() { return genNormalUInt(1ULL << 33, 1ULL << 31); },
      "point query, unsigned integer, normal distribution");

  runPointExperiments<uint64_t>(
      []() { return genUniformUInt(0, std::numeric_limits<uint64_t>::max()); },
      []() { return genUniformUInt(0, std::numeric_limits<uint64_t>::max()); },
      "point query, unsigned integer, uniform distribution");

  runRangeExperiments<uint64_t>(
      1e8, []() { return genNormalUInt(1ULL << 33, 1ULL << 31); },
      []() { return genUniformUInt(0, std::numeric_limits<uint64_t>::max()); },
      "unsigned integers, normal distriubtion, uniform query "
      "keys.");

  // This experiment has severely degraded performance (a completely untenable
  // false positive rate).
  runRangeExperiments<uint64_t>(
      1e5, []() { return genNormalUInt(1ULL << 33, 1ULL << 31); },
      []() { return genNormalUInt(1ULL << 33, 1ULL << 31); },
      "unsigned integers, normal distriubtion, normal query "
      "keys.");

  runRangeExperiments<uint64_t>(
      1e8,
      []() { return genUniformUInt(0, std::numeric_limits<uint64_t>::max()); },
      []() { return genUniformUInt(0, std::numeric_limits<uint64_t>::max()); },
      "unsigned integers, uniform distribution for both inserts and queries.");

  runRangeExperiments<double>(
      10,
      []() { return genUniformDouble(0, std::numeric_limits<double>::max()); },
      []() { return genUniformDouble(0, std::numeric_limits<double>::max()); },
      "floats, uniform distribution");

  runRangeExperiments<double>(
      1, []() { return genNormalDouble(1ULL << 32, 1ULL << 31); },
      []() { return genNormalDouble(1ULL << 32, 1ULL << 31); },
      "floats, normal distribution");
}

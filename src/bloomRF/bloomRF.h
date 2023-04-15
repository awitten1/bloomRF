#include <_types/_uint16_t.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include "city/city.h"

namespace filters {

struct BloomFilterRFParameters {
  BloomFilterRFParameters(size_t filter_size_,
                          size_t filter_hashes_,
                          size_t seed_,
                          uint16_t delta_);

  /// size of filter in bytes.
  size_t filter_size;
  /// number of used hash functions.
  size_t filter_hashes;
  /// random seed for hash functions generation.
  size_t seed;
  /// Delta: distance between layers.
  uint16_t delta;
};

template <typename T>
class BloomRF {
 public:
  using UnderType = uint64_t;
  using Container = std::vector<UnderType>;

  explicit BloomRF(const BloomFilterRFParameters& params);

  void add(T data);

  bool find(T data);

  bool findRange(T low, T high);

  const Container& getFilter() const { return filter; }
  Container& getFilter() { return filter; }

 private:

  class Checks {
  public:
    Checks(T lkey_, T hkey_, uint8_t delta_, size_t num_hashes_) : lkey(lkey_),
          hkey(hkey_), delta(delta_), num_hashes(num_hashes_) {}
    const std::vector<std::pair<T, T>>& getChecks() {
      return checks;
    }

    void init_checks();

  private:
    struct Check {
      Check(T low_, T high_, bool is_covering_, bool on_left_) :
        low{low_}, high{high_}, is_covering{is_covering_}, on_left{on_left_} {}
      T low;
      T high;
      bool is_covering;
      bool on_left; // is this check on the left?
    };
    std::vector<Check> checks;
    size_t domain_width = 8 * sizeof(T);
    T lkey;
    T hkey;
    uint16_t delta;
    size_t num_hashes;
  };

  explicit BloomRF(size_t size_, size_t hashes_, size_t seed_, uint16_t delta_);

  /// Returns size in bits.
  size_t numBits() { return 8 * sizeof(UnderType) * filter.size(); }

  /// Computes the ith PMHF hash of data in terms of ith generic hash.
  size_t bloomRFHash(T data, size_t i);

  size_t hash(T data, size_t i);


  /// Number of hash functions.
  size_t hashes;

  /// Seed of hash functions.
  size_t seed;


  size_t words;

  /// Distance between layers.
  uint16_t delta;

  /// Size of the domain in bits.
  uint16_t domain_size = 8 * sizeof(T);

  Container filter;
};

/// Decompose the range [low, high] into dyadic intervals.
template<typename T>
std::vector<std::vector<T>> decomposeIntoDyadicIntervals(T low, T high);

template <typename T>
size_t bloomRFSize(const BloomRF<T>& bf) {
  std::cout << sizeof(T) << std::endl;
  std::cout << bf.getFilter().size() << std::endl;

  return (bf.getFilter().size() * sizeof(T));
}

}  // namespace filters

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
  /// size -- size of filter in bytes.
  /// hashes -- number of used hash functions.
  /// seed -- random seed for hash functions generation.
  /// delta -- distance between bloomRF layers.
  explicit BloomRF(size_t size_, size_t hashes_, size_t seed_, uint16_t delta_);

  /// Returns size in bits.
  size_t numBits() { return 8 * sizeof(UnderType) * filter.size(); }

  size_t bloomRFHash(T data, size_t i);

  size_t hash(T data, size_t i);

  size_t hashes;
  size_t seed;
  size_t words;
  uint16_t delta;

  Container filter;
};

template <typename T>
size_t bloomRFSize(const BloomRF<T>& bf) {
  std::cout << sizeof(T) << std::endl;
  std::cout << bf.getFilter().size() << std::endl;

  return (bf.getFilter().size() * sizeof(T));
}

}  // namespace filters

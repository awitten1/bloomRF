#include <cstddef>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

#include "city/city.h"

namespace filters {

struct BloomFilterRFParameters {
  BloomFilterRFParameters(size_t filter_size_,
                          size_t seed_,
                          std::vector<size_t> delta_);

  /// size of filter in bytes.
  size_t filter_size;
  /// random seed for hash functions generation.
  size_t seed;
  /// Distance between layers.
  std::vector<size_t> delta;
};

#ifdef __SIZEOF_INT128__
using uint128_t = __uint128_t;
#endif

template <typename T, typename UnderType = uint64_t>
class BloomRfImpl {
 public:
  using Container = std::vector<UnderType>;

  explicit BloomRfImpl(const BloomFilterRFParameters& params);

  void add(T data);

  bool find(T data) const;

  bool findRange(T low, T high) const;

  const Container& getFilter() const { return filter; }
  Container& getFilter() { return filter; }

  const std::vector<size_t>& getDelta() { return delta; }

 private:
  class Checks {
   public:
    enum class IntervalLocation { Left, Right, NotYetSplit };

    struct Check {
      Check(T low_, T high_, IntervalLocation loc_)
          : low{low_}, high{high_}, loc{loc_} {}

      T low;
      T high;
      IntervalLocation loc;
    };

    Checks(T lkey_, T hkey_, std::vector<Check>&& checks_)
        : checks{std::move(checks_)}, lkey(lkey_), hkey(hkey_) {}

    const auto& getChecks() { return checks; }

    void initChecks(size_t delta_sum, size_t delta_back);

    void advanceChecks(size_t times);

    void concatenateChecks(const Checks& other) {
      checks.insert(checks.end(), other.checks.begin(), other.checks.end());
    }

    void compressChecks(size_t total_shift);

   private:
    std::vector<Check> checks;
    size_t domain_width = 8 * sizeof(T);
    T lkey;
    T hkey;

    friend class BloomRfImpl;
  };

  UnderType buildBitMaskForRange(T low, T high, size_t i, int wordPos) const;

  explicit BloomRfImpl(size_t size_, size_t seed_, std::vector<size_t> delta);

  /// Returns size in bits.
  size_t numBits() const { return 8 * sizeof(UnderType) * filter.size(); }

  /// Computes the ith PMHF hash of data. Only returns the word
  /// to which the data maps to.  Use bloomRFRemainder to retrieve the
  /// offset.
  size_t bloomRFHashToWord(T data, size_t i) const;

  UnderType bloomRFRemainder(T data, size_t i, int wordPos) const;

  std::pair<size_t, UnderType> hashToIndexAndBitMask(T data, size_t i) const;

  size_t hash(T data, size_t i) const;

  bool checkDIOfDecomposition(T low, T high, int layer) const;

  /// Number of hash functions.
  size_t hashes;

  /// Seed of hash functions.
  size_t seed;

  size_t words;

  /// Distance between layers.
  std::vector<size_t> delta;

  /// Prefix sums of delta.
  std::vector<size_t> shifts;

  /// Size of the domain in bits.
  uint16_t domain_size = 8 * sizeof(T);

  Container filter;
};

template<typename Key, typename UnderType = uint64_t,
    typename = void>
class BloomRF : private BloomRfImpl<Key, UnderType> {
public:
  using BloomRfImpl<Key, UnderType>::add;
  using BloomRfImpl<Key, UnderType>::find;
  using BloomRfImpl<Key, UnderType>::findRange;
  using BloomRfImpl<Key, UnderType>::getDelta;
  using BloomRfImpl<Key, UnderType>::getFilter;
  using BloomRfImpl<Key, UnderType>::BloomRfImpl;
};

template<typename Key, typename UnderType>
class BloomRF<Key, UnderType, std::enable_if_t<std::is_signed_v<Key>>>
  : private BloomRfImpl<std::make_unsigned_t<Key>, UnderType> {

  using UnsignedKey = std::make_unsigned_t<Key>;

public:
  void add(Key data) {
    UnsignedKey unsignedData = static_cast<UnsignedKey>(data) - std::numeric_limits<Key>::min();
    BloomRfImpl<UnsignedKey, UnderType>::add(unsignedData);
  }

  bool find(Key data) {
    UnsignedKey unsignedData = static_cast<UnsignedKey>(data) - std::numeric_limits<Key>::min();
    return BloomRfImpl<UnsignedKey, UnderType>::find(unsignedData);
  }

  bool findRange(Key low, Key high) {
    UnsignedKey unsignedLow = static_cast<UnsignedKey>(low) - std::numeric_limits<Key>::min();
    UnsignedKey unsignedHigh = static_cast<UnsignedKey>(high) - std::numeric_limits<Key>::min();
    return BloomRfImpl<UnsignedKey, UnderType>::findRange(unsignedLow, unsignedHigh);
  }

  using BloomRfImpl<UnsignedKey, UnderType>::getDelta;
  using BloomRfImpl<UnsignedKey, UnderType>::getFilter;
  using BloomRfImpl<UnsignedKey, UnderType>::BloomRfImpl;

};



}  // namespace filters

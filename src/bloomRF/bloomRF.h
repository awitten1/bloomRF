#include <cstddef>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>
#include <bit>
#include <cstring>

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

namespace detail {

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

}

//
// A wrapper on BloomRFImpl.
//
template <typename Key, typename UnderType = uint64_t, typename = void>
class BloomRF : private detail::BloomRfImpl<Key, UnderType> {
 public:
  using detail::BloomRfImpl<Key, UnderType>::add;
  using detail::BloomRfImpl<Key, UnderType>::find;
  using detail::BloomRfImpl<Key, UnderType>::findRange;
  using detail::BloomRfImpl<Key, UnderType>::getDelta;
  using detail::BloomRfImpl<Key, UnderType>::getFilter;
  using detail::BloomRfImpl<Key, UnderType>::BloomRfImpl;
};


//
// BloomRF for signed integers is implmented in terms of BloomRF
// for unsigned integers.  It first converts signed integers to unsigned
// integers in an order preserving way, and then makes the calls to
// the BloomRF on the unsigned integers.
//
template <typename Key, typename UnderType>
class BloomRF<Key, UnderType, std::enable_if_t<std::is_signed_v<Key> && !std::is_floating_point_v<Key>>>
    : private detail::BloomRfImpl<std::make_unsigned_t<Key>, UnderType> {
  using UnsignedKey = std::make_unsigned_t<Key>;

 public:
  void add(Key data) {
    UnsignedKey unsignedData =
        static_cast<UnsignedKey>(data) - static_cast<UnsignedKey>(std::numeric_limits<Key>::min());
    detail::BloomRfImpl<UnsignedKey, UnderType>::add(unsignedData);
  }

  bool find(Key data) {
    UnsignedKey unsignedData =
        static_cast<UnsignedKey>(data) - static_cast<UnsignedKey>(std::numeric_limits<Key>::min());
    return detail::BloomRfImpl<UnsignedKey, UnderType>::find(unsignedData);
  }

  bool findRange(Key low, Key high) {
    UnsignedKey unsignedLow =
        static_cast<UnsignedKey>(low) - static_cast<UnsignedKey>(std::numeric_limits<Key>::min());
    UnsignedKey unsignedHigh =
        static_cast<UnsignedKey>(high) - static_cast<UnsignedKey>(std::numeric_limits<Key>::min());
    return detail::BloomRfImpl<UnsignedKey, UnderType>::findRange(unsignedLow,
                                                          unsignedHigh);
  }

  using detail::BloomRfImpl<UnsignedKey, UnderType>::getDelta;
  using detail::BloomRfImpl<UnsignedKey, UnderType>::getFilter;
  using detail::BloomRfImpl<UnsignedKey, UnderType>::BloomRfImpl;
};


namespace detail {

template<typename FloatType>
struct FloatToUInt;

template<>
struct FloatToUInt<float> {
  using uint_type = uint32_t;
};

template<>
struct FloatToUInt<double> {
  using uint_type = uint64_t;
};

} // detail



//
// Floating point support.  We need to map floats to unsigned integers in a way that preserves
// ordering (if x < y then converted(x) < converted(y)).
//
template <typename FloatKey, typename UnderType>
class BloomRF<FloatKey, UnderType, std::enable_if_t<std::is_floating_point_v<FloatKey>>>
    : private detail::BloomRfImpl<typename detail::FloatToUInt<FloatKey>::uint_type, UnderType> {

  using UnsignedKey = typename detail::FloatToUInt<FloatKey>::uint_type;
  using SignedKey = std::make_signed_t<UnsignedKey>;

  static UnsignedKey order_preserving_repr(FloatKey x)
  {
      static_assert(sizeof(FloatKey) == sizeof(UnsignedKey));
      static_assert(std::numeric_limits<FloatKey>::is_iec559);

      // floats and ints must have the same endianness.  On platforms for which the endianness
      // of floats and ints are different, std::endian::native should be neither big, nor little.
      // Therefore, here we assert that the endianness is either big or little.
      static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little);

      SignedKey k;
      std::memcpy(&k, &x, sizeof k);
      if (k<0) k ^= std::numeric_limits<SignedKey>::max();
      return static_cast<UnsignedKey>(k) - static_cast<UnsignedKey>(std::numeric_limits<SignedKey>::min());
  }

public:
  void add(FloatKey data) {
    UnsignedKey unsignedData = order_preserving_repr(data);
    detail::BloomRfImpl<UnsignedKey, UnderType>::add(unsignedData);
  }

  bool find(FloatKey data) {
    UnsignedKey unsignedData = order_preserving_repr(data);
    return detail::BloomRfImpl<UnsignedKey, UnderType>::find(unsignedData);
  }

  bool findRange(FloatKey low, FloatKey high) {
    UnsignedKey unsignedLow = order_preserving_repr(low);
    UnsignedKey unsignedHigh = order_preserving_repr(high);
    return detail::BloomRfImpl<UnsignedKey, UnderType>::findRange(unsignedLow,
                                                          unsignedHigh);
  }

  using detail::BloomRfImpl<UnsignedKey, UnderType>::getDelta;
  using detail::BloomRfImpl<UnsignedKey, UnderType>::getFilter;
  using detail::BloomRfImpl<UnsignedKey, UnderType>::BloomRfImpl;

};

}  // namespace filters

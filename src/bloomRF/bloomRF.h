#include <_types/_uint16_t.h>
#include <_types/_uint64_t.h>
#include <cstddef>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

#include "city/city.h"

namespace filters {

struct BloomFilterRFParameters {
  BloomFilterRFParameters(size_t filter_size_,
                          size_t filter_hashes_,
                          size_t seed_);

  /// size of filter in bytes.
  size_t filter_size;
  /// number of used hash functions.
  size_t filter_hashes;
  /// random seed for hash functions generation.
  size_t seed;
};

#ifdef __SIZEOF_INT128__
namespace {
using uint128_t = __uint128_t;
}
#endif

template <typename T, typename UnderType = uint64_t, size_t Delta = 7>
class BloomRF {
  // A current limitation of the implementation. This should eventually
  // support other UnderTypes (in particular uint128_t).  The implementation
  // requires that a "bloomRF word" of width 2^(delta - 1) be the same size as
  // the UnderType.
  static_assert(std::is_integral_v<UnderType> && std::is_unsigned_v<UnderType>,
                "Must be integral and unsigned.");

  // This check follows from the above requriement; that a bloomRF word be the
  // same size as an UnderType.
  static_assert(8 * sizeof(UnderType) == 1 << (Delta - 1),
                "Delta is inconsistent with the size of UnderType.");

 public:
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
    enum class IntervalLocation { Left, Right, NotYetSplit };

    struct Check {
      Check(T low_, T high_, bool is_covering_, IntervalLocation loc_)
          : low{low_}, high{high_}, is_covering{is_covering_}, loc{loc_} {}

      T low;
      T high;
      bool is_covering;
      IntervalLocation loc;
    };

    Checks(T lkey_, T hkey_, std::vector<Check>&& checks_)
        : checks{std::move(checks_)}, lkey(lkey_), hkey(hkey_) {}

    const auto& getChecks() { return checks; }

    void initChecks(size_t num_hashes, uint16_t delta);

    void advanceChecks(size_t times);

    void concatenateChecks(const Checks& other) {
      checks.insert(checks.end(), other.checks.begin(), other.checks.end());
    }

   private:
    std::vector<Check> checks;
    size_t domain_width = 8 * sizeof(T);
    T lkey;
    T hkey;

    friend class BloomRF;
  };

  UnderType buildBitMaskForRange(T low, T high, size_t i) {
    UnderType lowOffset = 1 << ((low >> (i * Delta)) & remainderMask);
    UnderType highOffset = 1 << ((high >> (i * Delta)) & remainderMask);

    UnderType ret = ~0;
    ret &= ~((1 << lowOffset) - 1);
    ret &= ((1 << (highOffset + 1)) - 1);

    return ~ret;
  }

  explicit BloomRF(size_t size_, size_t hashes_, size_t seed_);

  /// Returns size in bits.
  size_t numBits() { return 8 * sizeof(UnderType) * filter.size(); }

  /// Computes the ith PMHF hash of data in terms of ith generic hash.
  /// Only returns the word that the data item maps to.
  /// Does not return offset.
  size_t bloomRFHashToWord(T data, size_t i);

  UnderType bloomRFRemainder(T data, size_t i);

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

  inline static constexpr UnderType remainderMask = (1 << Delta) - 1;

  Container filter;
};

}  // namespace filters

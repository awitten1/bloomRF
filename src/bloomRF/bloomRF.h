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
class BloomRF {
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
      Check(T low_, T high_, IntervalLocation loc_)
          : low{low_}, high{high_}, loc{loc_} {}

      T low;
      T high;
      IntervalLocation loc;
    };

    Checks(T lkey_, T hkey_, std::vector<Check>&& checks_)
        : checks{std::move(checks_)}, lkey(lkey_), hkey(hkey_) {}

    const auto& getChecks() { return checks; }

    void initChecks(size_t delta_sum);

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

  UnderType buildBitMaskForRange(T low, T high, size_t i, int wordPos);

  explicit BloomRF(size_t size_, size_t seed_, std::vector<size_t> delta);

  /// Returns size in bits.
  size_t numBits() { return 8 * sizeof(UnderType) * filter.size(); }

  /// Computes the ith PMHF hash of data. Only returns the word
  /// to which the data maps to.  Use bloomRFRemainder to retrieve the
  /// offset.
  size_t bloomRFHashToWord(T data, size_t i);

  UnderType bloomRFRemainder(T data, size_t i, int wordPos);

  std::ldiv_t getFilterPosAndOffset(size_t pos, size_t i);

  size_t hash(T data, size_t i);

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

}  // namespace filters

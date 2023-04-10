#include "bloomRF.h"

namespace filters {

namespace {

constexpr uint64_t SEED_GEN_A = 845897321;
constexpr uint64_t SEED_GEN_B = 217728422;

constexpr uint64_t MAX_BLOOM_FILTER_SIZE = 1 << 30;

} // namespace

BloomFilterRFParameters::BloomFilterRFParameters(size_t filter_size_,
                                                 size_t filter_hashes_,
                                                 size_t seed_, uint16_t delta_)
    : filter_size(filter_size_), filter_hashes(filter_hashes_), seed(seed_),
      delta(delta_) {
  if (filter_size == 0)
    throw "The size of bloom filter cannot be zero";
  if (filter_hashes == 0)
    throw "The number of hash functions for bloom filter cannot be zero";
}

template <typename T>
BloomRF<T>::BloomRF(const BloomFilterRFParameters &params)
    : BloomRF<T>(params.filter_size, params.filter_hashes, params.seed,
                 params.delta) {}

template <typename T> size_t BloomRF<T>::bloomRFHash(T data, size_t i) {

  auto hash = this->hash(data >> (delta - 1), i);
  return (hash % (numBits() >> (delta - 1)))
         << (delta - 1) + (data & (1 << (delta - 1)) - 1);
}

template <typename T> size_t BloomRF<T>::hash(T data, size_t i) {

    size_t hash1 = CityHash64WithSeed(reinterpret_cast<const char *>(&data),
                                      sizeof(data), seed);
    size_t hash2 =
        CityHash64WithSeed(reinterpret_cast<const char *>(&data), sizeof(data),
                           SEED_GEN_A * seed + SEED_GEN_B);
    return hash1 + i * hash2 + i * i;
}

template <typename T> void BloomRF<T>::add(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = hash(data, i) % numBits();
    filter[pos / (8 * sizeof(UnderType))] |=
        (1ULL << (pos % (8 * sizeof(UnderType))));
    data >>= delta;
  }
}

template <typename T> bool BloomRF<T>::find(T data) {
  for (size_t i = 0; i < hashes; ++i) {
    size_t pos = hash(data, i) % numBits();
    if (!(filter[pos / (8 * sizeof(UnderType))] &
          (1ULL << (pos % (8 * sizeof(UnderType)))))) {
      return false;
    }
    data >>= delta;
  }
  return true;
}

template <typename T>
BloomRF<T>::BloomRF(size_t size_, size_t hashes_, size_t seed_, uint16_t delta_)
    : hashes(hashes_), seed(seed_),
      words((size_ + sizeof(UnderType) - 1) / sizeof(UnderType)), delta(delta_),
      filter(words, 0) {}

template class BloomRF<uint16_t>;
template class BloomRF<uint32_t>;
template class BloomRF<uint64_t>;

} // namespace filters

#include "bloomRF.h"

namespace filters {


namespace {

constexpr uint64_t SEED_GEN_A = 845897321;
constexpr uint64_t SEED_GEN_B = 217728422;

constexpr uint64_t MAX_BLOOM_FILTER_SIZE = 1 << 30;

}

BloomFilterRFParameters::BloomFilterRFParameters(size_t filter_size_, size_t filter_hashes_, size_t seed_, uint16_t delta_)
    : filter_size(filter_size_), filter_hashes(filter_hashes_), seed(seed_), delta(delta_)
{
    if (filter_size == 0)
        throw "The size of bloom filter cannot be zero";
    if (filter_hashes == 0)
        throw "The number of hash functions for bloom filter cannot be zero";
}

template<typename T>
BloomRF<T>::BloomRF(const BloomFilterRFParameters& params) :
    BloomRF<T>(params.filter_size, params.filter_hashes, params.seed, params.delta) {}


template<typename T>
size_t BloomRF<T>::bloomRFHash(T data, size_t i) {
    size_t hash1 = CityHash64WithSeed(data >> (delta - 1), sizeof(data), seed);
    size_t hash2 = CityHash64WithSeed(data >> (delta - 1), sizeof(data), SEED_GEN_A * seed + SEED_GEN_B);
    return (hash1 + i * hash2 + i * i);
}

template<typename T>
BloomRF<T>::BloomRF(size_t size_, size_t hashes_, size_t seed_, uint16_t delta_)
  : size(size_),
    hashes(hashes_),
    seed(seed_),
    words((size + sizeof(UnderType) - 1) / sizeof(UnderType)),
    delta(delta_),
    filter(words, 0) {}





} // namespace filters


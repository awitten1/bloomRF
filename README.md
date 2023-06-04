# bloomRF

This is an implementation of the EDBT 2023 best paper [bloomRF](https://openproceedings.org/2023/conf/edbt/paper-190.pdf).
BloomRF is a probabilistic filter that supports range queries, in addition to point queries.  Traditional
[bloom filters](https://en.wikipedia.org/wiki/Bloom_filter) support only point queries.

Probabilistic filters have many applications.  One well known application is to [LSM-trees](https://www.cs.umb.edu/~poneil/lsmtree.pdf).
An LSM-tree stores sorted runs of keys in files on disk.  A probabilistic filter allows us to ask "are we storing this key on disk?" without requiring us to make the trek to disk (because the compact filter can fit in memory).  This is a tremendous opimization because disk is slow!  EBS gp3 volumes have been measured as having an average latency of roughly 5ms (see [here](https://www.percona.com/blog/performance-of-various-ebs-storage-types-in-aws/)).  This is an eternity compared to processor speeds and the speed of memory. BloomRF allows us to ask "are we storing any key in this range on disk?" and therefore can dramatically improve the performance of LSM-Trees for serving range queries.  A good overview of LSM-trees can be found [here](https://cs-people.bu.edu/mathan/publications/icde23-tutorial.pdf).

This project is a work in progress.  In particular, the memory management optimization described in section 7,
and also the tuning advisor are not yet implemented.

## Compiling

```
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build
```

That builds unit tests, experiments, and the standalone bloomRF library.

To run the unit tests, run `./build/src/test/test_bloomrf`.  To run the experiments, run `./build/src/test/experiments`.

## Usage Example

`BloomRF` supports floats and integers.  The interface that `BloomRF` supports is `BloomRF<T>::add(T key)`,
`BloomRF<T>::find(T key)`, and `BloomRF<T>::findRange(T low, T high)`.

```
#include "bloomrf.h"

int main() {
    // TODO: simplify parameter tuning.  For now, these are reasonable defaults.
    // The first parameter is the memory allotment in bytes. For 2000000 keys, this is
    // 16-bits per key.
    BloomRF<uint64_t> ubf{BloomFilterRFParameters{4000000, 0, {7, 7, 7, 4, 4, 2, 2, 2}}};

    // Insert random keys.
    for (int i = 0; i < 2000000; ++i) {
        // rand is a bad random number generator.
        ubf.add(rand());
    }

    for (int i = 0; i < 100000000; ++i) {
        // Do random range queries of size 1e8.

        uint64_t low = rand();
        bool stored = ubf.findRange(low, low + 1e8);
        if (stored) {
            std::cout << "There may be a key in range [" << low << "," << high << "]" << std::endl;
        } else {
            std::cout << "We have definitely not stored a key in range [" << low << "," << high << "]" << std::endl;
        }
    }

    return 0;
}

```
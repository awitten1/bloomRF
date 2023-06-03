# bloomRF

This project is a work in progress.

This is an implementation of the EDBT best paper [bloomRF](https://openproceedings.org/2023/conf/edbt/paper-190.pdf).
BloomRF is a probabilistic filter that supports range queries, in addition to point queries.  Traditional
[bloom filters](https://en.wikipedia.org/wiki/Bloom_filter) support only point queries.

## Compiling

```
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build
```

That builds unit tests, experiments, and the standalone bloomRF library.

## Usage

BloomRF supports floats and integers.

```
int main() {
    BloomRF<uint64_t> ubf;

    // Insert keys.
    for (int i = 0; i < 100000) {
        ubf.add(rand() % 100000);
    }

    for ()
}

```
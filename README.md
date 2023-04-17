# bloomRF

This project is a work in progress.

This is an implementation of [bloomRF](https://openproceedings.org/2023/conf/edbt/paper-190.pdf), a probabilistic
filter that supports range queries, in addition to point queries.

## Compiling

```
mkdir build
cd build
cmake ..
make
```

That builds unit tests, experiments, and the standalone bloomRF library.
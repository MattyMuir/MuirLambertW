# MuirLambertW

Fast and accurate implementations of the Lambert-W Function for a variety of different contexts. These implementations are the first to provide robust accuracy guarantees for all possible inputs, without sacrificing performance.

Implementations are provided for two different precisions (`float` and `double`) for both function branches ($W_{-1}$ and $W_{0}$) with versions for both serial and SIMD evaluation using the AVX2 instruction set, for a total of 8 implementations:

```C++
float MuirW0(float x);
double MuirW0(double x);
__m256 MuirW0(__m256 x);
__m256d MuirW0(__m256d x);

float MuirWm1(float x);
double MuirWm1(double x);
__m256 MuirWm1(__m256 x);
__m256d MuirWm1(__m256d x);
```

### Accuracy
All functions guarantee results within a fixed number of Units in the Last Place (ULPs) from the exact answers. Currently the following bounds have been verified with Clang 19.1.0 on Windows

| Function                 | Maximum ULP Error |
| ------------------------ | ----------------- |
| $W_0$ `float`            | 1                 |
| $W_0$ `double`           | 4                 |
| $W_0$ `float` (SIMD)     | 4                 |
| $W_0$ `double` (SIMD)    | 4                 |
| $W_{-1}$ `float`         | 1                 |
| $W_{-1}$ `double`        | 4                 |
| $W_{-1}$ `float` (SIMD)  | 4                 |
| $W_{-1}$ `double` (SIMD) | 4                 |

## Dependencies
The main implementation has no dependencies, but the benchmarks and tests require
- Boost
- MPFR (And PkgConfig to manage)
- ReferenceLambertW
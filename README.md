# MuirLambertW

A fast and accurate AVX2 implementation of the real Lambert-W Function. Serial versions are also provided,
but I would generally recommend using Boost instead if accuracy is not a major concern and you only
need serial evaluation.

## Dependencies
The main implementation has no dependencies, but the benchmarks and tests require
- Boost
- MPFR (And PkgConfig to manage)
- ReferenceLambertW
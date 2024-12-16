#pragma once
#include <immintrin.h>

float MuirW0(float x);
double MuirW0(double x);
__m256 MuirW0(__m256 x);
__m256 MuirW0v2(__m256 x);
__m256d MuirW0(__m256d x);

float MuirWm1(float x);
double MuirWm1(double x);
__m256 MuirWm1(__m256 x);
__m256d MuirWm1(__m256d x);
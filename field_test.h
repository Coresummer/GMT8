#pragma once

#ifndef FIELD_TEST_H
#define FIELD_TEST_H

// #include "count.h"
// #include "fp6.h"
#include "./time.h"
#include "miller.h"

int test_fp(int fp_n);
int test_fp2(int fp2_n);
int test_fp6(int fp6_n);
int test_field(int fp, int fp2, int fp6, int sqr);
int test_fp_montgomery(int fp_n);

void check_fp_with_montgomery();
void check_fp2_with_montgomery();
void check_fp4_with_montgomery();
void check_fp8_with_montgomery();

void check_fp2_count();
void check_fp4_count();
void check_fp8_count();

void BENCH_fp2_fp4_fp8_mul_lazy_montgomery(int LOOP);
void BENCH_miller_lazy_montgomery(int LOOP);
void BENCH_miller_coordinates(int LOOP);
void BENCH_finalexp_lazy_montgomery(int LOOP);
void BENCH_Pairing_jac_lazy_montgomery(int LOOP);
void BENCH_Pairing_proj_lazy_montgomery(int LOOP);
void BENCH_Pairing_lazy_montgomery2(int LOOP);

void check_finalexp_pow_cost_count_2NAF();
void check_finalexp_pow_cost_count_2NAF_montgomery();
#endif

#pragma once
#ifndef fp8_H
#define fp8_H

#include "fp4.h"

void fp8_init(fp8_t *A);
void fp8_printf(std::string str ,fp8_t *A);
void fp8_println(std::string str ,fp8_t *A);
void fpd8_println(std::string str ,fpd8_t *A);
void fp8_printf_montgomery(std::string str ,fp8_t *A);
void fp8_println_montgomery(std::string str, fp8_t *A);
void fp8_set(fp8_t *ANS,fp8_t *A);
void fp8_set_ui(fp8_t *ANS,unsigned long int UI);
void fp8_set_ui_ui(fp8_t *ANS,unsigned long int UI);
void fp8_set_mpn(fp8_t *ANS,mp_limb_t *A);
void fp8_set_neg(fp8_t *ANS,fp8_t *A);
void fp8_set_conj(fp8_t *ANS,fp8_t *A);
void fp8_set_conj_montgomery(fp8_t *ANS,fp8_t *A);
void fp8_set_conj_montgomery_fpd(fpd8_t *ANS,fp8_t *A);
void fp8_to_montgomery(fp8_t *ANS,fp8_t *A);
void fp8_mod_montgomery(fp8_t *ANS,fp8_t *A);
void fp8_mod_montgomery_double(fp8_t *ANS,fpd8_t *A);
void fp8_r1shift(fp8_t *ANS, fp8_t *A);
void fp8_lshift(fp8_t *ANS, fp8_t *A, unsigned long int UI);
void fp8_l1shift(fp8_t *ANS, fp8_t *A);
void fp8_l1shift_nonmod_single(fp8_t *ANS, fp8_t *A);
void fp8_l1shift_nonmod_double(fpd8_t *ANS, fpd8_t *A);
void fp8_hlv(fp8_t *ANS,fp8_t *A);
void fp8_set_random(fp8_t *ANS,gmp_randstate_t state);
void fp8_mul(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_mul_lazy(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_mul_lazy_montgomery(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_mul_nonmod_montgomery(fpd8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_mul_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI);
void fp8_mul_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B);
void fp8_mul_mpn_montgomery(fp8_t *ANS,fp8_t *A,mp_limb_t *B);
void fp8_sqr(fp8_t *ANS,fp8_t *A);
void fp8_sqr_final(fp8_t *ANS,fp8_t *A);
void fp8_sqr_lazy(fp8_t *ANS,fp8_t *A);
void fp8_sqr_lazy_montgomery(fp8_t *ANS,fp8_t *A);
void fp8_sqr_nonmod_montgomery(fpd8_t *ANS, fp8_t *A);
void fp8_add(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_add_nonmod_single(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_add_nonmod_double(fpd8_t *ANS,fpd8_t *A,fpd8_t *B);
void fp8_add_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI);
void fp8_add_ui_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI);
void fp8_add_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B);
void fp8_sub(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_sub_nonmod_single(fp8_t *ANS,fp8_t *A,fp8_t *B);
void fp8_sub_nonmod_double(fpd8_t *ANS,fpd8_t *A,fpd8_t *B);
void fp8_sub_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI);
void fp8_sub_ui_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI);
void fp8_sub_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B);
void fp8_inv(fp8_t *ANS,fp8_t *A);
void fp8_inv_lazy(fp8_t *ANS, fp8_t *A);
void fp8_inv_lazy_montgomery(fp8_t *ANS,fp8_t *A);
int fp8_legendre(fp8_t *A);
void fp8_sqrt(fp8_t *ANS,fp8_t *A);
void fp8_pow(fp8_t *ANS,fp8_t *A,mpz_t scalar);
void fp8_pow_montgomery(fp8_t *ANS, fp8_t *A, mpz_t scalar);
int fp8_cmp(fp8_t *A,fp8_t *B);
int fp8_cmp_ui(fp8_t *A,unsigned long int UI);
int fp8_cmp_mpn(fp8_t *A,mp_limb_t *B);
int fp8_cmp_zero(fp8_t *A);
int fp8_cmp_one(fp8_t *A);
int fp8_montgomery_trick_montgomery(fp8_t *A_inv,fp8_t *A,int n);

void fp8_frobenius_map_p1(fp8_t *ANS,fp8_t *A);
void fp8_frobenius_map_p2(fp8_t *ANS,fp8_t *A);
void fp8_frobenius_map_p3(fp8_t *ANS,fp8_t *A);

void fp8_mul_base(fp8_t *ANS,fp8_t *A);
void fp8_mul_base_nonmod_single(fp8_t *ANS,fp8_t *A);
void fp8_mul_base_nonmod_double(fpd8_t *ANS,fpd8_t *A);
#endif

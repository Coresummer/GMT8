#pragma once
#ifndef fp4_H
#define fp4_H

#include "fp2.h"

void fp4_init(fp4_t *A);
void fp4_printf(std::string str ,fp4_t *A);
void fp4_println(std::string str ,fp4_t *A);
void fpd4_printf(std::string str ,fpd4_t *A);
void fpd4_println(std::string str ,fpd4_t *A);
void fp4_printf_montgomery(std::string str ,fp4_t *A);
void fp4_println_montgomery(std::string str, fp4_t *A);
void fp4_set(fp4_t *ANS,fp4_t *A);
void fp4_set_ui(fp4_t *ANS,unsigned long int UI);
void fp4_set_ui_ui(fp4_t *ANS,unsigned long int UI);
void fp4_set_mpn(fp4_t *ANS,mp_limb_t *A);
void fp4_set_neg(fp4_t *ANS,fp4_t *A);
void fp4_set_neg_montgomery(fp4_t *ANS,fp4_t *A);
void fp4_set_conj(fp4_t *ANS,fp4_t *A);
void fp4_set_conj_montgomery(fp4_t *ANS,fp4_t *A);
void fp4_set_conj_montgomery_fpd(fpd4_t *ANS,fp4_t *A);
void fpd4_set(fpd4_t *ANS, fpd4_t *A) ;
void fp4_to_montgomery(fp4_t *ANS,fp4_t *A);
void fp4_mod_montgomery(fp4_t *ANS,fp4_t *A);
void fp4_mod_montgomery_double(fp4_t *ANS,fpd4_t *A);
void fp4_r1shift(fp4_t *ANS, fp4_t *A);
void fp4_lshift(fp4_t *ANS, fp4_t *A, unsigned long int UI);
void fp4_l1shift(fp4_t *ANS, fp4_t *A);
void fp4_l1shift_nonmod_single(fp4_t *ANS, fp4_t *A);
void fp4_l1shift_nonmod_double(fpd4_t *ANS, fpd4_t *A);
void fp4_hlv(fp4_t *ANS,fp4_t *A);
void fp4_set_random(fp4_t *ANS,gmp_randstate_t state);
void fp4_mul(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_mul_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_mul_lazy_montgomery(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_mul_nonmod_montgomery(fpd4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_mul_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);
void fp4_mul_ui_nonmod_single(fp4_t *ANS, fp4_t *A, unsigned long int UI);
void fp4_mul_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);
void fp4_mul_mpn_montgomery(fp4_t *ANS,fp4_t *A,mp_limb_t *B);
void fp4_sqr(fp4_t *ANS,fp4_t *A);
void fp4_sqr_GS(fp4_t *ANS,fp4_t *A);
void fp4_sqr_lazy(fp4_t *ANS,fp4_t *A);
void fp4_sqr_lazy_montgomery(fp4_t *ANS,fp4_t *A);
void fp4_sqr_nonmod_montgomery(fpd4_t *ANS, fp4_t *A);
void fp4_add(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_add_nonmod_single(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_add_nonmod_double(fpd4_t *ANS,fpd4_t *A,fpd4_t *B);
void fp4_add_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);
void fp4_add_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);
void fp4_add_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);
void fp4_sub(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_sub_nonmod_single(fp4_t *ANS,fp4_t *A,fp4_t *B);
void fp4_sub_nonmod_double(fpd4_t *ANS,fpd4_t *A,fpd4_t *B);
void fp4_sub_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);
void fp4_sub_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);
void fp4_sub_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);
void fp4_inv(fp4_t *ANS,fp4_t *A);
void fp4_inv_lazy(fp4_t *ANS, fp4_t *A);
void fp4_inv_lazy_montgomery(fp4_t *ANS,fp4_t *A);
int fp4_legendre(fp4_t *A);
void fp4_sqrt(fp4_t *ANS,fp4_t *A);
void fp4_pow(fp4_t *ANS,fp4_t *A,mpz_t scalar);
void fp4_pow_montgomery(fp4_t *ANS, fp4_t *A, mpz_t scalar);
int fp4_cmp(fp4_t *A,fp4_t *B);
int fp4_cmp_ui(fp4_t *A,unsigned long int UI);
int fp4_cmp_mpn(fp4_t *A,mp_limb_t *B);
int fp4_cmp_zero(fp4_t *A);
int fp4_cmp_one(fp4_t *A);
int fp4_montgomery_trick_montgomery(fp4_t *A_inv,fp4_t *A,int n);

void fp4_frobenius_map_p1(fp4_t *ANS,fp4_t *A);
void fp4_frobenius_map_p2(fp4_t *ANS,fp4_t *A);
void fp4_frobenius_map_p3(fp4_t *ANS,fp4_t *A);

void fp4_mul_base(fp4_t *ANS,fp4_t *A);
void fp4_mul_base_inv(fp4_t *ANS,fp4_t *A);
void fp4_mul_base_nonmod_single(fp4_t *ANS,fp4_t *A);
void fp4_mul_base_nonmod_double(fpd4_t *ANS,fpd4_t *A);
#endif

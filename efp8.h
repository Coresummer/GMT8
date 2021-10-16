#pragma once
#ifndef efp8_H
#define efp8_H

#include "efp4.h"

void efp8_init(efp8_t *P);
void efp8_projective_init(efp8_projective_t *P);
void efp8_jacobian_init(efp8_jacobian_t *P);
void efp8_printf(std::string str ,efp8_t *P);
void efp8_println(std::string str ,efp8_t *P);
void efp8_projective_printf(std::string str ,efp8_projective_t *P);
void efp8_jacobian_printf(std::string str ,efp8_jacobian_t *P);
void efp8_printf_montgomery(std::string str ,efp8_t *P);
void efp8_jacobian_printf_montgomery(std::string str ,efp8_jacobian_t *P);
void efp8_projective_printf_montgomery(std::string str ,efp8_projective_t *P);
void efp8_projective_printf_affine(std::string str ,efp8_projective_t *P);
void efp8_projective_printf_affine_montgomery(std::string str ,efp8_projective_t *P);
void efp8_set(efp8_t *ANS,efp8_t *A);
void efp8_projective_set(efp8_projective_t *ANS,efp8_projective_t *A);
void efp8_jacobian_set(efp8_jacobian_t *ANS,efp8_jacobian_t *A);
void efp8_affine_to_projective(efp8_projective_t *ANS,efp8_t *A);
void efp8_affine_to_jacobian(efp8_jacobian_t *ANS,efp8_t *A);
void efp8_affine_to_projective_montgomery(efp8_projective_t *ANS,efp8_t *A);
void efp8_affine_to_jacobian_montgomery(efp8_jacobian_t *ANS,efp8_t *A);
void efp8_jacobian_to_affine(efp8_t *ANS,efp8_jacobian_t *A);
void efp8_projective_to_affine(efp8_t *ANS,efp8_projective_t *A);
void efp8_jacobian_to_affine_montgomery(efp8_t *ANS,efp8_jacobian_t *A);
void efp8_projective_to_affine_montgomery(efp8_t *ANS,efp8_projective_t *A);
void efp8_mix(efp8_jacobian_t *ANS,efp8_jacobian_t *A,fp2_t *Zi);
void efp8_mix_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *A,fp2_t *Zi);
void efp8_set_ui(efp8_t *ANS,unsigned long int UI);
void efp8_to_montgomery(efp8_t *ANS,efp8_t *A);
void efp8_projective_to_montgomery(efp8_projective_t *ANS,efp8_projective_t *A);
void efp8_mod_montgomery(efp8_t *ANS,efp8_t *A);
void efp8_projective_mod_montgomery(efp8_projective_t *ANS,efp8_projective_t *A);
void efp8_set_mpn(efp8_t *ANS,mp_limb_t *A);
void efp8_set_neg(efp8_t *ANS,efp8_t *A);
void efp8_jacobian_set_neg(efp8_jacobian_t *ANS,efp8_jacobian_t *A);
int efp8_cmp(efp8_t *A,efp8_t *B);
void efp8_rational_point(efp8_t *P);
void generate_g1(efp8_t *P);
void generate_g2(efp8_t *Q);
void efp8_ecd(efp8_t *ANS,efp8_t *P);
void efp8_ecd_jacobian_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P);
void efp8_eca(efp8_t *ANS,efp8_t *P1,efp8_t *P2);
void efp8_eca_jacobian_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P1,efp8_jacobian_t *P2);
void efp8_eca_mixture_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P1,efp8_jacobian_t *P2);
void efp8_scm(efp8_t *ANS,efp8_t *P,mpz_t scalar);
void efp8_scm_dash(efp8_t *ANS,efp8_t *P,mpz_t scalar);

void efp8_frobenius_map_p1(efp8_t *ANS,efp8_t *A);
void efp8_checkOnCurve(efp8_t *A);
void efp8_checkOnTwsitCurve(efp8_t *A);
#endif

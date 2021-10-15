#pragma once
#ifndef efp4_H
#define efp4_H

#include "efp2.h"

void efp4_init(efp4_t *P);
void efp4_projective_init(efp4_projective_t *P);
void efp4_jacobian_init(efp4_jacobian_t *P);
void efp4_printf(std::string str ,efp4_t *P);
void efp4_println(std::string str ,efp4_t *P);
void efp4_projective_printf(std::string str ,efp4_projective_t *P);
void efp4_jacobian_printf(std::string str ,efp4_jacobian_t *P);
void efp4_printf_montgomery(std::string str ,efp4_t *P);
void efp4_jacobian_printf_montgomery(std::string str ,efp4_jacobian_t *P);
void efp4_projective_printf_montgomery(std::string str ,efp4_projective_t *P);
void efp4_projective_printf_affine(std::string str ,efp4_projective_t *P);
void efp4_projective_printf_affine_montgomery(std::string str ,efp4_projective_t *P);
void efp4_set(efp4_t *ANS,efp4_t *A);
void efp4_projective_set(efp4_projective_t *ANS,efp4_projective_t *A);
void efp4_jacobian_set(efp4_jacobian_t *ANS,efp4_jacobian_t *A);
void efp4_affine_to_projective(efp4_projective_t *ANS,efp4_t *A);
void efp4_affine_to_jacobian(efp4_jacobian_t *ANS,efp4_t *A);
void efp4_affine_to_projective_montgomery(efp4_projective_t *ANS,efp4_t *A);
void efp4_affine_to_jacobian_montgomery(efp4_jacobian_t *ANS,efp4_t *A);
void efp4_jacobian_to_affine(efp4_t *ANS,efp4_jacobian_t *A);
void efp4_projective_to_affine(efp4_t *ANS,efp4_projective_t *A);
void efp4_jacobian_to_affine_montgomery(efp4_t *ANS,efp4_jacobian_t *A);
void efp4_projective_to_affine_montgomery(efp4_t *ANS,efp4_projective_t *A);
void efp4_mix(efp4_jacobian_t *ANS,efp4_jacobian_t *A,fp2_t *Zi);
void efp4_mix_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *A,fp2_t *Zi);
void efp4_set_ui(efp4_t *ANS,unsigned long int UI);
void efp4_to_montgomery(efp4_t *ANS,efp4_t *A);
void efp4_projective_to_montgomery(efp4_projective_t *ANS,efp4_projective_t *A);
void efp4_mod_montgomery(efp4_t *ANS,efp4_t *A);
void efp4_projective_mod_montgomery(efp4_projective_t *ANS,efp4_projective_t *A);
void efp4_set_mpn(efp4_t *ANS,mp_limb_t *A);
void efp4_set_neg(efp4_t *ANS,efp4_t *A);
void efp4_jacobian_set_neg(efp4_jacobian_t *ANS,efp4_jacobian_t *A);
int efp4_cmp(efp4_t *A,efp4_t *B);
void efp4_rational_point(efp4_t *P);
void efp4_ecd(efp4_t *ANS,efp4_t *P);
void efp4_ecd_jacobian_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P);
void efp4_eca(efp4_t *ANS,efp4_t *P1,efp4_t *P2);
void efp4_eca_jacobian_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P1,efp4_jacobian_t *P2);
void efp4_eca_mixture_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P1,efp4_jacobian_t *P2);
void efp4_scm(efp4_t *ANS,efp4_t *P,mpz_t scalar);


#endif

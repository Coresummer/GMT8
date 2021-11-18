
#pragma once
#ifndef MILLER_H
#define MILLER_H

#include "create.h"
#include "efp8.h"

void efp8_to_Jacefp2(efp2_jacobian_t *ANS,efp8_t *A);
void efp8_to_Jacefp2_montgomery(efp2_jacobian_t *ANS,efp8_t *A);
void efp8_to_efp2(efp2_t *ANS,efp8_t *A);
void efp8_to_efp2_montgomery(efp2_t *ANS,efp8_t *A);
void ff_lttp(fp8_t *f, efp2_jacobian_t *S, efp_t *P);
void ff_lttp_Costello(fp8_t *f, efp2_jacobian_t *U, efp_t *S);
void ff_ltqp(fp8_t *f, efp2_jacobian_t *S, efp2_t *Q,efp_t *P);
void ff_ltqp_Costello_mixed(fp8_t *f, efp2_jacobian_t *U, efp2_t *R,efp_t *S);
void ff_lttp_lazy_montgomery(fp8_t *f, efp2_jacobian_t *S, efp_t *P);
void ff_ltqp_lazy_montgomery(fp8_t *f, efp2_jacobian_t *S, efp2_t *Q,efp_t *P);

void miller_proj_precomp_Costello(efp_t *S,efp2_t *R);

void miller_opt_ate_jac(fp8_t *f,efp8_t *P,efp8_t *Q);
void miller_proj_ate_jac(fp8_t *f,efp8_t *P,efp8_t *Q);

void miller_opt_ate_jac_2NAF(fp8_t *f,efp8_t *P,efp8_t *Q);
void miller_opt_ate_proj_2NAF(fp8_t *f,efp8_t *P,efp8_t *Q);


void miller_opt_ate_jac_2NAF_lazy_montgomery(fp8_t *f,efp8_t *P,efp8_t *Q);
void miller_opt_ate_proj_2NAF_lazy_montgomery(fp8_t *f,efp8_t *P,efp8_t *Q);
#endif
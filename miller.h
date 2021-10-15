#pragma once
#ifndef MILLER_H
#define MILLER_H

#include "create.h"
#include "efp8.h"

void efp8_to_Jacefp(efp_jacobian_t *ANS,efp8_t *A);
void efp8_to_Jacefp_montgomery(efp_jacobian_t *ANS,efp8_t *A);
void efp8_to_efp(efp_t *ANS,efp8_t *A);
void efp8_to_efp_montgomery(efp_t *ANS,efp8_t *A);
void ff_lttp(fp8_t *f, efp_jacobian_t *S, efp_t *P);
void ff_ltqp(fp8_t *f, efp_jacobian_t *S, efp_t *Q,efp_t *P);
void ff_lttp_lazy_montgomery(fp8_t *f, efp_jacobian_t *S, efp_t *P);
void ff_ltqp_lazy_montgomery(fp8_t *f, efp_jacobian_t *S, efp_t *Q,efp_t *P);

void miller_opt_ate_proj(fp8_t *f,efp8_t *P,efp8_t *Q);
void miller_opt_ate_proj_2NAF(fp8_t *f,efp8_t *P,efp8_t *Q);
void miller_opt_ate_proj_2NAF_lazy_montgomery(fp8_t *f,efp8_t *P,efp8_t *Q);
#endif

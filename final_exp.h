#pragma once
#ifndef FINAL_EXP_H
#define FINAL_EXP_H

#include "miller.h"

// void final_exp_basic(fp8_t *ANS,fp8_t *A);
void final_exp_direct(fp8_t *ANS,fp8_t *A);
void final_exp(fp8_t *ANS,fp8_t *A);
void final_exp_lazy_montgomery(fp8_t *ANS,fp8_t *A);
void final_exp_lazy_montgomery2(fp8_t *ANS,fp8_t *A);
#endif

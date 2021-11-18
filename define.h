#pragma once
#ifdef TTT_INSTANCE_HERE
    #define TTT_EXTERN
#else
    #define TTT_EXTERN extern
#endif

#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <vector>
#include <string>
#include <iostream>
#include "./mcl.h"


#define ARCBIT 64  //64bit processor
//#define ARCBIT 32 //32bit processor

// #define DEBUG_COST_A
//#define DEBUG_ASSERT
#define CHECK_PAIRING_TIME_LOOP 100000

/**************Option**************/
#define X64
/**********************************/

#ifdef  X64
#define FPLIMB_BITS FPLIMB*ARCBIT
#define FPLIMB 9 //
#define FPLIMB2 FPLIMB*2   //?? + 1
#endif

#ifdef X32
#define FPLIMB_BITS FPLIMB* ARCBIT
#define FPLIMB 22
#define FPLIMB2 FPLIMB*2  //??
#endif

#define scalar_t mpz_t

#define k6_X_length 128//////37??

/*Field*/
typedef struct{
  mp_limb_t x0[FPLIMB];
}fp_t;
typedef struct{
  fp_t x0,x1;
}fp2_t;
typedef struct{
  fp2_t x0,x1;
}fp4_t;
typedef struct{
  fp4_t x0,x1;
}fp8_t;


/*dField*/
typedef struct{
  mp_limb_t x0[FPLIMB2];
}fpd_t;
typedef struct{
  fpd_t x0,x1;
}fpd2_t;
typedef struct{
  fpd2_t x0,x1;
}fpd4_t;
typedef struct{
  fpd4_t x0,x1;
}fpd8_t;
//tmp finite field
TTT_EXTERN mp_limb_t buf[FPLIMB];

/*Elliptic Curve*/
typedef struct{
  fp_t x,y;
  int infinity;
}efp_t;
typedef struct{
  fp2_t x,y;
  int infinity;
}efp2_t;
typedef struct{
  fp4_t x,y;
  int infinity;
}efp4_t;
typedef struct{
  fp8_t x,y;
  int infinity;
}efp8_t;

/*Projectiv Elliptic Curve*/
typedef struct{
  fp_t x,y,z;
  int infinity;
}efp_projective_t;
typedef struct{
  fp2_t x,y,z;
  int infinity;
}efp2_projective_t;
typedef struct{
  fp4_t x,y,z;
  int infinity;
}efp4_projective_t;
typedef struct{
  fp8_t x,y,z;
  int infinity;
}efp8_projective_t;

/*Jacobian Elliptic Curve*/
typedef struct{
  fp_t x,y,z;
  int infinity;
}efp_jacobian_t;
typedef struct{
  fp2_t x,y,z;
  int infinity;
}efp2_jacobian_t;
typedef struct{
  fp4_t x,y,z;
  int infinity;
}efp4_jacobian_t;
typedef struct{
  fp8_t x,y,z;
  int infinity;
}efp8_jacobian_t;


TTT_EXTERN gmp_randstate_t state;//for random
TTT_EXTERN int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_mul_ui,cost_mul_base,cost_mul_base_inv,cost_sqr,cost_inv,cost_mod,cost_set_neg;
TTT_EXTERN int cost_add_nonmod, cost_add_nonmod_double, cost_sub_nonmod, cost_sub_nonmod_double, cost_r1shift, cost_mod_nomal;
TTT_EXTERN mpz_t X_z,prime_z,order_z,trace_z;
TTT_EXTERN mp_limb_t X,prime[FPLIMB],few_prime[FPLIMB];
TTT_EXTERN mp_limb_t prime2[FPLIMB2];
TTT_EXTERN fp2_t base_c,base_cMR,base_c_inv, base_c_invMR,fp2_neg_1;
TTT_EXTERN fp_t oneMR, p_1,three_1,three_1MR;

TTT_EXTERN fp_t curve_a,curve_a_dash,curve_b;
TTT_EXTERN mpz_t sqrt_power_z;

TTT_EXTERN mpz_t efp_total,efp2_total,efp4_total,efp8_total,fp8_total_r;//#efp,#efp5,#efp10,#efp7,#efp14
TTT_EXTERN mpz_t miller_loop_s;
TTT_EXTERN std::vector<int> miller_loop_v, finalexp_pow_x,finalexp_pow_x_2, finalexp_pow_4hy, finalexp_pow_hy;
TTT_EXTERN mpz_t X_1_div2,X_1,X_2,X_2_1,four;//(kai +1)/2,(kai -1),(kai^2)をあらかじめ求めておく
TTT_EXTERN mpz_t hardpart,hy_neg,fourhy_neg,three;
//emb6
TTT_EXTERN fp2_t frobenius_1_4, frobenius_2_4, frobenius_3_4;     //c^((p-1)/10)の計算結果
TTT_EXTERN fp2_t frobenius_1_4MR,frobenius_2_4MR,frobenius_3_4MR; //c^((p-1)/10)の計算結果

//miller precomp
TTT_EXTERN fp2_t miller_precomp_X2Z2, miller_precomp_xSZ2_2,miller_precomp_ySZ2_2,miller_precomp_X2_Y2,miller_precomp_xS_Y2,miller_precomp_yS_Y2;

//montgomery
TTT_EXTERN mp_limb_t R[FPLIMB],Ri[FPLIMB],R1[FPLIMB],RR[FPLIMB],Ni[FPLIMB];
TTT_EXTERN int m;
TTT_EXTERN mp_limb_t u[FPLIMB+1];
TTT_EXTERN mp_limb_t N[FPLIMB2],R2[FPLIMB],R3[FPLIMB],RmodP[FPLIMB];
TTT_EXTERN mp_limb_t Ni_neg;  //Ni_neg=-N^(-1)

TTT_EXTERN struct timeval tv_start,tv_end;
TTT_EXTERN float MILLER_ATE_4SPARSE_TIME;
TTT_EXTERN float MILLER_ATE_5SPARSE_TIME;
TTT_EXTERN float MILLER_ATE_6SPARSE_TIME;
TTT_EXTERN float MILLER_ATE_7SPARSE_TIME;
TTT_EXTERN float FINAL_EXP_TIME;

typedef struct {
  int add;
  int add_ui;
  int add_nonmod;
  int add_nonmod_double;
  int sub;
  int sub_ui;
  int sub_nonmod;
  int sub_nonmod_double;
  int mul;
  int set_neg;
  int r1shift;
  int sqr;
  int inv;
  int mul_base;
  int mul_base_inv;
  int mod;
  int mod_nomal;
} cost;

#endif

#include "fp.h"
#include "fp2.h"
#include "fp4.h"
#include "fp8.h"
#include <cstdio>
#define CYBOZU_BENCH_USE_GETTIMEOFDAY
#include <cybozu/benchmark.hpp>

#include "field_test.h"

#include "define.h"
#include <gmp.h>
#include "final_exp.h"
/*----------------------------------------------------------------------------*/
//test
int test_fp(int fp_n) {
  int i, j, n = 0;
  float add_time = 0, add_lazy_time = 0;
  float sub_time = 0, sub_lazy_time = 0;
  float mul_time = 0, mul_lazy_time = 0;
  float inv_time = 0, inv_lazy_time = 0;
  float sqr_time = 0, sqr_lazy_time = 0;
  struct timeval tv_A, tv_B;
  mp_limb_t test1_mpn[FPLIMB2], test2_mpn[FPLIMB2];
  fp_t A_fp, B_fp, test1_fp, test2_fp;

  fp_init(&A_fp);
  fp_init(&B_fp);
  fp_init(&test1_fp);
  fp_init(&test2_fp);

  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_add test\n");
  add_time = 0, add_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_add(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    add_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    //for(j=0;j<n;j++)fp_add_lazy(test2_mpn,FPLIMB,A_fp.x0,FPLIMB,B_fp.x0,FPLIMB);
    gettimeofday(&tv_B, NULL);
    add_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp_mod(&test2_fp, test2_mpn, FPLIMB);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp add.      : %.7f[ms]\n", add_time / fp_n);
  printf("fp add lazy. : %.7f[ms]\n", add_lazy_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_sub test\n");
  sub_time = 0, sub_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sub(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    sub_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    //for(j=0;j<n;j++)fp_sub_lazy(test2_mpn,FPLIMB,A_fp.x0,FPLIMB,B_fp.x0,FPLIMB);
    gettimeofday(&tv_B, NULL);
    sub_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp_mod(&test2_fp, test2_mpn, FPLIMB);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp sub.      : %.7f[ms]\n", sub_time / fp_n);
  printf("fp sub lazy. : %.7f[ms]\n", sub_lazy_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_mul test\n");
  mul_time = 0, mul_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_mul(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    mul_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    //for(j=0;j<n;j++)fp_mul_lazy(test2_mpn,A_fp.x0,B_fp.x0);
    gettimeofday(&tv_B, NULL);
    mul_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp_mod(&test2_fp, test2_mpn, FPLIMB2);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }

  printf("fp mul.      : %.6f[ms]\n", mul_time / fp_n);
  printf("fp mul lazy. : %.6f[ms]\n", mul_lazy_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_sqr test\n");
  sqr_time = 0, sqr_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sqr(&test1_fp, &A_fp);
    gettimeofday(&tv_B, NULL);
    sqr_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    //for(j=0;j<n;j++)fp_sqr_lazy(test2_mpn,A_fp.x0);
    gettimeofday(&tv_B, NULL);
    sqr_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp_mod(&test2_fp, test2_mpn, FPLIMB2);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp sqr.      : %.6f[ms]\n", sqr_time / fp_n);
  printf("fp sqr lazy. : %.6f[ms]\n", sqr_lazy_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_inv test\n");
  inv_time = 0, inv_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_inv(&test1_fp, &A_fp);
    gettimeofday(&tv_B, NULL);
    inv_time += timedifference_msec(tv_A, tv_B) / n;
  }
  printf("fp inv.      : %.6f[ms]\n", inv_time / fp_n);

  return 0;
}

int test_fp2(int fp2_n) {
  int i, j, n = 0;
  float add_time = 0, add_lazy_time = 0;
  float sub_time = 0, sub_lazy_time = 0;
  float mul_time = 0, mul_lazy_time = 0;
  float inv_time = 0, inv_lazy_time = 0;
  float sqr_time = 0, sqr_lazy_time = 0;
  struct timeval tv_A, tv_B;
  fp2_t A_fp2, B_fp2, test1_fp2, test2_fp2;
  fpd2_t test2_fpd2;

  fp2_init(&A_fp2);
  fp2_init(&B_fp2);
  fp2_init(&test1_fp2);
  fp2_init(&test2_fp2);

  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  gmp_randseed_ui(state, 1);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp2_mul test\n");
  mul_time = 0, mul_lazy_time = 0;
  n = 1;
  for (i = 0; i < fp2_n; i++) {
    fp2_set_random(&A_fp2, state);
    fp2_set_random(&B_fp2, state);
    fp2_to_montgomery(&A_fp2, &A_fp2);
    fp2_to_montgomery(&B_fp2, &B_fp2);

    gettimeofday(&tv_A, NULL);
    fp2_mul_lazy_montgomery(&test1_fp2, &A_fp2, &B_fp2);
    gettimeofday(&tv_B, NULL);
    mul_time += timedifference_msec(tv_A, tv_B) / n;
    fp2_to_montgomery(&B_fp2, &B_fp2);

    gettimeofday(&tv_A, NULL);
    fp2_mul_nonmod_montgomery(&test2_fpd2, &A_fp2, &B_fp2);
    gettimeofday(&tv_B, NULL);
    mul_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp2_mod_montgomery_double(&test2_fp2, &test2_fpd2);
    fp2_mod_montgomery(&test2_fp2, &test2_fp2);

    if (fp2_cmp(&test1_fp2, &test2_fp2) != 0) {
      printf("failed!\n\n");
      fp2_printf("", &test1_fp2);
      fp2_printf("\n", &test2_fp2);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp2 sqr.      : %.6f[ms]\n", sqr_time / fp2_n);
  printf("fp2 sqr lazy. : %.6f[ms]\n", sqr_lazy_time / fp2_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp2_sqr test\n");
  sqr_time = 0, sqr_lazy_time = 0;
  n = 1;
  for (i = 0; i < fp2_n; i++) {
    fp2_set_random(&A_fp2, state);
    fp2_set_random(&B_fp2, state);
    fp2_to_montgomery(&A_fp2, &A_fp2);
    fp2_to_montgomery(&B_fp2, &B_fp2);

    gettimeofday(&tv_A, NULL);
    fp2_sqr_lazy_montgomery(&test1_fp2, &A_fp2);
    gettimeofday(&tv_B, NULL);
    sqr_time += timedifference_msec(tv_A, tv_B) / n;
    //fp2_to_montgomery(&B_fp2,&B_fp2);

    gettimeofday(&tv_A, NULL);
    fp2_sqr_nonmod_montgomery(&test2_fpd2, &A_fp2);
    gettimeofday(&tv_B, NULL);
    sqr_lazy_time += timedifference_msec(tv_A, tv_B) / n;
    fp2_mod_montgomery_double(&test2_fp2, &test2_fpd2);
    fp2_mod_montgomery(&test2_fp2, &test2_fp2);

    if (fp2_cmp(&test1_fp2, &test2_fp2) != 0) {
      printf("failed!\n\n");
      fp2_printf("", &test1_fp2);
      fp2_printf("\n", &test2_fp2);
      printf("\n\n");
      return 1;
    }
  }

  printf("fp2 sqr.      : %.6f[ms]\n", sqr_time / fp2_n);
  printf("fp2 sqr lazy. : %.6f[ms]\n", sqr_lazy_time / fp2_n);
  return 0;
}

int test_fp_montgomery(int fp_n) {
  int i, j, n = 0;
  float add_time = 0, add_nonmod_single = 0, add_nonmod_double = 0;
  float lshift_nonmod_time = 0, l1shift_nonmod_time = 0, lshift3_nonmod_time = 0;
  float sub_time = 0, sub_nonmod_single = 0, sub_nonmod_double = 0;
  float mul_time = 0, mul_nonmod_time = 0;
  float inv_time = 0, inv_lazy_time = 0;
  float sqr_time = 0, sqr_nonmod = 0;
  float mod_time = 0;
  struct timeval tv_A, tv_B;
  mp_limb_t test1_mpn[FPLIMB2], test2_mpn[FPLIMB2];
  fp_t A_fp, B_fp, test1_fp, test2_fp;
  fpd_t A_fpd, B_fpd, test1_fpd, test2_fpd;

  fp_init(&A_fp);
  fp_init(&B_fp);
  fp_init(&test1_fp);
  fp_init(&test2_fp);

  // fpd_init(&A_fpd);
  // fpd_init(&B_fpd);
  // fpd_init(&test1_fpd);
  // fpd_init(&test2_fpd);

  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_add test\n");
  add_time = 0, add_nonmod_single = 0, add_nonmod_double = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random_montgomery(&A_fp, state);
    fp_set_random_montgomery(&B_fp, state);
    fp_mul_nonmod(&A_fpd, &A_fp, &B_fp);
    fp_set_random_montgomery(&A_fp, state);
    fp_set_random_montgomery(&B_fp, state);
    fp_mul_nonmod(&B_fpd, &A_fp, &B_fp);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_add(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    add_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_add_nonmod_single(&test2_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    add_nonmod_single += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_add_nonmod_double(&A_fpd, &A_fpd, &B_fpd);
    gettimeofday(&tv_B, NULL);
    add_nonmod_double += timedifference_msec(tv_A, tv_B) / n;

    fp_mod_montgomery(&test1_fp, &test1_fp);
    fp_mod_montgomery(&test2_fp, &test2_fp);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp add.               : %.7f[ms]\n", add_time / fp_n);
  printf("fp add nonmod single. : %.7f[ms]\n", add_nonmod_single / fp_n);
  printf("fp add nonmod double. : %.7f[ms]\n", add_nonmod_double / fp_n);

  // printf("------------------------------------------------------------------------------------\n");
  //     printf("fp_lshift test\n");
  // lshift_nonmod_time=0,l1shift_nonmod_time=0,lshift3_nonmod_time=0;
  // n=1000;
  // for(i=0;i<fp_n;i++){
  //     fp_set_random_montgomery(&A_fp,state);
  //     fp_set_random_montgomery(&B_fp,state);

  //     gettimeofday(&tv_A,NULL);
  //     for(j=0;j<n;j++)fp_lshift_ui_nonmod_single(&test1_fp,&A_fp,1);
  //     //for(j=0;j<n;j++)fp_mul_ui_nonmod_single(&test1_fp,&A_fp,2);
  //     gettimeofday(&tv_B,NULL);
  //     lshift_nonmod_time+=timedifference_msec(tv_A,tv_B)/n;

  //     gettimeofday(&tv_A,NULL);
  //     for(j=0;j<n;j++)fp_lshift_ui_nonmod_single(&test2_fp,&A_fp,2);
  //     gettimeofday(&tv_B,NULL);
  //     l1shift_nonmod_time+=timedifference_msec(tv_A,tv_B)/n;

  //     // fp_mod_montgomery(&test1_fp,&test1_fp);
  //     // fp_mod_montgomery(&test2_fp,&test2_fp);

  //     // if(fp_cmp(&test1_fp,&test2_fp)!=0){
  //     //     printf("failed!\n\n");
  //     // 	fp_printf("",&test1_fp);
  // 	//     fp_printf("\n",&test2_fp);
  // 	//     printf("\n\n");
  // 	//     return 1;
  //     // }
  // }
  //     // printf("fp lshift1.        : %.7f[ms]\n",lshift_nonmod_time/fp_n);
  //     // printf("fp l1shift. : %.7f[ms]\n",l1shift_nonmod_time/fp_n);
  //     // printf("fp lshift3. : %.7f[ms]\n",lshift3_nonmod_time/fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_sub test\n");
  sub_time = 0, sub_nonmod_single = 0, sub_nonmod_double = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random_montgomery(&A_fp, state);
    fp_set_random_montgomery(&B_fp, state);
    fp_mul_nonmod(&A_fpd, &A_fp, &B_fp);
    fp_set_random_montgomery(&A_fp, state);
    fp_set_random_montgomery(&B_fp, state);
    fp_mul_nonmod(&B_fpd, &A_fp, &B_fp);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sub(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    sub_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sub_nonmod_single(&test2_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    sub_nonmod_single += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sub_nonmod_double(&A_fpd, &A_fpd, &B_fpd);
    gettimeofday(&tv_B, NULL);
    sub_nonmod_double += timedifference_msec(tv_A, tv_B) / n;

    fp_mod_montgomery(&test1_fp, &test1_fp);
    fp_mod_montgomery(&test2_fp, &test2_fp);

    if (fp_cmp(&test1_fp, &test2_fp) != 0) {
      printf("failed!\n\n");
      fp_printf("", &test1_fp);
      fp_printf("\n", &test2_fp);
      printf("\n\n");
      return 1;
    }
  }
  printf("fp sub.               : %.7f[ms]\n", sub_time / fp_n);
  printf("fp sub nonmod single. : %.7f[ms]\n", sub_nonmod_single / fp_n);
  printf("fp sub nonmod double. : %.7f[ms]\n", sub_nonmod_double / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_mul test\n");
  mul_time = 0, mul_nonmod_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random_montgomery(&A_fp, state);
    fp_set_random_montgomery(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_mulmod_montgomery(&test1_fp, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    mul_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_mul_nonmod(&test2_fpd, &A_fp, &B_fp);
    gettimeofday(&tv_B, NULL);
    mul_nonmod_time += timedifference_msec(tv_A, tv_B) / n;
    fp_mod_montgomery(&test2_fp, &test2_fp);

    // if(fp_cmp(&test1_fp,&test2_fp)!=0){
    //     printf("failed!\n\n");
    // 	fp_printf("",&test1_fp);
    //     fp_printf("\n",&test2_fp);
    //     printf("\n\n");
    //     return 1;
    // }
  }

  printf("fp mulmod montgomery.      : %.6f[ms]\n", mul_time / fp_n);
  printf("fp mul nonmod.             : %.6f[ms]\n", mul_nonmod_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_sqr test\n");
  sqr_time = 0, sqr_nonmod = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sqrmod_montgomery(&test1_fp, &A_fp);
    gettimeofday(&tv_B, NULL);
    sqr_time += timedifference_msec(tv_A, tv_B) / n;

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_sqr_nonmod(&test2_fpd, &A_fp);
    gettimeofday(&tv_B, NULL);
    sqr_nonmod += timedifference_msec(tv_A, tv_B) / n;
    fp_mod(&test2_fp, test2_mpn, FPLIMB2);

    // if(fp_cmp(&test1_fp,&test2_fp)!=0){
    //     printf("failed!\n\n");
    // 	fp_printf("",&test1_fp);
    //     fp_printf("\n",&test2_fp);
    //     printf("\n\n");
    //     return 1;
    // }
  }
  printf("fp sqrmod montgomery.      : %.6f[ms]\n", sqr_time / fp_n);
  printf("fp sqr nonmod.             : %.6f[ms]\n", sqr_nonmod / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_inv test\n");
  inv_time = 0, inv_lazy_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) fp_inv_montgomery(&test1_fp, &A_fp);
    gettimeofday(&tv_B, NULL);
    inv_time += timedifference_msec(tv_A, tv_B) / n;
  }
  printf("fp inv montgomery.      : %.6f[ms]\n", inv_time / fp_n);

  printf("------------------------------------------------------------------------------------\n");
  printf("fp_mod test\n");
  mod_time = 0;
  n = 1000;
  for (i = 0; i < fp_n; i++) {
    fp_set_random(&A_fp, state);
    fp_set_random(&B_fp, state);
    fp_mul_nonmod(&A_fpd, &A_fp, &B_fp);

    gettimeofday(&tv_A, NULL);
    for (j = 0; j < n; j++) mpn_mod_montgomery(A_fp.x0, FPLIMB, A_fpd.x0, FPLIMB2);

    gettimeofday(&tv_B, NULL);
    mod_time += timedifference_msec(tv_A, tv_B) / n;
  }
  printf("fp mod montgomery.      : %.6f[ms]\n", mod_time / fp_n);

  printf("fp rate\n");
  printf("fp add.               : %.7f\n", (add_time / fp_n) / (add_time / fp_n));
  printf("fp add nonmod.        : %.7f\n", (add_nonmod_single / fp_n) / (add_time / fp_n));
  printf("fp add nonmod double. : %.7f\n", (add_nonmod_double / fp_n) / (add_time / fp_n));

  printf("fp sub.               : %.7f\n", (sub_time / fp_n) / (add_time / fp_n));
  printf("fp sub nonmod.        : %.7f\n", (sub_nonmod_single / fp_n) / (add_time / fp_n));
  printf("fp sub nonmod double. : %.7f\n", (sub_nonmod_double / fp_n) / (add_time / fp_n));

  printf("fp mul mod.           : %.7f\n", (mul_time / fp_n) / (add_time / fp_n));
  printf("fp mul nonmod.        : %.7f[ms]\n", (mul_nonmod_time / fp_n) / (add_time / fp_n));
  printf("fp sqr mod.           : %.7f[ms]\n", (sqr_time / fp_n) / (add_time / fp_n));
  printf("fp sqr nonmod.        : %.7f[ms]\n", (sqr_nonmod / fp_n) / (add_time / fp_n));
  printf("fp inv montgomery.    : %.7f[ms]\n", (inv_time / fp_n) / (add_time / fp_n));
  printf("fp mod montgomery.    : %.7f[ms]\n", (mod_time / fp_n) / (add_time / fp_n));

  return 0;
}

void check_fp_with_montgomery(){
  printf("********************CHECK FP WITH MONTGOMERY*************************************************\n\n");
  printf("check_fp_with_montogomery() start...\n");
  fp_t A, ANS;
  fp_t Am, ANSm;

  fp_init(&A);
  fp_init(&ANS);
  fp_init(&Am);
  fp_init(&ANSm);
  printf("---------------------------------\n");

  fp_set_random(&A,state);
  fp_println("A =  ",&A);
  fp_to_montgomery(&Am, &A);
  fp_println_montgomery("Am = ",&Am);

  fp_inv(&ANS,&A);
  fp_inv_montgomery(&ANSm,&Am);

  fp_println("\nA^-1 =  ",&ANS);
  fp_println_montgomery("Am^-1 = ",&ANSm);

  fp_mul(&ANS,&ANS,&A);
  fp_mulmod_montgomery(&ANSm,&ANSm,&Am);

  fp_println("\nA * A^-1 =  ",&ANS);
  fp_println_montgomery("Am * Am^-1 = ",&ANSm);
  printf("---------------------------------\n");

  printf("Check sqrt algorithm\n");
  int flag=fp_legendre(&A);
  printf("fp_legendre(A) = %d\n",flag);
  if(flag==1){
    fp_sqrt(&ANS,&A);   
    fp_println("fp_sqrt(A) = ",&ANS);
    fp_to_montgomery(&ANSm, &ANS);
    fp_println_montgomery("fp_sqrt(Am) = ",&ANSm);

    fp_sqr(&ANS,&ANS);
    if(fp_cmp(&ANS,&A)==0){
      printf("(fp_sqrt(A))^2 = A\n\n");
    }
    else  printf("(fp_sqrt(A))^2 != A\n\n");

    fp_sqrmod_montgomery(&ANSm,&ANSm);
    if(fp_cmp(&ANSm,&Am)==0){
      printf("(fp_sqrt(Am))^2 = A\n\n");
    }
    else  printf("(fp_sqrt(Am))^2 != A\n\n");

  }
  printf("---------------------------------\n");

  printf("Check Fermert's little theorem\n");
  mpz_t tmp;
  mpz_init(tmp);
  mpz_sub_ui(tmp,prime_z,1);
  fp_pow(&ANS,&A,tmp);
  fp_println("A^(p-1) =  ",&ANS);

  fp_pow_montgomery(&ANSm,&Am,tmp);
  fp_println_montgomery("Am^(p-1) = ",&ANSm);
  printf("---------------------------------\n");

  fp_println("\nA =  ",&A);
  fp_println_montgomery("Am = ",&Am);
  printf("---------------------------------\n");
  mpz_clear(tmp);

  printf("*********************************************************************************************\n\n");
}

void check_fp2_with_montgomery(){
  printf("********************CHECK FP2 WITH MONTGOMERY*************************************************\n\n");
  printf("check_fp2_with_montogomery() start...\n");
  fp2_t A, ANS;
  fp2_t Am, ANSm;

  fp2_init(&A);
  fp2_init(&ANS);
  fp2_init(&Am);
  fp2_init(&ANSm);
  printf("---------------------------------\n");

  fp2_set_random(&A,state);
  fp2_println("A =  ",&A);
  fp2_to_montgomery(&Am, &A);
  fp2_println_montgomery("Am = ",&Am);

  fp2_inv(&ANS,&A);
  fp2_inv_lazy_montgomery(&ANSm,&Am);

  fp2_println("\nA^-1 =  ",&ANS);
  fp2_println_montgomery("Am^-1 = ",&ANSm);

  fp2_mul(&ANS,&ANS,&A);
  fp2_mul_lazy_montgomery(&ANSm,&ANSm,&Am);

  fp2_println("\nA * A^-1 =  ",&ANS);
  fp2_println_montgomery("Am * Am^-1 = ",&ANSm);
  printf("---------------------------------\n");

  printf("Check sqrt algorithm\n");
  int flag=fp2_legendre(&A);
  printf("fp_legendre(A) = %d\n",flag);
  if(flag==1){
    fp2_sqrt(&ANS,&A);
    fp2_println("fp_sqrt(A) = ",&ANS);
    fp2_to_montgomery(&ANSm, &ANS);
    fp2_println_montgomery("fp_sqrt(Am) = ",&ANSm);

    fp2_sqr(&ANS,&ANS);
    if(fp2_cmp(&ANS,&A)==0){
      printf("(fp_sqrt(A))^2 = A\n\n");
    }
    else  printf("(fp_sqrt(A))^2 != A\n\n");

    fp2_sqr_lazy_montgomery(&ANSm,&ANSm);
    if(fp2_cmp(&ANSm,&Am)==0){
      printf("(fp_sqrt(Am))^2 = A\n\n");
    }
    else  printf("(fp_sqrt(Am))^2 != A\n\n");

  }
  printf("---------------------------------\n");

  printf("Check Fermert's little theorem\n");
  mpz_t tmp;
  mpz_init(tmp);
  mpz_pow_ui(tmp,prime_z,2);
  mpz_sub_ui(tmp,tmp,1);
  fp2_pow(&ANS,&A,tmp);
  fp2_println("A^(p^2-1) =  ",&ANS);

  fp2_pow_montgomery(&ANSm,&Am,tmp);
  fp2_println_montgomery("Am^(p^2-1) = ",&ANSm);
  printf("---------------------------------\n");

  fp2_println("\nA =  ",&A);
  fp2_println_montgomery("Am = ",&Am);
  printf("---------------------------------\n");
  printf("*********************************************************************************************\n\n");
}


void check_fp4_with_montgomery(){
  printf("********************CHECK FP4 WITH MONTGOMERY*************************************************\n\n");
  printf("check_fp4_with_montogomery() start...\n");
  fp4_t A, ANS;
  fp4_t Am, ANSm;

  fp4_init(&A);
  fp4_init(&ANS);
  fp4_init(&Am);
  fp4_init(&ANSm);
  printf("---------------------------------\n");

  fp4_set_random(&A,state);
  fp4_println("A =  ",&A);
  fp4_to_montgomery(&Am, &A);
  fp4_println_montgomery("Am = ",&Am);

  fp4_inv(&ANS,&A);
  fp4_inv_lazy_montgomery(&ANSm,&Am);

  fp4_println("\nA^-1 =  ",&ANS);
  fp4_println_montgomery("Am^-1 = ",&ANSm);

  fp4_mul(&ANS,&ANS,&A);
  fp4_mul_lazy_montgomery(&ANSm,&ANSm,&Am);

  fp4_println("\nA * A^-1 =  ",&ANS);
  fp4_println_montgomery("Am * Am^-1 = ",&ANSm);
  printf("---------------------------------\n");

  printf("Check sqrt algorithm\n");
  int flag=fp4_legendre(&A);
  printf("fp4_legendre(A) = %d\n",flag);
  if(flag==1){
    fp4_sqrt(&ANS,&A);
    fp4_println("fp4_sqrt(A) = ",&ANS);
    fp4_to_montgomery(&ANSm, &ANS);
    fp4_println_montgomery("fp4_sqrt(Am) = ",&ANSm);

    fp4_sqr(&ANS,&ANS);
    if(fp4_cmp(&ANS,&A)==0){
      printf("(fp4_sqrt(A))^4 = A\n\n");
    }
    else  printf("(fp4_sqrt(A))^4 != A\n\n");

    fp4_sqr_lazy_montgomery(&ANSm,&ANSm);
    if(fp4_cmp(&ANSm,&Am)==0){
      printf("(fp4_sqrt(Am))^4 = A\n\n");
    }
    else  printf("(fp4_sqrt(Am))^4 != A\n\n");

  }
  printf("---------------------------------\n");

  printf("Check Fermert's little theorem\n");
  mpz_t tmp;
  mpz_init(tmp);
  mpz_pow_ui(tmp,prime_z,4);
  mpz_sub_ui(tmp,tmp,1);
  fp4_pow(&ANS,&A,tmp);
  fp4_println("A^(p^4-1) =  ",&ANS);

  fp4_pow_montgomery(&ANSm,&Am,tmp);
  fp4_println_montgomery("Am^(p^4-1) = ",&ANSm);
  printf("---------------------------------\n");

  fp4_println("\nA =  ",&A);
  fp4_println_montgomery("Am = ",&Am);
  printf("---------------------------------\n");
  printf("*********************************************************************************************\n\n");
}

void check_fp8_with_montgomery(){
  printf("********************CHECK FP8 WITH MONTGOMERY*************************************************\n\n");
  printf("check_fp8_with_montogomery() start...\n");
  fp8_t A, ANS;
  fp8_t Am, ANSm;

  fp8_init(&A);
  fp8_init(&ANS);
  fp8_init(&Am);
  fp8_init(&ANSm);
  printf("---------------------------------\n");

  fp8_set_random(&A,state);
  fp8_println("A =  ",&A);
  fp8_to_montgomery(&Am, &A);
  fp8_println_montgomery("Am = ",&Am);

  fp8_inv(&ANS,&A);
  fp8_inv_lazy_montgomery(&ANSm,&Am);

  fp8_println("\nA^-1 =  ",&ANS);
  fp8_println_montgomery("Am^-1 = ",&ANSm);

  fp8_mul(&ANS,&ANS,&A);
  fp8_mul_lazy_montgomery(&ANSm,&ANSm,&Am);

  fp8_println("\nA * A^-1 =  ",&ANS);
  fp8_println_montgomery("Am * Am^-1 = ",&ANSm);
  printf("---------------------------------\n");

  printf("Check sqrt algorithm\n");
  int flag=fp8_legendre(&A);
  printf("fp8_legendre(A) = %d\n",flag);
  if(flag==1){
    fp8_sqrt(&ANS,&A);
    fp8_println("fp8_sqrt(A) = ",&ANS);
    fp8_to_montgomery(&ANSm, &ANS);
    fp8_println_montgomery("fp8_sqrt(Am) = ",&ANSm);

    fp8_sqr(&ANS,&ANS);
    if(fp8_cmp(&ANS,&A)==0){
      printf("(fp8_sqrt(A))^8 = A\n\n");
    }
    else  printf("(fp8_sqrt(A))^8 != A\n\n");

    fp8_sqr_lazy_montgomery(&ANSm,&ANSm);
    if(fp8_cmp(&ANSm,&Am)==0){
      printf("(fp8_sqrt(Am))^8 = A\n\n");
    }
    else  printf("(fp8_sqrt(Am))^8 != A\n\n");

  }
  printf("---------------------------------\n");

  printf("Check Fermert's little theorem\n");
  mpz_t tmp;
  mpz_init(tmp);
  mpz_pow_ui(tmp,prime_z,8);
  mpz_sub_ui(tmp,tmp,1);
  fp8_pow(&ANS,&A,tmp);
  fp8_println("A^(p^8-1) =  ",&ANS);

  fp8_pow_montgomery(&ANSm,&Am,tmp);
  fp8_println_montgomery("Am^(p^8-1) = ",&ANSm);
  printf("---------------------------------\n");

  fp8_println("\nA =  ",&A);
  fp8_println_montgomery("Am = ",&Am);
  printf("---------------------------------\n");
  printf("*********************************************************************************************\n\n");
}

void BENCH_fp2_fp4_fp8_mul_lazy_montgomery(int LOOP){
  printf("============================================================================\n");
  printf("-------------------------------fp, fp2, fp6 BENCHMARK-----------------------\n");
  printf("============================================================================\n");

  fp_t ANS_fp, ANSm_fp, A_fp,Am_fp ,B_fp,Bm_fp;
  fp_init(&A_fp);
  fp_init(&Am_fp);
  fp_init(&B_fp);
  fp_init(&Bm_fp);
  fp_init(&ANS_fp);
  fp_init(&ANSm_fp);
  fp_set_random(&A_fp,state);
  fp_to_montgomery(&Am_fp, &A_fp);
  fp_set_random(&B_fp,state);
  fp_to_montgomery(&Bm_fp, &B_fp);

  CYBOZU_BENCH_C("fp_mul()", LOOP, fp_mul,&ANS_fp,&A_fp,&B_fp);
  CYBOZU_BENCH_C("fp_mul_lazy_montgomery()", LOOP, fp_mulmod_montgomery, &ANSm_fp,&Am_fp,&Bm_fp);
  fp_println("ANS:", &ANS_fp);
  fp_println_montgomery("ANSm:", &ANSm_fp);

  CYBOZU_BENCH_C("fp_sqr()", LOOP, fp_sqr,&ANS_fp,&A_fp);
  CYBOZU_BENCH_C("fp_sqr_lazy_montgomery()", LOOP, fp_sqrmod_montgomery, &ANSm_fp,&Am_fp);
  fp_println("ANS:", &ANS_fp);
  fp_println_montgomery("ANSm:", &ANSm_fp);

  fp2_t ANS_fp2, ANSm_fp2, A_fp2,Am_fp2 ,B_fp2,Bm_fp2;
  fp2_init(&A_fp2);
  fp2_init(&Am_fp2);
  fp2_init(&B_fp2);
  fp2_init(&Bm_fp2);
  fp2_init(&ANS_fp2);
  fp2_init(&ANSm_fp2);

  fp2_set_random(&A_fp2,state);
  fp2_to_montgomery(&Am_fp2, &A_fp2);
  fp2_set_random(&B_fp2,state);
  fp2_to_montgomery(&Bm_fp2, &B_fp2);

  CYBOZU_BENCH_C("fp2_mul()", LOOP, fp2_mul,&ANS_fp2,&A_fp2,&B_fp2);
  CYBOZU_BENCH_C("fp2_mul_lazy_montgomery()", LOOP, fp2_mul_lazy_montgomery, &ANSm_fp2,&Am_fp2,&Bm_fp2);
  fp2_println("ANS:", &ANS_fp2);
  fp2_println_montgomery("ANSm:", &ANSm_fp2);

  CYBOZU_BENCH_C("fp2_sqr()", LOOP, fp2_sqr,&ANS_fp2,&A_fp2);
  CYBOZU_BENCH_C("fp2_sqr_lazy_montgomery()", LOOP, fp2_sqr_lazy_montgomery, &ANSm_fp2,&Am_fp2);
  fp2_println("ANS:", &ANS_fp2);
  fp2_println_montgomery("ANSm:", &ANSm_fp2);

  fp4_t ANS_fp4, ANSm_fp4, A_fp4,Am_fp4 ,B_fp4,Bm_fp4;
  fp4_init(&A_fp4);
  fp4_init(&Am_fp4);
  fp4_init(&B_fp4);
  fp4_init(&Bm_fp4);
  fp4_init(&ANS_fp4);
  fp4_init(&ANSm_fp4);

  fp4_set_random(&A_fp4,state);
  fp4_to_montgomery(&Am_fp4, &A_fp4);
  fp4_set_random(&B_fp4,state);
  fp4_to_montgomery(&Bm_fp4, &B_fp4);

  CYBOZU_BENCH_C("fp4_mul()", LOOP, fp4_mul,&ANS_fp4,&A_fp4,&B_fp4);
  CYBOZU_BENCH_C("fp4_mul_lazy_montgomery()", LOOP, fp4_mul_lazy_montgomery, &ANSm_fp4,&Am_fp4,&Bm_fp4);
  fp4_println("ANS:", &ANS_fp4);
  fp4_println_montgomery("ANSm:", &ANSm_fp4);

  CYBOZU_BENCH_C("fp4_sqr()", LOOP, fp4_sqr,&ANS_fp4,&A_fp4);
  CYBOZU_BENCH_C("fp4_sqr_lazy_montgomery()", LOOP, fp4_sqr_lazy_montgomery, &ANSm_fp4,&Am_fp4);
  fp4_println("ANS:", &ANS_fp4);
  fp4_println_montgomery("ANSm:", &ANSm_fp4);

  fp8_t ANS_fp8, ANSm_fp8, A_fp8,Am_fp8 ,B_fp8,Bm_fp8;
  fp8_init(&A_fp8);
  fp8_init(&Am_fp8);
  fp8_init(&B_fp8);
  fp8_init(&Bm_fp8);
  fp8_init(&ANS_fp8);
  fp8_init(&ANSm_fp8);

  fp8_set_random(&A_fp8,state);
  fp8_to_montgomery(&Am_fp8, &A_fp8);
  fp8_set_random(&B_fp8,state);
  fp8_to_montgomery(&Bm_fp8, &B_fp8);

  CYBOZU_BENCH_C("fp8_mul()", LOOP, fp8_mul,&ANS_fp8,&A_fp8,&B_fp8);
  CYBOZU_BENCH_C("fp8_mul_lazy_montgomery()", LOOP, fp8_mul_lazy_montgomery, &ANSm_fp8,&Am_fp8,&Bm_fp8);
  fp8_println("ANS:", &ANS_fp8);
  fp8_println_montgomery("ANSm:", &ANSm_fp8);

  CYBOZU_BENCH_C("fp8_sqr()", LOOP, fp8_sqr,&ANS_fp8,&A_fp8);
  CYBOZU_BENCH_C("fp8_sqr_lazy_montgomery()", LOOP, fp8_sqr_lazy_montgomery, &ANSm_fp8,&Am_fp8);
  fp8_println("ANS:", &ANS_fp8);
  fp8_println_montgomery("ANSm:", &ANSm_fp8);

  CYBOZU_BENCH_C("fp8_sqr_GS()", LOOP, fp8_sqr_GS,&ANS_fp8,&A_fp8);
  CYBOZU_BENCH_C("fp8_sqr_GS_lazy_montgomery()", LOOP, fp8_sqr_GS_lazy_montgomery, &ANSm_fp8,&Am_fp8);
  fp8_println("ANS:", &ANS_fp8);
  fp8_println_montgomery("ANSm:", &ANSm_fp8);
  
  CYBOZU_BENCH_C("fp8_mul_sparse_add()", LOOP, fp8_mul_sparse_add,&ANS_fp8,&A_fp8,&B_fp8);
  CYBOZU_BENCH_C("fp8_mul_sparse_add_lazy_montgomery()", LOOP, fp8_mul_sparse_add_lazy_montgomery, &ANSm_fp8,&Am_fp8,&Bm_fp8);
  fp8_println("ANS:", &ANS_fp8);
  fp8_println_montgomery("ANSm:", &ANSm_fp8);

  CYBOZU_BENCH_C("fp8_mul_sparse_dbl()", LOOP, fp8_mul_sparse_dbl,&ANS_fp8,&A_fp8,&B_fp8);
  CYBOZU_BENCH_C("fp8_mul_sparse_dbl_lazy_montgomery()", LOOP, fp8_mul_sparse_dbl_lazy_montgomery, &ANSm_fp8,&Am_fp8,&Bm_fp8);
  fp8_println("ANS:", &ANS_fp8);
  fp8_println_montgomery("ANSm:", &ANSm_fp8);

}

void BENCH_miller_lazy_montgomery(int LOOP){
  printf("============================================================================\n");
  printf("--------------------------------------ltpq, lttp----------------------------\n");
  printf("============================================================================\n");

  efp8_t P,Q;
  fp8_t f,e1,e2;
  fp8_t fm,e1m,e2m;

  mpz_t a,b,ab;
  efp8_init(&P);
  efp8_init(&Q);

  fp8_init(&f);
  fp8_init(&e1);
  fp8_init(&e2);
  fp8_init(&fm);
  fp8_init(&e1m);
  fp8_init(&e2m);
  mpz_init(a);
  mpz_init(b);
  mpz_init(ab);

  generate_g1(&P);
  generate_g2(&Q);

  static efp_t mapped_P;
  static efp2_t mapped_Q,mapped_Q_neg;
  static efp_t mapped_Pm;
  static efp2_t mapped_Qm,mapped_Q_negm;
  static efp2_jacobian_t S;
  static efp2_jacobian_t Sm;
//------------------Regular-----------------------
  fp8_set_ui_ui(&f,0);
  fp8_set_ui(&f,1);

  fp_set(&mapped_P.x,&P.x.x0.x0.x0);
  fp_set(&mapped_P.y,&P.y.x0.x0.x0);
  mapped_P.infinity = 0;

  efp8_to_efp2(&mapped_Q,&Q);//twist
  efp8_to_Jacefp2(&S,&Q);
  efp2_set_neg(&mapped_Q_neg,&mapped_Q);
//------------------Montgomery-----------------------
  fp8_set_ui_ui(&fm,0);
  fp8_set_ui(&fm,1);
  fp_to_montgomery(&fm.x0.x0.x0, &fm.x0.x0.x0);

  fp_to_montgomery(&mapped_Pm.x,&P.x.x0.x0.x0);
  fp_to_montgomery(&mapped_Pm.y,&P.y.x0.x0.x0);
  mapped_Pm.infinity = 0;

  efp8_to_efp2_montgomery(&mapped_Qm,&Q);//twist
  efp8_to_Jacefp2_montgomery(&Sm,&Q);
  efp2_set_neg_montgomery(&mapped_Q_negm,&mapped_Q);//here
//-----------------------------------------

  CYBOZU_BENCH_C("ff_lttp()", LOOP, ff_lttp,&f,&S,&mapped_P);
  CYBOZU_BENCH_C("ff_lttp_lazy_montgomery()", LOOP, ff_lttp_lazy_montgomery, &fm,&Sm,&mapped_Pm);
  fp8_println("ANS:", &f);
  fp8_println_montgomery("ANSm:", &fm);
  fp8_set_ui_ui(&f,0);
  fp8_set_ui(&f,1);
  fp8_set_ui_ui(&fm,0);
  fp8_set_ui(&fm,1);
  fp_to_montgomery(&fm.x0.x0.x0, &fm.x0.x0.x0);

  CYBOZU_BENCH_C("ff_ltpq()", LOOP, ff_ltqp,&f,&S,&mapped_Q,&mapped_P);
  CYBOZU_BENCH_C("ff_ltpq_lazy_montgomery()", LOOP, ff_ltqp_lazy_montgomery, &fm,&Sm,&mapped_Qm,&mapped_Pm);
  fp8_println("ANS:", &f);
  fp8_println_montgomery("ANSm:", &fm);

  printf("============================================================================\n");
  printf("--------------------------------------Miller--------------------------------\n");
  printf("============================================================================\n");


  CYBOZU_BENCH_C("miller_opt_ate_jac_2NAF()", LOOP, miller_opt_ate_jac_2NAF,&f,&P,&Q);
  fp8_println("ANS:", &f);
  
  CYBOZU_BENCH_C("miller_opt_ate_jac_2NAF_lazy_montgomery()", LOOP, miller_opt_ate_jac_2NAF_lazy_montgomery,&fm,&P,&Q);
  fp8_println_montgomery("ANSm:", &fm);
}

void BENCH_miller_coordinates(int LOOP){
  printf("============================================================================\n");
  printf("--------------------------------------ltpq, lttp----------------------------\n");
  printf("============================================================================\n");

  efp8_t P,Q;
  fp8_t f,e1,e2;

  mpz_t a,b,ab;
  efp8_init(&P);
  efp8_init(&Q);

  fp8_init(&f);
  fp8_init(&e1);
  fp8_init(&e2);

  mpz_init(a);
  mpz_init(b);
  mpz_init(ab);

  generate_g1(&P);
  generate_g2(&Q);

  static efp_t mapped_P;
  static efp2_t mapped_Q,mapped_Q_neg;

  static efp2_jacobian_t S;
  static efp2_jacobian_t Sp;

  fp8_set_ui_ui(&f,0);
  fp8_set_ui(&f,1);

  fp_set(&mapped_P.x,&P.x.x0.x0.x0);
  fp_set(&mapped_P.y,&P.y.x0.x0.x0);
  mapped_P.infinity = 0;

  efp8_to_efp2(&mapped_Q,&Q);//twist
  efp8_to_Jacefp2(&S,&Q);
  efp8_to_Jacefp2(&Sp,&Q);
  efp2_set_neg(&mapped_Q_neg,&mapped_Q);
  // miller_proj_precomp_Costello(&mapped_P, &mapped_Q);
//-----------------------------------------

  CYBOZU_BENCH_C("ff_lttp()", LOOP, ff_lttp,&f,&S,&mapped_P);
  fp8_set_ui(&f,1);
  CYBOZU_BENCH_C("ff_lttp_Costello()", LOOP, ff_lttp_Costello, &f,&Sp,&mapped_P);
  fp8_println("ANS:", &f);
  fp8_println_montgomery("ANS:", &f);

  fp8_set_ui(&f,1);
  CYBOZU_BENCH_C("ff_ltpq()", LOOP, ff_ltqp,&f,&S,&mapped_Q,&mapped_P);
  fp8_set_ui(&f,1);
  CYBOZU_BENCH_C("ff_ltqp_Costello_mixed()", LOOP, ff_ltqp_Costello_mixed, &f,&Sp,&mapped_Q,&mapped_P);
  fp8_println("ANS:", &f);
  fp8_println_montgomery("ANS:", &f);

  printf("============================================================================\n");
  printf("--------------------------------------Miller--------------------------------\n");
  printf("============================================================================\n");

//------------------Jacobi(mixed)-----------------------

  CYBOZU_BENCH_C("miller_opt_ate_jac_2NAF()", LOOP, miller_opt_ate_jac_2NAF,&f,&P,&Q);
  fp8_println("ANS:", &f);
  //------------------Projective*[1:2](mixed)-----------------------
  fp8_set_ui(&f,1);

  CYBOZU_BENCH_C("miller_opt_ate_proj_2NAF()", LOOP, miller_opt_ate_proj_2NAF,&f,&P,&Q);
  fp8_println_montgomery("ANS:", &f);
}



void BENCH_finalexp_lazy_montgomery(int LOOP){
  printf("********************CHECK Finalexp WITH MONTGOMERY*************************************************\n\n");
  printf("check_finalexp_with_montogomery() start...\n");
  fp8_t A, ANS;
  fp8_t Am, ANSm;

  fp8_init(&A);
  fp8_init(&ANS);
  fp8_init(&Am);
  fp8_init(&ANSm);

  fp8_set_random(&A, state);
  fp8_to_montgomery(&Am, &A);

  CYBOZU_BENCH_C("fp8_sqr_GS()", LOOP, fp8_sqr_GS,&ANS, &A);
  CYBOZU_BENCH_C("fp8_sqr_GS_lazy_montgomery()", LOOP, fp8_sqr_GS_lazy_montgomery, &ANSm, &Am);

  fp8_println("fp8_sqr_GS\n", &ANS);
  fp8_println_montgomery("fp8_sqr_GS_lazy_montgomery\n", &ANSm);

  CYBOZU_BENCH_C("fp8_finalexpow_x_2NAF()", LOOP, fp8_finalexpow_x_2NAF,&ANS, &A);
  CYBOZU_BENCH_C("fp8_finalexpow_x_2NAF_lazy_montgomery()", LOOP, fp8_finalexpow_x_2NAF_lazy_montgomery, &ANSm, &Am);

  fp8_println("fp8_finalexpow_x_2NAF\n", &ANS);
  fp8_println_montgomery("fp8_finalexpow_x_2NAF_lazy_montgomery\n", &ANSm);

  CYBOZU_BENCH_C("final_exp()", LOOP, final_exp,&ANS, &A);
  CYBOZU_BENCH_C("final_exp_lazy_montgomery()", LOOP, final_exp_lazy_montgomery, &ANSm, &Am);

  fp8_println("final_exp\n", &ANS);
  fp8_println_montgomery("final_exp_lazy_montgomery\n", &ANSm);

}

void BENCH_Pairing_jac_lazy_montgomery(int LOOP){
  printf("check_pairing_jac() 開始\n");
  efp8_t P,Q,aP,bQ,tmp1;
  fp8_t f,e1,e2;
  mpz_t a,b,ab;
  efp8_init(&P);
  efp8_init(&Q);

  efp8_init(&aP);
  efp8_init(&bQ);
  efp8_init(&tmp1);
  fp8_init(&f);
  fp8_init(&e1);
  fp8_init(&e2);
  mpz_init(a);
  mpz_init(b);
  mpz_init(ab);

  generate_g1(&P);
  generate_g2(&Q);

  mpz_urandomm(a,state,prime_z);
  mpz_urandomm(b,state,prime_z);

  #if 1
  efp8_println("P = ",&P);
  efp8_println("Q = ",&Q);

  efp8_scm_fp(&tmp1,&P,order_z);
  efp8_println("[r]P = ",&tmp1);
  efp8_scm(&tmp1,&Q,order_z);
  efp8_println("[r]Q = ",&tmp1);

  gmp_printf("a = %Zd\n",a);
  gmp_printf("b = %Zd\n",b);
  printf("---------------------------------\n");
  #endif

  printf("---------------------------------\n");
  printf("check regular pairing()\n");
  printf("---------------------------------\n");
  //e([a]P,[b]Q) を求める
  efp8_scm_fp(&aP,&P,a);
  efp8_scm(&bQ,&Q,b);
  miller_opt_ate_jac_2NAF(&f,&aP,&bQ);
  final_exp(&e1,&f);
  //e(P,Q)^(a*b) を求める
  miller_opt_ate_jac_2NAF(&f,&P,&Q);
  final_exp(&e2,&f);
  mpz_mul(ab,a,b);
  fp8_pow(&e2,&e2,ab);
  fp8_println("e([a]P,[b]Q) = ",&e1);
  fp8_println("e(P,Q)^(a*b) = ",&e2);
  if(fp8_cmp(&e1,&e2)==0)  {
  printf("=====================================================\n");
  printf("------------------bilinear!!-------------------------\n");
  printf("=====================================================\n");
  }else{
    printf("e([a]P,[b]Q) != e(P,Q)^(a*b)\n\n");
  }

  printf("---------------------------------\n");
  printf("check lazy montgomery pairing()\n");
  printf("---------------------------------\n");
  miller_opt_ate_jac_2NAF_lazy_montgomery(&f,&aP,&bQ);
  final_exp_lazy_montgomery(&e1,&f);
  //e(P,Q)^(a*b) を求める
  miller_opt_ate_jac_2NAF_lazy_montgomery(&f,&P,&Q);
  final_exp_lazy_montgomery(&e2,&f);

  fp8_pow_montgomery(&e2,&e2,ab);
  fp8_println_montgomery("e([a]P,[b]Q) = ",&e1);
  fp8_println_montgomery("e(P,Q)^(a*b) = ",&e2);
  if(fp8_cmp(&e1,&e2)==0)  {
  printf("=====================================================\n");
  printf("------------------bilinear!!-------------------------\n");
  printf("=====================================================\n\n\n");
  }else{
    printf("e([a]P,[b]Q) != e(P,Q)^(a*b)\n\n");
  }

  // printf("--------Benching pairing()-------\n");

  // CYBOZU_BENCH_C("miller_opt_ate_proj_2NAF()", LOOP, miller_opt_ate_proj_2NAF,&f,&P,&Q);
  // CYBOZU_BENCH_C("final_exp()               ", LOOP, final_exp,&e2, &e1);
  // printf("---------------------------------\n");

  CYBOZU_BENCH_C("miller_opt_ate_jac_2NAF_lazy_montgomery()", LOOP, miller_opt_ate_jac_2NAF_lazy_montgomery,&f,&P,&Q);
  CYBOZU_BENCH_C("final_exp_lazy_montgomery()               ", LOOP, final_exp_lazy_montgomery, &e2, &e1);
  printf("---------------------------------\n");

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(ab);

  printf("*********************************************************************************************\n\n");
}


void BENCH_Pairing_proj_lazy_montgomery(int LOOP){
  printf("check_pairing_proj() 開始\n");
  efp8_t P,Q,aP,bQ,tmp1;
  fp8_t f,e1,e2;
  mpz_t a,b,ab;
  efp8_init(&P);
  efp8_init(&Q);

  efp8_init(&aP);
  efp8_init(&bQ);
  efp8_init(&tmp1);
  fp8_init(&f);
  fp8_init(&e1);
  fp8_init(&e2);
  mpz_init(a);
  mpz_init(b);
  mpz_init(ab);

  generate_g1(&P);
  generate_g2(&Q);

  mpz_urandomm(a,state,prime_z);
  mpz_urandomm(b,state,prime_z);

  #if 1
  efp8_println("P = ",&P);
  efp8_println("Q = ",&Q);

  efp8_scm_fp(&tmp1,&P,order_z);
  efp8_println("[r]P = ",&tmp1);
  efp8_scm(&tmp1,&Q,order_z);
  efp8_println("[r]Q = ",&tmp1);

  gmp_printf("a = %Zd\n",a);
  gmp_printf("b = %Zd\n",b);
  printf("---------------------------------\n");
  #endif

  printf("---------------------------------\n");
  printf("check regular proj pairing()\n");
  printf("---------------------------------\n");
  //e([a]P,[b]Q) を求める
  efp8_scm_fp(&aP,&P,a);
  efp8_scm(&bQ,&Q,b);
  miller_opt_ate_proj_2NAF(&f,&aP,&bQ);
  final_exp(&e1,&f);
  //e(P,Q)^(a*b) を求める
  miller_opt_ate_proj_2NAF(&f,&P,&Q);
  final_exp(&e2,&f);
  mpz_mul(ab,a,b);
  fp8_pow(&e2,&e2,ab);
  fp8_println("e([a]P,[b]Q) = ",&e1);
  fp8_println("e(P,Q)^(a*b) = ",&e2);
  if(fp8_cmp(&e1,&e2)==0)  {
  printf("=====================================================\n");
  printf("------------------bilinear!!-------------------------\n");
  printf("=====================================================\n");
  }else{
    printf("e([a]P,[b]Q) != e(P,Q)^(a*b)\n\n");
  }

  printf("---------------------------------\n");
  printf("check lazy montgomery pairing()\n");
  printf("---------------------------------\n");
  miller_opt_ate_proj_2NAF_lazy_montgomery(&f,&aP,&bQ);
  final_exp_lazy_montgomery(&e1,&f);
  //e(P,Q)^(a*b) を求める
  miller_opt_ate_proj_2NAF_lazy_montgomery(&f,&P,&Q);
  final_exp_lazy_montgomery(&e2,&f);

  fp8_pow_montgomery(&e2,&e2,ab);
  fp8_println_montgomery("e([a]P,[b]Q) = ",&e1);
  fp8_println_montgomery("e(P,Q)^(a*b) = ",&e2);
  if(fp8_cmp(&e1,&e2)==0)  {
  printf("=====================================================\n");
  printf("------------------bilinear!!-------------------------\n");
  printf("=====================================================\n\n\n");
  }else{
    printf("e([a]P,[b]Q) != e(P,Q)^(a*b)\n\n");
  }

  // printf("--------Benching pairing()-------\n");

  // CYBOZU_BENCH_C("miller_opt_ate_proj_2NAF()", LOOP, miller_opt_ate_proj_2NAF,&f,&P,&Q);
  // CYBOZU_BENCH_C("final_exp()               ", LOOP, final_exp,&e2, &e1);
  // printf("---------------------------------\n");

  CYBOZU_BENCH_C("miller_opt_ate_proj_2NAF_lazy_montgomery()", LOOP, miller_opt_ate_proj_2NAF_lazy_montgomery,&f,&P,&Q);
  CYBOZU_BENCH_C("final_exp_lazy_montgomery()               ", LOOP, final_exp_lazy_montgomery, &e2, &e1);
  printf("---------------------------------\n");

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(ab);

  printf("*********************************************************************************************\n\n");
}


void check_finalexp_pow_cost_count_2NAF(){
  printf("check_pairing_count_2NAF() 開始\n");
  cost final_exp_x_cost, final_exp_x_2_cost,final_exp_hy_cost,final_exp_4hy_cost;

  fp8_t f;
  fp8_init(&f);
  fp8_set_random(&f, state);

  printf("fp8_finalexpow_x_2NAF count\n");
  cost_zero();
  fp8_finalexpow_x_2NAF(&f,&f);
  cost_check(&final_exp_x_cost);
  cost_printf("fp8_finalexpow_x_2NAF cost",&final_exp_x_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_finalexpow_hy_neg_2NAF count\n");
  cost_zero();
  fp8_finalexpow_hy_neg_2NAF(&f,&f);
  cost_check(&final_exp_hy_cost);
  cost_printf("fp8_finalexpow_hy_neg_2NAF cost",&final_exp_hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_finalexpow_4hy_neg_2NAF count\n");
  cost_zero();
  fp8_finalexpow_4hy_neg_2NAF(&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_finalexpow_4hy_neg_2NAF cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_mul count\n");
  cost_zero();
  fp8_mul(&f,&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_mul cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_sqr count\n");
  cost_zero();
  fp8_sqr(&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_sqr cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("*********************************************************************************************\n\n");
}

void check_finalexp_pow_cost_count_2NAF_montgomery(){
  printf("check_pairing_count_2NAF() 開始\n");
  cost final_exp_x_cost, final_exp_x_2_cost,final_exp_hy_cost,final_exp_4hy_cost;

  fp8_t f;
  fp8_init(&f);
  fp8_set_random(&f, state);
  fp8_to_montgomery(&f, &f);

  printf("fp8_finalexpow_x_2NAF count\n");
  cost_zero();
  fp8_finalexpow_x_2NAF_lazy_montgomery(&f,&f);
  cost_check(&final_exp_x_cost);
  cost_printf("fp8_finalexpow_x_2NAF cost",&final_exp_x_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_finalexpow_hy_neg_2NAF count\n");
  cost_zero();
  fp8_finalexpow_hy_neg_2NAF_lazy_montgomery(&f,&f);
  cost_check(&final_exp_hy_cost);
  cost_printf("fp8_finalexpow_hy_neg_2NAF cost",&final_exp_hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_finalexpow_4hy_neg_2NAF count\n");
  cost_zero();
  fp8_finalexpow_4hy_neg_2NAF_lazy_montgomery(&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_finalexpow_4hy_neg_2NAF cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_mul count\n");
  cost_zero();
  fp8_mul_lazy_montgomery(&f,&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_mul cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_sqr count\n");
  cost_zero();
  fp8_sqr_lazy_montgomery(&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_sqr cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("fp8_inv count\n");
  cost_zero();
  fp8_inv(&f,&f);
  cost_check(&final_exp_4hy_cost);
  cost_printf("fp8_inv cost",&final_exp_4hy_cost,CHECK_PAIRING_TIME_LOOP);
  printf("---------------------------------\n");

  printf("*********************************************************************************************\n\n");
}

void BENCH_Pairing_proj_lazy_montgomery_static(int LOOP){
  printf("check_pairing_proj_static() 開始\n");
  efp8_t P,Q,aP,bQ,tmp1;
  fp8_t f,e1,e2;
  mpz_t a,b,ab;
  mpz_t px,py,qx1,qx2,qy1,qy2;
  efp8_init(&P);
  efp8_init(&Q);

  efp8_init(&aP);
  efp8_init(&bQ);
  efp8_init(&tmp1);
  fp8_init(&f);
  fp8_init(&e1);
  fp8_init(&e2);
  mpz_init_set_str(a,"11572949769166380585226343310304350276488912480716866848234549974660437578720144552171126138682194799744294811883399580170344708389502031284500218968890503062489504",10);
  mpz_init_set_str(b,"4504182117316004236696961332504379040964341668012094821966581709623693876978478622078946935261083156616644047362379323023674883729763837056167376954281761355584822",10);
  mpz_init(ab);
  mpz_mul(ab,a,b);

  mpz_init_set_str(px,"7467710015278176793660393103569331477529841715415080697012876328859729958066492978505442439617072015466752820723220839236534637183885440781110217707638903976782227",10);
  mpz_init_set_str(py,"6088279501888905280881173198966267404997289204083625466407527005926994968981179507755569183218572146929511055540237905378670023408321875057754423935572054116115708",10);
  mpz_init_set_str(qx1,"10671655066887411328397278463188876464639543705795716877800448659470679931618354899200260167799108043301013519469714763854149286903891803007745825431260689369853284",10);
  mpz_init_set_str(qx2,"11014797048893612536343310546760419762816767069339245261942010635594949509899080685714684950612622057818175867547251423637374137130693438744056222812034522532405935",10);
  mpz_init_set_str(qy1,"6441373089628608027325282797453628876799675726313687361132764629904185179891386524497545431211574477354766519515118678089875871374414426927121416252290295677827654",10);
  mpz_init_set_str(qy2,"10015102703601587049205701816758987818504702469415710818672176397800923925736678922978002589618956690638917227779287038173369027081108414140155264387247209652322518",10);
  mpn_set_mpz(P.x.x0.x0.x0.x0, px);
  mpn_set_mpz(P.y.x0.x0.x0.x0, py);
  mpn_set_mpz(Q.x.x0.x1.x0.x0, qx1);
  mpn_set_mpz(Q.x.x0.x1.x1.x0, qx2);
  mpn_set_mpz(Q.y.x1.x0.x0.x0, qy1);
  mpn_set_mpz(Q.y.x1.x0.x1.x0, qy2);

  //e([a]P,[b]Q) を求める
  efp8_scm_fp(&aP,&P,a);
  efp8_scm(&bQ,&Q,b);

  #if 0
  efp8_println("P = ",&P);
  efp8_println("Q = ",&Q);

  efp8_scm_fp(&tmp1,&P,order_z);
  efp8_println("[r]P = ",&tmp1);
  efp8_scm(&tmp1,&Q,order_z);
  efp8_println("[r]Q = ",&tmp1);

  gmp_printf("a = %Zd\n",a);
  gmp_printf("b = %Zd\n",b);
  printf("---------------------------------\n");
  #endif
  efp_t mapped_P;
  efp_init(&mapped_P);
  efp2_t mapped_Q,mapped_Q_neg;
  efp2_init(&mapped_Q);
  efp2_init(&mapped_Q_neg);
  efp2_jacobian_t S;
  efp2_jacobian_init(&S);

  printf("---------------------------------\n");
  printf("check lazy montgomery pairing()\n");
  printf("---------------------------------\n");
  pre_miller_opt_ate_proj_loop_2NAF_lazy_montgomery(&f, &mapped_P, &mapped_Q, &mapped_Q_neg, &S, &aP, &bQ);
  miller_opt_ate_proj_loop_2NAF_lazy_montgomery(&f,&mapped_P,&mapped_Q,&mapped_Q_neg,&S);
  final_exp_lazy_montgomery(&e1,&f);
  //e(P,Q)^(a*b) を求める
  pre_miller_opt_ate_proj_loop_2NAF_lazy_montgomery(&f, &mapped_P, &mapped_Q, &mapped_Q_neg, &S, &P, &Q);
  miller_opt_ate_proj_loop_2NAF_lazy_montgomery(&f,&mapped_P,&mapped_Q,&mapped_Q_neg,&S);
  final_exp_lazy_montgomery(&e2,&f);

  fp8_pow_montgomery(&e2,&e2,ab);
  fp8_println_montgomery("e([a]P,[b]Q) = ",&e1);
  fp8_println_montgomery("e(P,Q)^(a*b) = ",&e2);
  if(fp8_cmp(&e1,&e2)==0)  {
  printf("=====================================================\n");
  printf("------------------bilinear!!-------------------------\n");
  printf("=====================================================\n\n\n");
  }else{
    printf("e([a]P,[b]Q) != e(P,Q)^(a*b)\n\n");
  }

  pre_miller_opt_ate_proj_loop_2NAF_lazy_montgomery(&f, &mapped_P, &mapped_Q, &mapped_Q_neg, &S, &aP, &bQ);
  CYBOZU_BENCH_C("miller_opt_ate_proj_2NAF_lazy_montgomery()", LOOP, miller_opt_ate_proj_loop_2NAF_lazy_montgomery,&f,&mapped_P,&mapped_Q,&mapped_Q_neg,&S);
  CYBOZU_BENCH_C("final_exp_lazy_montgomery()               ", LOOP, final_exp_lazy_montgomery, &e2, &e1);
  printf("---------------------------------\n");

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(ab);

  printf("*********************************************************************************************\n\n");
}

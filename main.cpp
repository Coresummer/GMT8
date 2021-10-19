
// #include "mcl.h"

#include <gmp.h>
#define TTT_INSTANCE_HERE

#include <cstdio>
// #include "libmcl/mcl.h"
#include "define.h"
#include "./time.h"
#include "scalar.h"
#include "mpn.h"
#include "fp.h"
#include "fp2.h"
#include "fp4.h"
#include "fp8.h"
// #include "test_fp.h"
#include "field_test.h"
#include "efp.h"
#include "efp2.h"
#include "test_efp.h"
#include "create.h"
#include "miller.h"
#include "final_exp.h"
#include "test_pairing.h"


int main(){
  // mcl_init();

  gmp_randinit_default(state);
  gmp_randseed_ui(state,(unsigned long int)time(NULL));
  tmp_init();
  create_prt();
  check_base();
  pre_montgomery();
  frobenius_precalculation();
  curve_search();
  create_weil();

  printf("*********************************************************************************************\n\n");
  

  //各関数の動作確認、コスト計算、時間計測など
  // test_fp_montgomery(CHECK_PAIRING_TIME_LOOP);
  // test_field(0, CHECK_PAIRING_TIME_LOOP, CHECK_PAIRING_TIME_LOOP, CHECK_PAIRING_TIME_LOOP);

  // test_fp(CHECK_PAIRING_TIME_LOOP);
  // test_fp2(CHECK_PAIRING_TIME_LOOP);
  // test_fp6(CHECK_PAIRING_TIME_LOOP);

  // check_fp_with_montgomery();
  // check_fp2_with_montgomery();
  // check_fp4_with_montgomery();
  // check_fp8_with_montgomery();

  // BENCH_fp2_fp6_mul_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_miller_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_finalexp_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_Pairingn_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // check_efp();
  // check_efp2();
  // check_efp4();
  // check_efp8();
  // check_g1_g2();
  // check_g2_dash();

  // check_pairing();
  // check_pairing_2NAF();
  // check_pairing_static();
  // check_pairing_count_2NAF_lazy_montgomery();
  // check_pairing_count_2NAF();
  // check_pairing_time_2NAF();

  printf("*********************************************************************************************\n\n");
  //playground
  fp8_t a,b;
  fp8_init(&a);
  fp8_init(&b);
  fp8_set_random(&a, state);
  fp8_set(&b,&a);
  // printf("here"),getchar();
  mpz_t exp;
  mpz_init_set(exp,fp8_total_r);
  mpz_pow_ui(exp,exp,4);
  // mpz_init(exp);
  // mpz_sub_ui(exp,prime_z,1);

  fp8_printf("A",&a);
  fp8_pow(&a, &a, exp);
  fp8_printf("A^4(p^8-1)/r",&a);
  final_exp(&b, &b);
  fp8_printf("B^4(p^8-1)/r",&b);  

  return 0;
}

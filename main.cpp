#include <gmp.h>
#define TTT_INSTANCE_HERE

#include <cstdio>
#include "mcl.h"
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
  const char *pStr = "0x347111bfc75e57d130de7be68437c8d75455d209459d421455023bee14df9fe75aa4734686ca3d08c1fa594100d79421d56c53899ee0f066fad9eb45b0985dbdbba2dcc1";
  mcl_init(pStr);

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

  // check_efp();
  // check_efp2();
  // check_efp4();
  // check_efp8();
  // check_g1_g2();
  // check_g2_dash();

  // check_pairing();
  // check_pairing_2NAF();
  // check_pairing_2NAF_costello();
  // check_pairing_static();
  // check_pairing_time_2NAF();
  // check_pairing_count_2NAF();
  // check_pairing_jac_count_2NAF_lazy_montgomery();
  // check_pairing_proj_count_2NAF_lazy_montgomery();
  // check_finalexp_pow_cost_count_2NAF();
  // check_finalexp_pow_cost_count_2NAF_montgomery();

  // BENCH_fp2_fp4_fp8_mul_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_miller_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_miller_coordinates(CHECK_PAIRING_TIME_LOOP);
  // BENCH_finalexp_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_Pairing_jac_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_Pairing_proj_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  BENCH_Pairing_proj_lazy_montgomery_static(CHECK_PAIRING_TIME_LOOP);
  printf("*********************************************************************************************\n\n");

  return 0;
}

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
  const char* sPtr = "0x347111bfc75e57d130de7be68437c8d75455d209459d421455023bee14df9fe75aa4734686ca3d08c1fa594100d79421d56c53899ee0f066fad9eb45b0985dbdbba2dcc1";
  mcl_init(sPtr);

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
  // check_pairing_static();
  // check_pairing_time_2NAF();
  // check_pairing_count_2NAF();
  // check_pairing_count_2NAF_lazy_montgomery();
  // check_finalexp_pow_cost_count_2NAF();
  // check_finalexp_pow_cost_count_2NAF_montgomery();

  // BENCH_fp2_fp4_fp8_mul_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_miller_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  // BENCH_finalexp_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);
  BENCH_Pairingn_lazy_montgomery(CHECK_PAIRING_TIME_LOOP);

  printf("*********************************************************************************************\n\n");
  // mpz_t expo;
  // mpz_init(expo);
  // mpz_mul(expo,prime_z,prime_z);
  // mpz_mul(expo,expo,expo);
  // fp8_t ANS, A;
  // fp8_init(&A);
  // fp8_init(&ANS);
  // fp8_set_random(&A, state);
  // fp8_println("A=", &A);
  // fp8_pow(&ANS,&A,expo);
  // fp8_println("ANS^p^4=", &ANS);
  // fp8_frobenius_map_p4(&ANS, &A);
  // // fp8_println("fp^4(ANS)=", &ANS);
  // fp2_t ANS1,ANS2,A;
  // fp2_init(&A);
  // fp2_init(&ANS1);
  // fp2_init(&ANS2);

  // fp_set_ui(&A.x0,3);
  // fp_set_ui(&A.x0,1);
  // fp2_mul_base_inv(&ANS1, &A);
  // fp2_println("new:",&ANS1);
  // fp2_mul_base_inv_classic(&ANS2, &A);
  // fp2_println("cls:",&ANS2);

  return 0;
}

#include "test_efp.h"
#include "efp2.h"
#include "efp8.h"
#include "fp.h"
#include "fp2.h"

void check_efp(){
  printf("check_efp() 開始\n");
  efp_t P,ANS;
  efp_init(&P);
  efp_init(&ANS);

  efp_rational_point(&P);
  efp_println("P = ",&P);
  printf("---------------------------------\n");

  printf("weil定理の確認\n");
  efp_scm(&ANS,&P,efp_total);
  efp_println("[p +1 -t]P = ",&ANS);
  printf("---------------------------------\n");

  efp_println("P = ",&P);//Aが変わっていないことの確認
  printf("*********************************************************************************************\n\n");
}

void check_efp2(){
  printf("check_efp2() 開始\n");
  efp2_t P,ANS;
  efp2_init(&P);
  efp2_init(&ANS);

  efp2_rational_point(&P);
  efp2_println("P = ",&P);
  printf("---------------------------------\n");

  printf("weil定理の確認\n");
  efp2_scm(&ANS,&P,efp2_total);
  efp2_println("[p^2 +1 -t2]P = ",&ANS);
  printf("---------------------------------\n");

  efp2_println("P = ",&P);//Aが変わっていないことの確認

  printf("*********************************************************************************************\n\n");
}

void check_efp4(){
  printf("check_efp4() 開始\n");
  efp4_t P,ANS;
  efp4_init(&P);
  efp4_init(&ANS);
  efp4_rational_point(&P);
  efp4_println("P = ",&P);
  // efp4_checkOnCurve(&P);

  printf("---------------------------------\n");

  printf("weil定理の確認\n");
  efp4_scm(&ANS,&P,efp4_total);
  // efp4_checkOnCurve(&ANS);
  efp4_println("[p^4 +1 -t4]P = ",&ANS);
  printf("---------------------------------\n");

  efp4_println("P = ",&P);//Aが変わっていないことの確認

  printf("*********************************************************************************************\n\n");
}

void check_efp8(){
  printf("check_efp8() 開始\n");
  efp8_t P,ANS;
  efp8_init(&P);
  efp8_init(&ANS);
  efp8_rational_point(&P);
  efp8_println("P = ",&P);
  // efp8_checkOnCurve(&P);

  printf("---------------------------------\n");

  printf("weil定理の確認\n");
  efp8_scm(&ANS,&P,efp8_total);
  // efp8_checkOnCurve(&ANS);
  efp8_println("[p^8 +1 -t8]P = ",&ANS);
  printf("---------------------------------\n");

  efp8_println("P = ",&P);//Aが変わっていないことの確認

  printf("*********************************************************************************************\n\n");
}

void check_g1_g2(){
  printf("check_g1_g2() 開始\n");
  efp8_t P,Q;
  efp_t _P;
  efp8_init(&P);
  efp8_init(&Q);
  efp_init(&_P);

  generate_g1(&P);
  fp_set(&_P.x, &P.x.x0.x0.x0);
  fp_set(&_P.y, &P.y.x0.x0.x0);
  efp_println("P in G1 = ",&_P);
  efp_scm(&_P,&_P,order_z);
  efp_println("[r]P = ",&_P);
  printf("---------------------------------\n");
  generate_g2(&Q);
  efp8_println("Q in G2 = ",&Q);
  efp8_scm(&Q,&Q,order_z);
  efp8_println("[r]Q = ",&Q);

  printf("*********************************************************************************************\n\n");
}

void check_g2_dash(){
  printf("check_g1_g2() 開始\n");
  efp8_t Q;
  efp2_t Q_dash;
  efp8_init(&Q);
  efp2_init(&Q_dash);

  generate_g2(&Q);

  efp8_println("Q in G2 = ",&Q);

  fp2_mul_base(&Q_dash.x,&Q.x.x0.x1);
  fp2_mul_base(&Q_dash.y,&Q.y.x1.x0);
  efp2_println("Q_dash in G2_dash = ",&Q_dash);
  efp2_checkOnTwsitCurve(&Q_dash);
  efp2_scm_dash(&Q_dash,&Q_dash,order_z);
  efp2_println("[r]Q_dash = ",&Q_dash);

  printf("*********************************************************************************************\n\n");
}
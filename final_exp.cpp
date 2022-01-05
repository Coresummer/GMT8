#include "final_exp.h"
#include "define.h"
#include "fp8.h"
#include <ELiPS/define.h>
#include <cstdio>

void final_exp_direct(fp8_t *ANS, fp8_t *A){
  fp8_pow(ANS, A, fp8_total_r);
}

void final_exp(fp8_t *ANS, fp8_t *A){
  static fp8_t tmp1_fp8, tmp2_fp8,tmp3_fp8,tmp4_fp8,tmp5_fp8,tmp6_fp8,tmp7_fp8;
  //easy
  fp8_inv(&tmp1_fp8, A);
  fp8_frobenius_map_p4(&tmp2_fp8, A);
  fp8_mul(&tmp1_fp8,&tmp1_fp8,&tmp2_fp8); //A^(p^4-1) = M

  //hard
  fp8_sqr_GS(ANS, &tmp1_fp8);
  fp8_sqr_GS(ANS, ANS);  //M^4

  fp8_frobenius_map_p1(&tmp2_fp8, &tmp1_fp8);  //M^p
  fp8_finalexpow_x_2NAF(&tmp3_fp8, &tmp1_fp8); //M^x
  fp8_mul(&tmp2_fp8,&tmp2_fp8,  &tmp3_fp8);    //M^(p+x) = M'

  fp8_frobenius_map_p2(&tmp3_fp8,&tmp2_fp8);   //M'^(p^2)
  fp8_finalexpow_x_2NAF(&tmp2_fp8,&tmp2_fp8);  //M'^(x^2)
  fp8_finalexpow_x_2NAF(&tmp2_fp8,&tmp2_fp8);  //M'^(x^2)
  fp8_mul(&tmp2_fp8,&tmp3_fp8,&tmp2_fp8);      //M" = M'(p^2 + x^2) = A(p^4-1)(p+x)(p^2+x^2)

  //L2 M"^4c = ((((4ℎ_𝑦+1)𝜒+〖4ℎ〗_𝑦 )𝜒−4ℎ_𝑦+1)𝜒)𝜒+4+〖4ℎ〗_𝑦^2
  fp8_sqr_GS(&tmp3_fp8, &tmp2_fp8);                  //M"^2
  fp8_sqr_GS(&tmp3_fp8, &tmp3_fp8);                  //M"^4
  fp8_finalexpow_hy_neg_2NAF(&tmp4_fp8, &tmp3_fp8);  //M"^-4hy
  fp8_frobenius_map_p4(&tmp3_fp8,&tmp3_fp8);         //M"^-4
  fp8_frobenius_map_p4(&tmp5_fp8,&tmp4_fp8);         //M"^4hy
  fp8_finalexpow_hy_neg_2NAF(&tmp6_fp8, &tmp4_fp8);  //M"^-4hy*-hy = M"^4hy^2 //not sure..

  fp8_mul(&tmp7_fp8,&tmp6_fp8,&tmp2_fp8);     //M"^(4hy^2+1)
  fp8_finalexpow_x_2NAF(&tmp7_fp8, &tmp7_fp8);//M"^(4hy^2+1)x
  fp8_mul(&tmp7_fp8,&tmp7_fp8,&tmp4_fp8);     //M"^(4hy^2+1)x -4hy
  fp8_finalexpow_x_2NAF(&tmp7_fp8, &tmp7_fp8);//M"^((4hy^2+1)x -4hy)x
  fp8_mul(&tmp7_fp8,&tmp7_fp8,&tmp5_fp8);     //M"^((4hy^2+1)x -4hy)x + 4hy
  fp8_mul(&tmp7_fp8,&tmp7_fp8,&tmp2_fp8);     //M"^((4hy^2+1)x -4hy)x + 4hy + 1
  fp8_finalexpow_x_2NAF(&tmp7_fp8, &tmp7_fp8);//M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x

  fp8_mul(&tmp7_fp8,&tmp7_fp8,&tmp3_fp8);     //M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4 
  fp8_finalexpow_x_2NAF(&tmp7_fp8, &tmp7_fp8);//M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4)x
  fp8_mul(&tmp7_fp8,&tmp7_fp8,&tmp6_fp8);     //M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4)x+ 4hy^2
  
  fp8_mul(ANS,&tmp7_fp8,ANS);           //(p^4-1)(4 + (p+x)(p^2+x^2)(((4hy^2+1)x -4hy + 1)x^2 + 4hy^2 + 4)
}


void final_exp_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp1_fp8, tmp2_fp8,tmp3_fp8,tmp4_fp8,tmp5_fp8,tmp6_fp8,tmp7_fp8;
  //easy
  fp8_inv_lazy_montgomery(&tmp1_fp8, A);
  fp8_frobenius_map_p4_montgomery(&tmp2_fp8, A);
  fp8_mul_lazy_montgomery(&tmp1_fp8,&tmp1_fp8,&tmp2_fp8); //A^(p^4-1) = M

  //hard
  fp8_sqr_GS_lazy_montgomery(ANS, &tmp1_fp8);
  fp8_sqr_GS_lazy_montgomery(ANS, ANS);  //M^4

  fp8_frobenius_map_p1_montgomery(&tmp2_fp8, &tmp1_fp8);  //M^p
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp3_fp8, &tmp1_fp8); //M^x
  fp8_mul_lazy_montgomery(&tmp2_fp8,&tmp2_fp8,  &tmp3_fp8);    //M^(p+x) = M'

  fp8_frobenius_map_p2_montgomery(&tmp3_fp8,&tmp2_fp8);   //M'^(p^2)
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp2_fp8,&tmp2_fp8);  //M'^(x^2)
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp2_fp8,&tmp2_fp8);  //M'^(x^2)
  fp8_mul_lazy_montgomery(&tmp2_fp8,&tmp3_fp8,&tmp2_fp8);      //M" = M'(p^2 + x^2) = A(p^4-1)(p+x)(p^2+x^2)

  //L2 M"^4c = ((((4ℎ_𝑦+1)𝜒+〖4ℎ〗_𝑦 )𝜒−4ℎ_𝑦+1)𝜒)𝜒+4+〖4ℎ〗_𝑦^2
  fp8_sqr_GS_lazy_montgomery(&tmp3_fp8, &tmp2_fp8);                  //M"^2
  fp8_sqr_GS_lazy_montgomery(&tmp3_fp8, &tmp3_fp8);                  //M"^4
  fp8_finalexpow_hy_neg_2NAF_lazy_montgomery(&tmp4_fp8, &tmp3_fp8);  //M"^-4hy
  fp8_frobenius_map_p4_montgomery(&tmp3_fp8,&tmp3_fp8);         //M"^-4
  fp8_frobenius_map_p4_montgomery(&tmp5_fp8,&tmp4_fp8);         //M"^4hy
  fp8_finalexpow_hy_neg_2NAF_lazy_montgomery(&tmp6_fp8, &tmp4_fp8);  //M"^-4hy*-hy = M"^4hy^2 //not sure..

  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp6_fp8,&tmp2_fp8);     //M"^(4hy^2+1)
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp7_fp8, &tmp7_fp8);//M"^(4hy^2+1)x
  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp7_fp8,&tmp4_fp8);     //M"^(4hy^2+1)x -4hy
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp7_fp8, &tmp7_fp8);//M"^((4hy^2+1)x -4hy)x
  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp7_fp8,&tmp5_fp8);     //M"^((4hy^2+1)x -4hy)x + 4hy
  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp7_fp8,&tmp2_fp8);     //M"^((4hy^2+1)x -4hy)x + 4hy + 1
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp7_fp8, &tmp7_fp8);//M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x

  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp7_fp8,&tmp3_fp8);     //M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4 
  fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp7_fp8, &tmp7_fp8);//M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4)x
  fp8_mul_lazy_montgomery(&tmp7_fp8,&tmp7_fp8,&tmp6_fp8);     //M"^(((4hy^2+1)x -4hy)x + 4hy + 1)x - 4)x + 4hy^2
  
  fp8_mul_lazy_montgomery(ANS,&tmp7_fp8,ANS);           //(p^4-1)(4 + (p+x)(p^2+x^2)(((4hy^2+1)x -4hy + 1)x^2 + 4hy^2 + 4)
}
#include "final_exp.h"
#include "define.h"
#include "fp8.h"
#include <ELiPS/define.h>
#include <cstdio>

void final_exp_direct(fp8_t *ANS, fp8_t *A){
  fp8_pow(ANS, A, fp8_total_r);
}

void final_exp(fp8_t *ANS, fp8_t *A){
  static fp8_t tmp1_fp8, tmp2_fp8,tmp3_fp8,tmp4_fp8,tmp5_fp8,tmp6_fp8;
  //easy
  fp8_inv(&tmp1_fp8, A);
  fp8_frobenius_map_p4(&tmp2_fp8, A);
  fp8_mul(&tmp1_fp8,&tmp1_fp8,&tmp2_fp8); //A^(p^4-1) = M
  //hard
  fp8_l1shift(&tmp2_fp8, &tmp1_fp8);
  fp8_l1shift(&tmp2_fp8, &tmp2_fp8);  //M^4

  fp8_frobenius_map_p1(&tmp3_fp8, &tmp1_fp8); //M^p
  fp8_pow(&tmp4_fp8, &tmp1_fp8,X_z);          //M^x
  fp8_mul(&tmp3_fp8,&tmp3_fp8,  &tmp4_fp8);   //M^(p+x) = M'

  fp8_frobenius_map_p2(&tmp4_fp8,&tmp3_fp8);  //M'^(p^2)
  fp8_pow(&tmp5_fp8,&tmp3_fp8,X_2);           //M'^(x^2)
  fp8_mul(&tmp3_fp8,&tmp4_fp8,&tmp5_fp8);      //M" = M'(p^2 + x^2) = A(p^4-1)(p+x)(p^2+x^2)

  //L2 M"^4c = ((((4ℎ_𝑦+1)𝜒+〖4ℎ〗_𝑦 )𝜒−4ℎ_𝑦+1)𝜒)𝜒+4+〖4ℎ〗_𝑦^2
  fp8_l1shift(&tmp1_fp8, &tmp3_fp8); //M"^2
  fp8_l1shift(&tmp1_fp8, &tmp1_fp8); //M"^4

  fp8_pow(&tmp4_fp8, &tmp3_fp8, fourhy_neg);  //M"^-4hy
  fp8_frobenius_map_p4(&tmp6_fp8, &tmp4_fp8); //M"^4hy
  fp8_mul(&tmp5_fp8,&tmp6_fp8,&tmp3_fp8);     //M"^(4hy+1)
  fp8_pow(&tmp5_fp8, &tmp5_fp8, X_z);         //M"^(4hy+1)^(x)
  fp8_mul(&tmp5_fp8,&tmp5_fp8,&tmp4_fp8);     //M"^(4hy+1)^(x) -4hy
  fp8_mul(&tmp5_fp8,&tmp5_fp8,&tmp3_fp8);     //M"^(4hy+1)^(x) -4hy + 1
  fp8_pow(&tmp5_fp8, &tmp5_fp8, X_2);         //M"^((4hy+1)^(x) -4hy + 1)x^2
  fp8_pow(&tmp6_fp8,&tmp4_fp8,hy_neg);        //M"^-4hy*-hy = M"^(4hy^2)
  fp8_mul(&tmp5_fp8,&tmp5_fp8,&tmp6_fp8);     //M"^((4hy+1)^(x) -4hy + 1)x^2 + (4hy^2)
  fp8_mul(&tmp5_fp8,&tmp5_fp8,&tmp1_fp8);     //M"^((4hy+1)^(x) -4hy + 1)x^2 + (4hy^2) + 4

  fp8_mul(ANS,&tmp5_fp8,&tmp2_fp8);           //M(4 + (p+x)(p^2+x^2)(((4hy)(x+1)x + 1 -4hy)x -2)x+1 +4hy^2))
}

// void final_exp(fp8_t *ANS,fp8_t *A){
//   static fp8_t tmp1_fp8, tmp2_fp8,tmp3_fp8,tmp4_fp8,tmp5_fp8,tmp6_fp8,tmp7_fp8,tmp8_fp8,tmp9_fp8,tmp10_fp8,tmp11_fp8,tmp12_fp8;

//   fp8_inv(&tmp1_fp8,A);
//   fp8_frobenius_map_p3(ANS,A);//(p^3)
//   fp8_mul(ANS,ANS,&tmp1_fp8);//(p^3-1)

//   fp8_frobenius_map_p1(&tmp1_fp8,ANS);//(p^3-1)(p)
//   fp8_mul(ANS,ANS,&tmp1_fp8);         //(p^3-1)(p+1) = M
//   //Second part
// //L1
//   fp8_frobenius_map_p3(&tmp2_fp8,ANS);                     //M^-1                                       //I
//   fp8_sqr_GS(&tmp3_fp8,&tmp2_fp8);               //M^-2                                       //S
  
//   //M^x                                        //128
//   fp8_finalexpow_x_2NAF(&tmp4_fp8, ANS);

//   fp8_mul(&tmp5_fp8,&tmp4_fp8,&tmp3_fp8);     //M^(x-2) = -t0/x                            //M
//   fp8_frobenius_map_p3(&tmp6_fp8,&tmp5_fp8);               //M^(-x+2) = t0/x                            //I
  
//   //M^(x^2-2x) = -t0                           //128
//   fp8_finalexpow_x_2NAF(&tmp7_fp8, &tmp5_fp8);

//   fp8_mul(&tmp8_fp8,&tmp6_fp8,&tmp7_fp8);     //M^(x^2-3x+2) = M^(-t0+t0/x) //reuse        //M
//   fp8_mul(&tmp8_fp8,&tmp8_fp8,ANS);           //M^(x^2-3x+3) = M^(-t0+t0/x+1)              //M
//   fp8_sqr_GS(&tmp9_fp8,&tmp8_fp8);
//   fp8_mul(&tmp8_fp8,&tmp9_fp8,&tmp8_fp8);

//   fp8_frobenius_map_p1(&tmp9_fp8,ANS);        //M^p                                        //p
//   fp8_frobenius_map_p3(&tmp12_fp8,&tmp7_fp8);              //M^t0                                       //I
//   fp8_mul(&tmp9_fp8,&tmp9_fp8,&tmp12_fp8);    //M^*(p+t0)                                  //M
//   fp8_mul(&tmp9_fp8,&tmp9_fp8,&tmp3_fp8);     //M^*(p+t0-2) = M' right L1                  //M

// //L2
//   //M'^3w                                       //80
//   fp8_finalexpow_3w_2NAF(&tmp5_fp8, &tmp9_fp8);

//   fp8_mul(&tmp6_fp8,&tmp5_fp8,&tmp9_fp8);    //M'(3w+1)                                    //M
  
//   //M'^9w^2 +3w                                 //80
//   fp8_finalexpow_3w_2NAF(&tmp7_fp8, &tmp6_fp8);

//   fp8_mul(&tmp11_fp8,&tmp7_fp8,&tmp9_fp8);   //M'^(9w^2+3w+1)                                //M

//   //M'^(9w^2+3w+1)(x-1)                           //128
//   fp8_finalexpow_x_1_2NAF(&tmp11_fp8, &tmp11_fp8);

//   fp8_mul(&tmp11_fp8,&tmp11_fp8,&tmp7_fp8);  //M'^(9w^2+3w+1)(x-1)+(9w^2+3w)                 //M
//   fp8_mul(&tmp11_fp8,&tmp11_fp8,&tmp5_fp8);  //M'^(9w^2+3w+1)(x-1)+(9w^2+6w)                 //M

//   //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1)            //128
//   fp8_finalexpow_x_1_2NAF(&tmp11_fp8, &tmp11_fp8);

//   fp8_mul(&tmp11_fp8,&tmp11_fp8,&tmp7_fp8);  //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1) +(9w^2+3w) //M

//   fp8_sqr_GS(&tmp12_fp8,&tmp6_fp8);             //M'2*(3w+1) = M'^6w+2                          //S
//   fp8_mul(&tmp12_fp8,&tmp12_fp8,&tmp9_fp8);  //M'^(6w+3)                                     //M
//   fp8_mul(&tmp11_fp8,&tmp11_fp8,&tmp12_fp8);  //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1) +(9w^2+3w) + (6w+3) //M

//   fp8_mul(ANS,&tmp8_fp8,&tmp11_fp8);         //left * right = M^(3(-t0-x+3) + (p+t0-2)((9w^2+3w+1)(x-1)+9w^2+6w)(x-1)+9w^2+9w + 3) //M
// }

// void final_exp_lazy_montgomery(fp8_t *ANS,fp8_t *A){
//   static fp8_t tmp1_fp8, tmp2_fp8,tmp3_fp8,tmp4_fp8,tmp5_fp8,tmp6_fp8,tmp7_fp8,tmp8_fp8,tmp9_fp8,tmp10_fp8,tmp11_fp8,tmp12_fp8;

//   fp8_inv_lazy_montgomery(&tmp1_fp8,A);
//   fp8_frobenius_map_p3_montgomery(ANS,A);//(p^3)
//   fp8_mul_lazy_montgomery(ANS,ANS,&tmp1_fp8);//(p^3-1)

//   fp8_frobenius_map_p1_montgomery(&tmp1_fp8,ANS);//(p^3-1)(p)
//   fp8_mul_lazy_montgomery(ANS,ANS,&tmp1_fp8);         //(p^3-1)(p+1) = M
//   fp8_sqr_lazy_montgomery(&tmp12_fp8, ANS);
//   //Second part
// //L1
//   fp8_frobenius_map_p3_montgomery(&tmp2_fp8,ANS);           //M^-1                                       //I
//   fp8_frobenius_map_p1_montgomery(&tmp3_fp8, ANS);          //M^p
//   fp8_finalexpow_x_2NAF_lazy_montgomery(&tmp4_fp8, ANS);    //M^x
//   fp8_mul_lazy_montgomery(&tmp9_fp8, &tmp2_fp8, &tmp3_fp8); // M^(p-1)
//   fp8_mul_lazy_montgomery(&tmp9_fp8, &tmp9_fp8, &tmp4_fp8); // M^(p-1+x) = M'

// //L2
//   //M'^3w                                       //80
//   fp8_finalexpow_3w_2NAF_lazy_montgomery(&tmp5_fp8, &tmp9_fp8);

//   fp8_mul_lazy_montgomery(&tmp6_fp8,&tmp5_fp8,&tmp9_fp8);    //M'(3w+1)                                    //M
  
//   //M'^9w^2 +3w                                 //80
//   fp8_finalexpow_3w_2NAF_lazy_montgomery(&tmp7_fp8, &tmp6_fp8);

//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp7_fp8,&tmp9_fp8);   //M'^(9w^2+3w+1)                                //M

//   //M'^(9w^2+3w+1)(x-1)                           //128
//   fp8_finalexpow_x_1_2NAF_lazy_montgomery(&tmp11_fp8, &tmp11_fp8);

//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp11_fp8,&tmp7_fp8);  //M'^(9w^2+3w+1)(x-1)+(9w^2+3w)                 //M
//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp11_fp8,&tmp5_fp8);  //M'^(9w^2+3w+1)(x-1)+(9w^2+6w)                 //M

//   //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1)            //128
//   fp8_finalexpow_x_1_2NAF_lazy_montgomery(&tmp11_fp8, &tmp11_fp8);

//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp11_fp8,&tmp7_fp8);  //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1) +(9w^2+3w) //M
//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp11_fp8,&tmp5_fp8);  //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1) +(9w^2+6w) //M
//   fp8_mul_lazy_montgomery(&tmp11_fp8,&tmp11_fp8,&tmp5_fp8);  //M'^((9w^2+3w+1)(x-1)+9w^2+6w)(x-1) +(9w^2+9w) //M

//   fp8_mul_lazy_montgomery(ANS,ANS,&tmp11_fp8);         //left * right = M^(1 + (p-1+x)((9w^2+3w+1)(x-1)+9w^2+6w)(x-1)+9w^2+9w ) //M
//   fp8_mul_lazy_montgomery(ANS,ANS,&tmp12_fp8);         //left * right = M^(3 + (p-1+x)((9w^2+3w+1)(x-1)+9w^2+6w)(x-1)+9w^2+9w ) //M
// }
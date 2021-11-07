

#include "fp8.h"
#include "define.h"
#include "fp.h"
#include "fp2.h"
#include "fp4.h"

void fp8_init(fp8_t *A){
  fp4_init(&A->x0);
  fp4_init(&A->x1);
}

void fp8_printf(std::string str ,fp8_t *A){
  gmp_printf("%s(",str.c_str());
  fp4_printf("",&A->x0);
  gmp_printf(",");
  fp4_printf("",&A->x1);
  gmp_printf(")\n");
}

void fp8_println(std::string str ,fp8_t *A){
  gmp_printf("%s(",str.c_str());
  fp4_printf("",&A->x0);
  gmp_printf(",");
  fp4_printf("",&A->x1);
  gmp_printf(")\n");
}

void fpd8_printf(std::string str ,fpd8_t *A){
  gmp_printf("%s(",str.c_str());
  fpd4_printf("",&A->x0);
  gmp_printf(",");
  fpd4_printf("",&A->x1);
  gmp_printf(")\n");
}

void fpd8_println(std::string str ,fpd8_t *A){
  gmp_printf("%s(",str.c_str());
  fpd4_printf("",&A->x0);
  gmp_printf(",");
  fpd4_printf("",&A->x1);
  gmp_printf(")\n");
}

void fp8_printf_montgomery(std::string str, fp8_t *A) {
  gmp_printf("%s(", str.c_str());
  fp4_printf_montgomery("", &A->x0);
  gmp_printf(",");
  fp4_printf_montgomery("", &A->x1);
  gmp_printf(")\n");
}

void fp8_println_montgomery(std::string str, fp8_t *A) {
  gmp_printf("%s(", str.c_str());
  fp4_printf_montgomery("", &A->x0);
  gmp_printf(",");
  fp4_printf_montgomery("", &A->x1);
  gmp_printf(")\n");
}

void fp8_set(fp8_t *ANS,fp8_t *A){
  fp4_set(&ANS->x0,&A->x0);
  fp4_set(&ANS->x1,&A->x1);
}

void fpd8_set(fpd8_t *ANS, fpd8_t *A) {
  fpd4_set(&ANS->x0, &A->x0);
  fpd4_set(&ANS->x1, &A->x1);
}

void fp8_set_ui(fp8_t *ANS,unsigned long int UI){
  fp4_set_ui(&ANS->x0,UI);
  fp4_set_ui(&ANS->x1,0);
}

void fp8_set_ui_ui(fp8_t *ANS,unsigned long int UI){
  fp4_set_ui(&ANS->x0,UI);
  fp4_set_ui(&ANS->x1,UI);
}

void fp8_set_mpn(fp8_t *ANS,mp_limb_t *A){
  fp4_set_mpn(&ANS->x0,A);
  fp4_set_ui(&ANS->x1,0);
}

void fp8_set_neg(fp8_t *ANS,fp8_t *A){
  fp4_set_neg(&ANS->x0,&A->x0);
  fp4_set_neg(&ANS->x1,&A->x1);
}

void fp8_set_neg_montgomery(fp8_t *ANS,fp8_t *A){
  fp4_set_neg_montgomery(&ANS->x0,&A->x0);
  fp4_set_neg_montgomery(&ANS->x1,&A->x1);
}

void fp8_set_conj(fp8_t *ANS,fp8_t *A){
  fp4_set(&ANS->x0,&A->x0);
  fp4_set_neg(&ANS->x1,&A->x1);
}

void fp8_set_conj_montgomery(fp8_t *ANS,fp8_t *A){
  fp4_set(&ANS->x0,&A->x0);
  fp4_set_neg_montgomery(&ANS->x1,&A->x1);
}

void fp8_to_montgomery(fp8_t *ANS,fp8_t *A){
  fp4_to_montgomery(&ANS->x0,&A->x0);
  fp4_to_montgomery(&ANS->x1,&A->x1);
}

void fp8_mod_montgomery(fp8_t *ANS,fp8_t *A){
  mpn_mod_montgomery(ANS->x0.x0.x0.x0,FPLIMB,A->x0.x0.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x0.x0.x1.x0,FPLIMB,A->x0.x0.x1.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x0.x1.x0.x0,FPLIMB,A->x0.x1.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x0.x1.x1.x0,FPLIMB,A->x0.x1.x1.x0,FPLIMB);

  mpn_mod_montgomery(ANS->x1.x0.x0.x0,FPLIMB,A->x1.x0.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x1.x0.x1.x0,FPLIMB,A->x1.x0.x1.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x1.x1.x0.x0,FPLIMB,A->x1.x1.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x1.x1.x1.x0,FPLIMB,A->x1.x1.x1.x0,FPLIMB);
}

void fp8_mod_montgomery_double(fp8_t *ANS,fpd8_t *A){
  mpn_mod_montgomery(ANS->x0.x0.x0.x0,FPLIMB,A->x0.x0.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x0.x0.x1.x0,FPLIMB,A->x0.x0.x1.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x0.x1.x0.x0,FPLIMB,A->x0.x1.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x0.x1.x1.x0,FPLIMB,A->x0.x1.x1.x0,FPLIMB2);

  mpn_mod_montgomery(ANS->x1.x0.x0.x0,FPLIMB,A->x1.x0.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x1.x0.x1.x0,FPLIMB,A->x1.x0.x1.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x1.x1.x0.x0,FPLIMB,A->x1.x1.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x1.x1.x1.x0,FPLIMB,A->x1.x1.x1.x0,FPLIMB2);
}

void fp8_lshift(fp8_t *ANS, fp8_t *A, unsigned long int UI) {
  fp4_lshift(&ANS->x0, &A->x0, UI);
  fp4_lshift(&ANS->x1, &A->x1, UI);
}
void fp8_l1shift(fp8_t *ANS, fp8_t *A) {
  fp4_l1shift(&ANS->x0, &A->x0);
  fp4_l1shift(&ANS->x1, &A->x1);
}

void fp8_l1shift_nonmod_single(fp8_t *ANS, fp8_t *A) {
  fp4_l1shift_nonmod_single(&ANS->x0, &A->x0);
  fp4_l1shift_nonmod_single(&ANS->x1, &A->x1);
}

void fp8_l1shift_nonmod_double(fpd8_t *ANS, fpd8_t *A) {
  fp4_l1shift_nonmod_double(&ANS->x0, &A->x0);
  fp4_l1shift_nonmod_double(&ANS->x1, &A->x1);
}


void fp8_r1shift(fp8_t *ANS, fp8_t *A) {
  fp4_r1shift(&ANS->x0, &A->x0);
  fp4_r1shift(&ANS->x1, &A->x1);
}

void fp8_set_random(fp8_t *ANS,gmp_randstate_t state){
  fp4_set_random(&ANS->x0,state);
  fp4_set_random(&ANS->x1,state);
}

void fp8_mul(fp8_t *ANS,fp8_t *A,fp8_t *B){ 
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);
  fp8_set(&tmp_B,B);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  fp4_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //ab+bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
} 

// void fp8_mul_lazy(fp8_t *ANS,fp8_t *A,fp8_t *B){ 
//   static fp8_t tmp_A,tmp_B;
//   fp8_set(&tmp_A,A);
//   fp8_set(&tmp_B,B);

//   static fp4_t tmp3_fp,tmp4_fp;
//   static fp4_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
//   fp4_mul_lazy(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
//   fp4_mul_lazy(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
//   fp4_l1shift(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
//   fp4_add_nonmod_single(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2
//   fp4_mod(&ANS->x0,tmp4_fpd.x0,FPLIMB);

//   fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
//   fp4_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
//   fp4_mul_lazy(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
//   fp4_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
//   fp4_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
//   fp4_mod(&ANS->x1,tmp4_fpd.x0,FPLIMB2);
// } 

void fp8_mul_lazy_montgomery(fp8_t *ANS,fp8_t *A,fp8_t *B){ 
  // static fp8_t tmp_A,tmp_B;
  // fp8_set(&tmp_A,A);
  // fp8_set(&tmp_B,B);

  // static fp4_t tmp3_fp,tmp4_fp;
  // static fpd8_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
  // fp4_mul_nonmod(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
  // fp4_mul_nonmod(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
  // fp4_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
  // fp4_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2
  // mpn_mod_montgomery(ANS->x0.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  // fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  // fp4_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  // fp4_mul_nonmod(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  // fp4_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
  // fp4_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
  // mpn_mod_montgomery(ANS->x1.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);
  fp8_set(&tmp_B,B);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul_lazy_montgomery(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  fp4_mul_lazy_montgomery(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //ab+bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul_lazy_montgomery(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
} 


void fp8_mul_nonmod_montgomery(fpd8_t *ANS, fp8_t *A, fp8_t *B) {
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);
  fp8_set(&tmp_B,B);

  static fp4_t tmp3_fp,tmp4_fp;
  static fpd4_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
  fp4_mul_nonmod_montgomery(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
  fp4_mul_nonmod_montgomery(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
  fp4_l1shift_nonmod_double(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
  fp4_add_nonmod_double(&ANS->x0, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2

  fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul_nonmod_montgomery(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  fp4_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
  fp4_sub_nonmod_double(&ANS->x1,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
}

void fp8_mul_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI){
  fp4_mul_ui(&ANS->x0,&A->x0,UI);
  fp4_mul_ui(&ANS->x1,&A->x1,UI);
}
void fp8_mul_ui_nonmod_single(fp8_t *ANS, fp8_t *A, unsigned long int UI) {
  fp4_mul_ui_nonmod_single(&ANS->x0, &A->x0, UI);
  fp4_mul_ui_nonmod_single(&ANS->x1, &A->x1, UI);
}
void fp8_mul_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B){
  fp4_mul_mpn(&ANS->x0,&A->x0,B);
  fp4_mul_mpn(&ANS->x1,&A->x1,B);
}

void fp8_mul_mpn_montgomery(fp8_t *ANS,fp8_t *A,mp_limb_t *B){
  mpn_mulmod_montgomery(ANS->x0.x0.x0.x0,FPLIMB,A->x0.x0.x0.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x0.x0.x1.x0,FPLIMB,A->x0.x0.x1.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x0.x1.x0.x0,FPLIMB,A->x0.x1.x0.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x0.x1.x1.x0,FPLIMB,A->x0.x1.x1.x0,FPLIMB,B,FPLIMB);

  mpn_mulmod_montgomery(ANS->x1.x0.x0.x0,FPLIMB,A->x1.x0.x0.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x1.x0.x1.x0,FPLIMB,A->x1.x0.x1.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x1.x1.x0.x0,FPLIMB,A->x1.x1.x0.x0,FPLIMB,B,FPLIMB);
  mpn_mulmod_montgomery(ANS->x1.x1.x1.x0,FPLIMB,A->x1.x1.x1.x0,FPLIMB,B,FPLIMB);
}

void fp8_sqr(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp4_sqr(&tmp1_fp, &tmp_A.x0);  //a^2
  fp4_sqr(&tmp2_fp, &tmp_A.x1);  //b^2
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //b^2@^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp); //a^2+b^2@^

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp4_sqr(&tmp3_fp, &tmp3_fp);  //(a+b)^2
  fp4_sub(&tmp3_fp,&tmp3_fp, &tmp1_fp); //(a+b)^2 - a^2 
  fp4_sub(&ANS->x1,&tmp3_fp, &tmp2_fp); //(a+b)^2 - a^2 - b^2 
}

// void fp8_sqr_lazy(fp8_t *ANS,fp8_t *A){
//   static fp8_t tmp_A;
//   fp8_set(&tmp_A,A);
//   static fp4_t tmp3_fp;
//   static fpd8_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
//   fp4_sqr_nonmod_montgomery(&tmp1_fpd, &tmp_A.x0);  //a^2
//   fp4_sqr_nonmod_montgomery(&tmp2_fpd, &tmp_A.x1);  //b^2

//   fp4_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
//   fp4_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^
//   fp4_mod(&ANS->x0,tmp4_fpd.x0,FPLIMB2);

//   fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
//   fp4_sqr_nonmod(&tmp3_fpd, &tmp3_fp);  //(a+b)^2

//   fp4_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
//   fp4_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 
//   fp4_mod(&ANS->x1,tmp4_fpd.x0,FPLIMB2);
// }

void fp8_sqr_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  // static fp8_t tmp_A;
  // fp8_set(&tmp_A,A);
  // static fp4_t tmp3_fp;
  // static fpd8_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
  // fp4_sqr_nonmod(&tmp1_fpd, &tmp_A.x0);  //a^2
  // fp4_sqr_nonmod(&tmp2_fpd, &tmp_A.x1);  //b^2
  // fp4_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
  // fp4_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^
  // mpn_mod_montgomery(ANS->x0.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  // fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  // fp4_sqr_nonmod(&tmp3_fpd, &tmp3_fp);  //(a+b)^2
  // fp4_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
  // fp4_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 
  // mpn_mod_montgomery(ANS->x1.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp4_sqr_lazy_montgomery(&tmp1_fp, &tmp_A.x0);  //a^2
  fp4_sqr_lazy_montgomery(&tmp2_fp, &tmp_A.x1);  //b^2
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //b^2@^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp); //a^2+b^2@^

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp4_sqr_lazy_montgomery(&tmp3_fp, &tmp3_fp);  //(a+b)^2
  fp4_sub(&tmp3_fp,&tmp3_fp, &tmp1_fp); //(a+b)^2 - a^2 
  fp4_sub(&ANS->x1,&tmp3_fp, &tmp2_fp); //(a+b)^2 - a^2 - b^2 
}

void fp8_sqr_nonmod_montgomery(fpd8_t *ANS, fp8_t *A) {

  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);
  static fp4_t tmp3_fp;
  static fpd4_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
  fp4_sqr_nonmod_montgomery(&tmp1_fpd, &tmp_A.x0);  //a^2
  fp4_sqr_nonmod_montgomery(&tmp2_fpd, &tmp_A.x1);  //b^2
  fp4_l1shift_nonmod_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
  fp4_add_nonmod_double(&ANS->x0, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^

  fp4_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp4_sqr_nonmod_montgomery(&tmp3_fpd, &tmp3_fp);  //(a+b)^2
  fp4_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
  fp4_sub_nonmod_double(&ANS->x1,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 

}

void fp8_sqr_GS(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp4_sqr(&tmp1_fp, &tmp_A.x0);  //a^2
  fp4_add(&tmp2_fp, &tmp_A.x0, &tmp_A.x1);  //(a+b)
  fp4_sqr(&tmp2_fp,&tmp2_fp);  //(a+b)^2

  fp4_sub_ui(&tmp3_fp,&tmp1_fp,1);//a^2-1
  fp4_add(&ANS->x0,&tmp3_fp,&tmp1_fp);//2a^2-1

  fp4_mul_base_inv(&tmp3_fp,&tmp3_fp);//(a^2-1)/2
  
  fp4_sub(&ANS->x1,&tmp2_fp,&tmp1_fp);//(a+b)^2 - a^2
  fp4_sub(&ANS->x1,&ANS->x1, &tmp3_fp);//(a+b)^2 - a^2 - (a^2-1)/2
}

void fp8_sqr_GS_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp4_sqr_lazy_montgomery(&tmp1_fp, &tmp_A.x0);  //a^2
  fp4_add(&tmp2_fp, &tmp_A.x0, &tmp_A.x1);  //(a+b)
  fp4_sqr_lazy_montgomery(&tmp2_fp,&tmp2_fp);  //(a+b)^2

  fp2_sub_mpn(&tmp3_fp.x0, &tmp1_fp.x0, oneMR.x0);//a^2-1
  fp2_set(&tmp3_fp.x1, &tmp1_fp.x1);
  fp4_add(&ANS->x0,&tmp3_fp,&tmp1_fp);//2a^2-1

  fp4_mul_base_inv_montgomery(&tmp3_fp,&tmp3_fp);//(a^2-1)/z

  fp4_sub(&ANS->x1,&tmp2_fp,&tmp1_fp);//(a+b)^2 - a^2
  fp4_sub(&ANS->x1,&ANS->x1, &tmp3_fp);//(a+b)^2 - a^2 - (a^2-1)/z
}

void fp8_add(fp8_t *ANS,fp8_t *A,fp8_t *B){
  fp4_add(&ANS->x0,&A->x0,&B->x0);
  fp4_add(&ANS->x1,&A->x1,&B->x1);

}

void fp8_add_nonmod_single(fp8_t *ANS,fp8_t *A,fp8_t *B){
  fp4_add_nonmod_single(&ANS->x0,&A->x0,&B->x0);
  fp4_add_nonmod_single(&ANS->x1,&A->x1,&B->x1);

}

void fp8_add_nonmod_double(fpd8_t *ANS,fpd8_t *A,fpd8_t *B){
  fp4_add_nonmod_double(&ANS->x0,&A->x0,&B->x0);
  fp4_add_nonmod_double(&ANS->x1,&A->x1,&B->x1);

}

void fp8_add_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI){
  fp4_add_ui(&ANS->x0,&A->x0,UI);
  fp4_set(&ANS->x1,&A->x1);

}

void fp8_add_ui_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI){
  fp4_add_ui(&ANS->x0,&A->x0,UI);
  fp4_add_ui(&ANS->x1,&A->x1,UI);

}

void fp8_add_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B){
  fp4_add_mpn(&ANS->x0,&A->x0,B);
  fp4_add_mpn(&ANS->x1,&A->x1,B);

}

void fp8_sub(fp8_t *ANS,fp8_t *A,fp8_t *B){
  fp4_sub(&ANS->x0,&A->x0,&B->x0);
  fp4_sub(&ANS->x1,&A->x1,&B->x1);

}

void fp8_sub_nonmod_single(fp8_t *ANS,fp8_t *A,fp8_t *B){
  fp4_sub_nonmod_single(&ANS->x0,&A->x0,&B->x0);
  fp4_sub_nonmod_single(&ANS->x1,&A->x1,&B->x1);

}

void fp8_sub_nonmod_double(fpd8_t *ANS,fpd8_t *A,fpd8_t *B){
  fp4_sub_nonmod_double(&ANS->x0,&A->x0,&B->x0);
  fp4_sub_nonmod_double(&ANS->x1,&A->x1,&B->x1);

;}

void fp8_sub_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI){
  fp4_sub_ui(&ANS->x0,&A->x0,UI);
  fp4_set(&ANS->x1,&A->x1);

}

void fp8_sub_ui_ui(fp8_t *ANS,fp8_t *A,unsigned long int UI){
  fp4_sub_ui(&ANS->x0,&A->x0,UI);
  fp4_sub_ui(&ANS->x1,&A->x1,UI);

}

void fp8_sub_mpn(fp8_t *ANS,fp8_t *A,mp_limb_t *B){
  fp4_sub_mpn(&ANS->x0,&A->x0,B);
  fp4_sub_mpn(&ANS->x1,&A->x1,B);

}

void fp8_inv(fp8_t *ANS,fp8_t *A){  
	static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
  fp4_set(&tmp1_fp,&A->x0);  //a
  fp4_set_neg(&tmp2_fp,&A->x1); //b

  fp4_sqr(&tmp3_fp,&tmp1_fp);  //a^2
  fp4_mul(&tmp4_fp,&tmp2_fp,&A->x1);  //b^2*c
  fp4_mul_base(&tmp4_fp,&tmp4_fp);
  fp4_add(&tmp3_fp,&tmp3_fp,&tmp4_fp); //a^2 - b^2c //mabe need mul_basis?

  fp4_inv(&tmp3_fp,&tmp3_fp);  //  (a^2 - b^2c)^-1
  fp4_mul(&ANS->x0,&tmp1_fp,&tmp3_fp); // a*(a^2- b^2c)^-1
  fp4_mul(&ANS->x1,&tmp2_fp,&tmp3_fp); //-b*(a^2 - b^2c)^-1
}

void fp8_inv_lazy(fp8_t *ANS, fp8_t *A) {
	static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
  fp4_set(&tmp1_fp,&A->x0);  //a
  fp4_set_neg(&tmp2_fp,&A->x1); //b

  fp4_sqr(&tmp3_fp,&tmp1_fp);  //a^2
  fp4_mul(&tmp4_fp,&tmp2_fp,&A->x1);  //b^2*c
  fp4_mul_base(&tmp4_fp,&tmp4_fp);
  fp4_add(&tmp3_fp,&tmp3_fp,&tmp4_fp); //a^2 - b^2c //mabe need mul_basis?

  fp4_inv(&tmp3_fp,&tmp3_fp);  //  (a^2 - b^2c)^-1
  fp4_mul(&ANS->x0,&tmp1_fp,&tmp3_fp); // a*(a^2- b^2c)^-1
  fp4_mul(&ANS->x1,&tmp2_fp,&tmp3_fp); //-b*(a^2 - b^2c)^-1
}

void fp8_inv_lazy_montgomery(fp8_t *ANS, fp8_t *A) {
  static fp4_t tmp1_fp, tmp2_fp;
  static fp4_t tmp3;
  static fp4_t tmp1, tmp2;
  fp4_set(&tmp1_fp, &A->x0);
  fp4_set_neg(&tmp2_fp, &A->x1);

  fp4_sqr_lazy_montgomery(&tmp1, &tmp1_fp); //a^2
  fp4_mul_lazy_montgomery(&tmp2, &tmp2_fp, &A->x1);//b^2
  fp4_mul_base(&tmp2,&tmp2);//b^2*c
  fp4_add(&tmp3, &tmp1, &tmp2);//a^2-b^2*c

  fp4_inv_lazy_montgomery(&tmp3, &tmp3);
  fp4_mul_lazy_montgomery(&ANS->x0, &tmp1_fp, &tmp3);
  fp4_mul_lazy_montgomery(&ANS->x1, &tmp2_fp, &tmp3);
}

void fp8_pow(fp8_t *ANS,fp8_t *A,mpz_t scalar){
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);
  fp8_t tmp;
  fp8_init(&tmp);
  fp8_set(&tmp,A);

  for(i=1;i<length; i++){
    fp8_sqr(&tmp,&tmp);
    if(binary[i]=='1')  fp8_mul(&tmp,A,&tmp);
  }
  fp8_set(ANS,&tmp);
}

void fp8_pow_GS(fp8_t *ANS,fp8_t *A,mpz_t scalar){
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);
  fp8_t tmp;
  fp8_init(&tmp);
  fp8_set(&tmp,A);

  for(i=1;i<length; i++){
    fp8_sqr_GS(&tmp,&tmp);
    if(binary[i]=='1')  fp8_mul(&tmp,A,&tmp);
  }
  fp8_set(ANS,&tmp);
}


void fp8_finalexpow_x_2NAF(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_x.size()-2);i!=-1;i--){
    switch(finalexp_pow_x[i]){
      case 0:
        fp8_sqr_GS(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS, ANS,&A_inv);

        break;
      default:
        break;
    }
  }
}

void fp8_finalexpow_hy_neg_2NAF(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_hy.size()-2);i!=-1;i--){
    switch(finalexp_pow_hy[i]){
      case 0:
        fp8_sqr_GS(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS, ANS,&A_inv);
        break;
      default:
        break;
    }
  }
}

void fp8_finalexpow_4hy_neg_2NAF(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_4hy.size()-2);i!=-1;i--){
    switch(finalexp_pow_4hy[i]){
      case 0:
        fp8_sqr_GS(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS(ANS, ANS);
        fp8_mul(ANS, ANS,&A_inv);
        break;
      default:
        break;
    }
  }
}


void fp8_finalexpow_x_2NAF_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4_montgomery(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_x.size()-2);i!=-1;i--){
    switch(finalexp_pow_x[i]){
      case 0:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS, ANS,&A_inv);

        break;
      default:
        break;
    }
  }
}

void fp8_finalexpow_hy_neg_2NAF_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4_montgomery(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_hy.size()-2);i!=-1;i--){
    switch(finalexp_pow_hy[i]){
      case 0:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS, ANS,&A_inv);
        break;
      default:
        break;
    }
  }
}

void fp8_finalexpow_4hy_neg_2NAF_lazy_montgomery(fp8_t *ANS,fp8_t *A){
  static fp8_t tmp_A;
  fp8_set(&tmp_A,A);

  static fp8_t A_inv;
  fp8_frobenius_map_p4_montgomery(&A_inv, A);

  fp8_set(ANS,&tmp_A);
  for(int i=(finalexp_pow_4hy.size()-2);i!=-1;i--){
    switch(finalexp_pow_4hy[i]){
      case 0:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        break;
      case 1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS,ANS,&tmp_A);
        break;
      case -1:
        fp8_sqr_GS_lazy_montgomery(ANS, ANS);
        fp8_mul_lazy_montgomery(ANS, ANS,&A_inv);
        break;
      default:
        break;
    }
  }
}

int fp8_legendre(fp8_t *A){
  fp8_t tmp;
  mpz_t expo;
  fp8_init(&tmp);
  mpz_init(expo);

  //(p^2 -1)/2 を計算
  mpz_pow_ui(expo,prime_z,8);
  mpz_sub_ui(expo,expo,1);
  mpz_tdiv_q_ui(expo,expo,2);
  fp8_pow(&tmp,A,expo);

  if(fp8_cmp_one(&tmp)==0){
    mpz_clear(expo);
    return 1;
  }
  else{
     mpz_clear(expo);
    return -1;
  }
}

void fp8_sqrt(fp8_t *ANS,fp8_t *A){
  fp8_t x,y,t,k,n,tmp;
  fp8_init(&x);
  fp8_init(&y);
  fp8_init(&t);
  fp8_init(&k);
  fp8_init(&n);
  fp8_init(&tmp);
  unsigned long int e,m;
  mpz_t exp,q,z,result;
  mpz_init(exp);
  mpz_init(q);
  mpz_init(z);
  mpz_init(result);
  //gmp_randstate_t state;
  //gmp_randinit_default(state);
  //gmp_randseed_ui(state,(unsigned long)time(NULL));

  fp8_set_random(&n,state);
  while(fp8_legendre(&n)!=-1){
    fp8_set_random(&n,state);
  }
  mpz_pow_ui(q,prime_z,8);
  mpz_sub_ui(q,q,1);
  mpz_mod_ui(result,q,2);
  e=0;
  while(mpz_cmp_ui(result,0)==0){
    mpz_tdiv_q_ui(q,q,2);
    mpz_mod_ui(result,q,2);
    e++;
  }
  fp8_pow(&y,&n,q);
  mpz_set_ui(z,e);
  mpz_sub_ui(exp,q,1);
  mpz_tdiv_q_ui(exp,exp,2);
  fp8_pow(&x,A,exp);
  fp8_mul(&tmp,&x,&x);
  fp8_mul(&k,&tmp,A);
  fp8_mul(&x,&x,A);
  while(fp8_cmp_one(&k)!=0){
    m=1;
    mpz_ui_pow_ui(exp,2,m);
    fp8_pow(&tmp,&k,exp);
    while(fp8_cmp_one(&tmp)!=0){
      m++;
      mpz_ui_pow_ui(exp,2,m);
      fp8_pow(&tmp,&k,exp);
    }
    mpz_sub_ui(exp,z,m);
    mpz_sub_ui(exp,exp,1);
    mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
    fp8_pow(&t,&y,result);
    fp8_mul(&y,&t,&t);
    mpz_set_ui(z,m);
    fp8_mul(&x,&x,&t);
    fp8_mul(&k,&k,&y);
  }
  fp8_set(ANS,&x);

  mpz_clear(exp);
  mpz_clear(q);
  mpz_clear(z);
  mpz_clear(result);
}


void fp8_pow_montgomery(fp8_t *ANS, fp8_t *A, mpz_t scalar) {
  int length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);
  fp8_t tmp;
  fp8_init(&tmp); // not need?
  fp8_set(&tmp, A);

  for (int i = 1; i < length; i++) {
    fp8_sqr_lazy_montgomery(&tmp, &tmp);
    if (binary[i] == '1') {
      fp8_mul_lazy_montgomery(&tmp, A, &tmp);
    }
  }
  fp8_set(ANS, &tmp);
}

int fp8_cmp(fp8_t *A,fp8_t *B){
  if(fp4_cmp(&A->x0,&B->x0)==0 && fp4_cmp(&A->x1,&B->x1)==0){
    return 0;
  }
  return 1;
}

int fp8_cmp_ui(fp8_t *A,unsigned long int UI){
  if(fp4_cmp_ui(&A->x0,UI)==0 && fp4_cmp_ui(&A->x1,UI)==0){
    return 0;
  }
  return 1;
}

int fp8_cmp_mpn(fp8_t *A,mp_limb_t *B){
  if(fp4_cmp_mpn(&A->x0,B)==0 && fp4_cmp_mpn(&A->x1,B)==0){
    return 0;
  }
  return 1;
}

int fp8_cmp_zero(fp8_t *A){
  if(fp4_cmp_zero(&A->x0)==0 && fp4_cmp_zero(&A->x1)==0 ){
    return 0;
  }
  return 1;
}

int fp8_cmp_one(fp8_t *A){
  if(fp4_cmp_one(&A->x0)==0 && fp4_cmp_zero(&A->x1)==0 ){
    return 0;
  }
  return 1;
}

void fp8_frobenius_map_p1(fp8_t *ANS, fp8_t *A){
  fp2_swap(&ANS->x0.x0,&A->x0.x0);
  
  fp2_swap(&ANS->x0.x1,&A->x0.x1);
  fp2_mul(&ANS->x0.x1,&ANS->x0.x1,&frobenius_2_4);

  fp2_swap(&ANS->x1.x0,&A->x1.x0);
  fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_1_4);

  fp2_swap(&ANS->x1.x1,&A->x1.x1);
  fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&frobenius_3_4);
}

void fp8_frobenius_map_p2(fp8_t *ANS, fp8_t *A){
  fp2_set(&ANS->x0.x0,&A->x0.x0);
  
  fp2_set_neg(&ANS->x0.x1,&A->x0.x1);

  fp2_mul(&ANS->x1.x0,&A->x1.x0,&frobenius_2_4);

  fp2_mul(&ANS->x1.x1,&A->x1.x1,&frobenius_2_4);
  fp2_set_neg(&ANS->x1.x1,&ANS->x1.x1);
}

void fp8_frobenius_map_p4(fp8_t *ANS,fp8_t *A){
  fp4_set(&ANS->x0,&A->x0);
  fp4_set_neg(&ANS->x1,&A->x1);
}

void fp8_frobenius_map_p1_montgomery(fp8_t *ANS, fp8_t *A){
  fp2_swap(&ANS->x0.x0,&A->x0.x0);
  
  fp2_swap(&ANS->x0.x1,&A->x0.x1);
  fp2_mul_lazy_montgomery(&ANS->x0.x1,&ANS->x0.x1,&frobenius_2_4MR);

  fp2_swap(&ANS->x1.x0,&A->x1.x0);
  fp2_mul_lazy_montgomery(&ANS->x1.x0,&ANS->x1.x0,&frobenius_1_4MR);

  fp2_swap(&ANS->x1.x1,&A->x1.x1);
  fp2_mul_lazy_montgomery(&ANS->x1.x1,&ANS->x1.x1,&frobenius_3_4MR);
}

void fp8_frobenius_map_p2_montgomery(fp8_t *ANS, fp8_t *A){
  fp2_set(&ANS->x0.x0,&A->x0.x0);
  
  fp2_set_neg_montgomery(&ANS->x0.x1,&A->x0.x1);

  fp2_mul_lazy_montgomery(&ANS->x1.x0,&A->x1.x0,&frobenius_2_4MR);

  fp2_mul_lazy_montgomery(&ANS->x1.x1,&A->x1.x1,&frobenius_2_4MR);
  fp2_set_neg_montgomery(&ANS->x1.x1,&ANS->x1.x1);
  }

void fp8_frobenius_map_p4_montgomery(fp8_t *ANS,fp8_t *A){
  fp4_set(&ANS->x0,&A->x0);
  fp4_set_neg_montgomery(&ANS->x1,&A->x1);
}

void fp8_mul_sparse_add(fp8_t *ANS,fp8_t *A,fp8_t *B){  //?000?? * ??????
  //fp2_mul*3 fp_mul*4 = 13 m
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);//?? ?? ?? 00 a+bβ+cγ+dβγ d=0 tmpA.x1.x1
  fp8_set(&tmp_B,B);//?? ?? ?? ?? e+fβ+gγ+hβγ

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  // fp4_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd       //here
//---------------------------------------------------------
  //(a + 0) * c + d
  fp2_mul(&tmp2_fp.x0, &tmp_A.x1.x0, &tmp_B.x1.x0);
  fp2_mul(&tmp2_fp.x1, &tmp_A.x1.x0, &tmp_B.x1.x1);
//---------------------------------------------------------
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd

}

void fp8_mul_sparse_dbl(fp8_t *ANS,fp8_t *A,fp8_t *B){  //??0?00 * ??????
  //fpmul*2 fp2mul*4 = 14 m prob because 2->6
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);//?? ?? 00 ?? a+b0+c0^2 tmpA.x1.x0
  fp8_set(&tmp_B,B);//?? ?? ?? ?? d+e0+f0^2
  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac //here
  // fp4_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  //---------------------------------------------------------
  //(0 + b)* (c + d)
  fp2_mul(&tmp2_fp.x0, &tmp_A.x1.x1, &tmp_B.x1.x0);
  fp2_mul(&tmp2_fp.x1, &tmp_A.x1.x1, &tmp_B.x1.x1);
  fp4_mul_base(&tmp2_fp, &tmp2_fp);
//---------------------------------------------------------
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd

}

void fp8_mul_sparse_add_lazy_montgomery(fp8_t *ANS,fp8_t *A,fp8_t *B){  //?000?? * ??????
  //fp2_mul*3 fp_mul*4 = 13 m
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);//?? ?? ?? 00 a+bβ+cγ+dβγ d=0 tmpA.x1.x1
  fp8_set(&tmp_B,B);//?? ?? ?? ?? e+fβ+gγ+hβγ

  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul_lazy_montgomery(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  // fp4_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd       //here
//---------------------------------------------------------
  //(a + 0) * c + d
  fp2_mul_lazy_montgomery(&tmp2_fp.x0, &tmp_A.x1.x0, &tmp_B.x1.x0);
  fp2_mul_lazy_montgomery(&tmp2_fp.x1, &tmp_A.x1.x0, &tmp_B.x1.x1);
//---------------------------------------------------------
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul_lazy_montgomery(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
}

void fp8_mul_sparse_dbl_lazy_montgomery(fp8_t *ANS,fp8_t *A,fp8_t *B){  //??0?00 * ??????
  //fpmul*2 fp2mul*4 = 14 m prob because 2->6
  static fp8_t tmp_A,tmp_B;
  fp8_set(&tmp_A,A);//?? ?? 00 ?? a+b0+c0^2 tmpA.x1.x0
  fp8_set(&tmp_B,B);//?? ?? ?? ?? d+e0+f0^2
  static fp4_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp4_mul_lazy_montgomery(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac //here
  // fp4_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  //---------------------------------------------------------
  //(0 + b)* (c + d)
  fp2_mul_lazy_montgomery(&tmp2_fp.x0, &tmp_A.x1.x1, &tmp_B.x1.x0);
  fp2_mul_lazy_montgomery(&tmp2_fp.x1, &tmp_A.x1.x1, &tmp_B.x1.x1);
  fp4_mul_base(&tmp2_fp, &tmp2_fp);
//---------------------------------------------------------
  fp4_mul_base(&tmp3_fp, &tmp2_fp);  //bdΘ^2
  fp4_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp4_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp4_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp4_mul_lazy_montgomery(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp4_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp4_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
}
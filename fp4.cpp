#include "fp4.h"
#include "fp2.h"

void fp4_init(fp4_t *A){
  fp2_init(&A->x0);
  fp2_init(&A->x1);
}

void fp4_printf(std::string str ,fp4_t *A){
  gmp_printf("%s(",str.c_str());
  fp2_printf("",&A->x0);
  gmp_printf(",");
  fp2_printf("",&A->x1);
  gmp_printf(")\n");
}

void fp4_println(std::string str ,fp4_t *A){
  gmp_printf("%s(",str.c_str());
  fp2_printf("",&A->x0);
  gmp_printf(",");
  fp2_printf("",&A->x1);
  gmp_printf(")\n");
}

void fpd4_printf(std::string str ,fpd4_t *A){
  gmp_printf("%s(",str.c_str());
  fpd2_printf("",&A->x0);
  gmp_printf(",");
  fpd2_printf("",&A->x1);
  gmp_printf(")");
}

void fpd4_println(std::string str ,fpd4_t *A){
  gmp_printf("%s(",str.c_str());
  fpd2_printf("",&A->x0);
  gmp_printf(",");
  fpd2_printf("",&A->x1);
  gmp_printf(")\n");
}

void fp4_printf_montgomery(std::string str, fp4_t *A) {
  gmp_printf("%s(", str.c_str());
  fp2_printf_montgomery("", &A->x0);
  gmp_printf(",");
  fp2_printf_montgomery("", &A->x1);
  gmp_printf(")\n");
}

void fp4_println_montgomery(std::string str, fp4_t *A) {
  gmp_printf("%s(", str.c_str());
  fp2_printf_montgomery("", &A->x0);
  gmp_printf(",");
  fp2_printf_montgomery("", &A->x1);
  gmp_printf(")\n");
}

void fp4_set(fp4_t *ANS,fp4_t *A){
  fp2_set(&ANS->x0,&A->x0);
  fp2_set(&ANS->x1,&A->x1);
}

void fpd4_set(fpd4_t *ANS, fpd4_t *A) {
  fpd2_set(&ANS->x0, &A->x0);
  fpd2_set(&ANS->x1, &A->x1);
}

void fp4_set_ui(fp4_t *ANS,unsigned long int UI){
  fp2_set_ui(&ANS->x0,UI);
  fp2_set_ui(&ANS->x1,0);
}

void fp4_set_ui_ui(fp4_t *ANS,unsigned long int UI){
  fp2_set_ui(&ANS->x0,UI);
  fp2_set_ui(&ANS->x1,UI);
}

void fp4_set_mpn(fp4_t *ANS,mp_limb_t *A){
  fp2_set_mpn(&ANS->x0,A);
  fp2_set_ui(&ANS->x1,0);
}

void fp4_set_neg(fp4_t *ANS,fp4_t *A){
  fp2_set_neg(&ANS->x0,&A->x0);
  fp2_set_neg(&ANS->x1,&A->x1);
}


void fp4_set_neg_montgomery(fp4_t *ANS,fp4_t *A){
  fp2_set_neg_montgomery(&ANS->x0,&A->x0);
  fp2_set_neg_montgomery(&ANS->x1,&A->x1);
}

void fp4_set_conj(fp4_t *ANS,fp4_t *A){
  fp2_set(&ANS->x0,&A->x0);
  fp2_set_neg(&ANS->x1,&A->x1);
}

void fp4_set_conj_montgomery(fp4_t *ANS,fp4_t *A){
  fp2_set(&ANS->x0,&A->x0);
  fp2_set_neg_montgomery(&ANS->x1,&A->x1);
}

void fp4_to_montgomery(fp4_t *ANS,fp4_t *A){
  fp2_to_montgomery(&ANS->x0,&A->x0);
  fp2_to_montgomery(&ANS->x1,&A->x1);
}

void fp4_mod_montgomery(fp4_t *ANS,fp4_t *A){
  mpn_mod_montgomery(ANS->x0.x0.x0,FPLIMB,A->x0.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x0.x1.x0,FPLIMB,A->x0.x1.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x1.x0.x0,FPLIMB,A->x1.x0.x0,FPLIMB);
  mpn_mod_montgomery(ANS->x1.x1.x0,FPLIMB,A->x1.x1.x0,FPLIMB);
}

void fp4_mod_montgomery_double(fp4_t *ANS,fpd4_t *A){
  mpn_mod_montgomery(ANS->x0.x0.x0,FPLIMB,A->x0.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x0.x1.x0,FPLIMB,A->x0.x1.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x1.x0.x0,FPLIMB,A->x1.x0.x0,FPLIMB2);
  mpn_mod_montgomery(ANS->x1.x1.x0,FPLIMB,A->x1.x1.x0,FPLIMB2);
}

void fp4_lshift(fp4_t *ANS, fp4_t *A, unsigned long int UI) {
  fp2_lshift(&ANS->x0, &A->x0, UI);
  fp2_lshift(&ANS->x1, &A->x1, UI);
}
void fp4_l1shift(fp4_t *ANS, fp4_t *A) {
  fp2_l1shift(&ANS->x0, &A->x0);
  fp2_l1shift(&ANS->x1, &A->x1);
}

void fp4_l1shift_nonmod_single(fp4_t *ANS, fp4_t *A) {
  fp2_l1shift_nonmod_single(&ANS->x0, &A->x0);
  fp2_l1shift_nonmod_single(&ANS->x1, &A->x1);
}

void fp4_l1shift_nonmod_double(fpd4_t *ANS, fpd4_t *A) {
  fp2_l1shift_nonmod_double(&ANS->x0, &A->x0);
  fp2_l1shift_nonmod_double(&ANS->x1, &A->x1);
}

void fp4_r1shift(fp4_t *ANS, fp4_t *A) {
  fp2_r1shift(&ANS->x0, &A->x0);
  fp2_r1shift(&ANS->x1, &A->x1);
}

void fp4_set_random(fp4_t *ANS,gmp_randstate_t state){
  fp2_set_random(&ANS->x0,state);
  fp2_set_random(&ANS->x1,state);
}

void fp4_mul(fp4_t *ANS,fp4_t *A,fp4_t *B){ 
  static fp4_t tmp_A,tmp_B;
  fp4_set(&tmp_A,A);
  fp4_set(&tmp_B,B);

  static fp2_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp2_mul(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  fp2_mul(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  fp2_mul_base(&tmp3_fp, &tmp2_fp);  //ab+bdΘ^2
  fp2_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp2_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp2_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp2_mul(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp2_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp2_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
} 

// void fp4_mul_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B){ 
//   static fp4_t tmp_A,tmp_B;
//   fp4_set(&tmp_A,A);
//   fp4_set(&tmp_B,B);

//   static fp2_t tmp3_fp,tmp4_fp;
//   static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
//   fp2_mul_nonmod(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
//   fp2_mul_nonmod(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
//   fp2_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
//   fp2_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2
//   fp2_mod(&ANS->x0,tmp4_fpd.x0,FPLIMB2);

//   fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
//   fp2_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
//   fp2_mul_nonmod(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
//   fp2_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
//   fp2_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
//   fp2_mod(&ANS->x1,tmp4_fpd.x0,FPLIMB2);
// } 

void fp4_mul_lazy_montgomery(fp4_t *ANS,fp4_t *A,fp4_t *B){ 
  // static fp4_t tmp_A,tmp_B;
  // fp4_set(&tmp_A,A);
  // fp4_set(&tmp_B,B);

  // static fp2_t tmp3_fp,tmp4_fp;
  // static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
  // fp2_mul_nonmod(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
  // fp2_mul_nonmod(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
  // fp2_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
  // fp2_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2
  // mpn_mod_montgomery(ANS->x0.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  // fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  // fp2_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  // fp2_mul_nonmod(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  // fp2_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
  // fp2_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
  // mpn_mod_montgomery(ANS->x1.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  static fp4_t tmp_A,tmp_B;
  fp4_set(&tmp_A,A);
  fp4_set(&tmp_B,B);

  static fp2_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp,tmp5_fp;
  fp2_mul_lazy_montgomery(&tmp1_fp,&tmp_A.x0,&tmp_B.x0); //ac
  fp2_mul_lazy_montgomery(&tmp2_fp,&tmp_A.x1,&tmp_B.x1); //bd
  fp2_mul_base(&tmp3_fp, &tmp2_fp);  //ab+bdΘ^2
  fp2_add(&ANS->x0, &tmp1_fp, &tmp3_fp);  //ab+bdΘ^2

  fp2_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp2_add(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp2_mul_lazy_montgomery(&tmp5_fp,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  
  fp2_sub(&tmp3_fp,&tmp5_fp,&tmp1_fp);//(a+b)(c+d) - ac
  fp2_sub(&ANS->x1,&tmp3_fp,&tmp2_fp);//(a+b)(c+d) - ac -bd
} 


void fp4_mul_nonmod_montgomery(fpd4_t *ANS, fp4_t *A, fp4_t *B) {
  static fp4_t tmp_A,tmp_B;
  fp4_set(&tmp_A,A);
  fp4_set(&tmp_B,B);

  static fp2_t tmp3_fp,tmp4_fp;
  static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd,tmp5_fpd;
  
  fp2_mul_nonmod_montgomery(&tmp1_fpd,&tmp_A.x0,&tmp_B.x0); //ac
  fp2_mul_nonmod_montgomery(&tmp2_fpd,&tmp_A.x1,&tmp_B.x1); //bd
  fp2_l1shift_nonmod_double(&tmp3_fpd, &tmp2_fpd);  //ab+bdΘ^2
  fp2_add_nonmod_double(&ANS->x0, &tmp1_fpd, &tmp3_fpd);  //ab+bdΘ^2

  fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1);//a+b
  fp2_add_nonmod_single(&tmp4_fp,&tmp_B.x0,&tmp_B.x1);//c+d
  fp2_mul_nonmod_montgomery(&tmp5_fpd,&tmp3_fp,&tmp4_fp); //(a+b)(c+d)
  fp2_sub_nonmod_double(&tmp3_fpd,&tmp5_fpd,&tmp1_fpd);//(a+b)(c+d) - ac
  fp2_sub_nonmod_double(&ANS->x1,&tmp3_fpd,&tmp2_fpd);//(a+b)(c+d) - ac -bd
}

void fp4_mul_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
  fp2_mul_ui(&ANS->x0,&A->x0,UI);
  fp2_mul_ui(&ANS->x1,&A->x1,UI);
}
void fp4_mul_ui_nonmod_single(fp4_t *ANS, fp4_t *A, unsigned long int UI) {
  fp2_mul_ui_nonmod_single(&ANS->x0, &A->x0, UI);
  fp2_mul_ui_nonmod_single(&ANS->x1, &A->x1, UI);
}
void fp4_mul_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
  fp2_mul_mpn(&ANS->x0,&A->x0,B);
  fp2_mul_mpn(&ANS->x1,&A->x1,B);
}

void fp4_mul_mpn_montgomery(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
  fp2_mul_mpn_montgomery(&ANS->x0, &A->x0, B);
  fp2_mul_mpn_montgomery(&ANS->x1, &A->x1, B);

}

void fp4_sqr(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  static fp2_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp2_sqr(&tmp1_fp, &tmp_A.x0);  //a^2
  fp2_sqr(&tmp2_fp, &tmp_A.x1);  //b^2
  fp2_mul_base(&tmp3_fp, &tmp2_fp);  //b^2@^2
  fp2_add(&ANS->x0, &tmp1_fp, &tmp3_fp); //a^2+b^2@^

  fp2_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp2_sqr(&tmp3_fp, &tmp3_fp);  //(a+b)^2
  fp2_sub(&tmp3_fp,&tmp3_fp, &tmp1_fp); //(a+b)^2 - a^2 
  fp2_sub(&ANS->x1,&tmp3_fp, &tmp2_fp); //(a+b)^2 - a^2 - b^2 

}

// void fp4_sqr_lazy(fp4_t *ANS,fp4_t *A){
//   static fp4_t tmp_A;
//   fp4_set(&tmp_A,A);
//   static fp2_t tmp3_fp;
//   static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
//   fp2_sqr_nonmod(&tmp1_fpd, &tmp_A.x0);  //a^2
//   fp2_sqr_nonmod(&tmp2_fpd, &tmp_A.x1);  //b^2

//   fp2_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
//   fp2_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^
//   fp2_mod(&ANS->x0,tmp4_fpd.x0,FPLIMB2);

//   fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
//   fp2_sqr_nonmod(&tmp3_fpd, &tmp3_fp);  //(a+b)^2

//   fp2_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
//   fp2_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 
//   fp2_mod(&ANS->x1,tmp4_fpd.x0,FPLIMB2);
// }

void fp4_sqr_lazy_montgomery(fp4_t *ANS,fp4_t *A){
  // static fp4_t tmp_A;
  // fp4_set(&tmp_A,A);
  // static fp2_t tmp3_fp;
  // static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
  // fp2_sqr_nonmod(&tmp1_fpd, &tmp_A.x0);  //a^2
  // fp2_sqr_nonmod(&tmp2_fpd, &tmp_A.x1);  //b^2
  // fp2_l1shift_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
  // fp2_add_nonmod_double(&tmp4_fpd, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^
  // mpn_mod_montgomery(ANS->x0.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);

  // fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  // fp2_sqr_nonmod(&tmp3_fpd, &tmp3_fp);  //(a+b)^2
  // fp2_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
  // fp2_sub_nonmod_double(&tmp4_fpd,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 
  // mpn_mod_montgomery(ANS->x1.x0,FPLIMB,tmp4_fpd.x0,FPLIMB2);
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  static fp2_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp2_sqr_lazy_montgomery(&tmp1_fp, &tmp_A.x0);  //a^2
  fp2_sqr_lazy_montgomery(&tmp2_fp, &tmp_A.x1);  //b^2
  fp2_mul_base(&tmp3_fp, &tmp2_fp);  //b^2@^2
  fp2_add(&ANS->x0, &tmp1_fp, &tmp3_fp); //a^2+b^2@^

  fp2_add(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp2_sqr_lazy_montgomery(&tmp3_fp, &tmp3_fp);  //(a+b)^2
  fp2_sub(&tmp3_fp,&tmp3_fp, &tmp1_fp); //(a+b)^2 - a^2 
  fp2_sub(&ANS->x1,&tmp3_fp, &tmp2_fp); //(a+b)^2 - a^2 - b^2 

}

void fp4_sqr_nonmod_montgomery(fpd4_t *ANS, fp4_t *A) {

  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);
  static fp2_t tmp3_fp;
  static fpd2_t tmp1_fpd,tmp2_fpd,tmp3_fpd,tmp4_fpd;
  
  fp2_sqr_nonmod_montgomery(&tmp1_fpd, &tmp_A.x0);  //a^2
  fp2_sqr_nonmod_montgomery(&tmp2_fpd, &tmp_A.x1);  //b^2
  fp2_l1shift_nonmod_double(&tmp3_fpd, &tmp2_fpd);  //b^2@^2
  fp2_add_nonmod_double(&ANS->x0, &tmp1_fpd, &tmp3_fpd); //a^2+b^2@^

  fp2_add_nonmod_single(&tmp3_fp,&tmp_A.x0,&tmp_A.x1); //(a+b)
  fp2_sqr_nonmod_montgomery(&tmp3_fpd, &tmp3_fp);  //(a+b)^2
  fp2_sub_nonmod_double(&tmp3_fpd,&tmp3_fpd, &tmp1_fpd); //(a+b)^2 - a^2 
  fp2_sub_nonmod_double(&ANS->x1,&tmp3_fpd, &tmp2_fpd); //(a+b)^2 - a^2 - b^2 

}

void fp4_sqr_GS(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  static fp2_t tmp1_fp,tmp2_fp,tmp3_fp;

  fp2_sqr(&tmp1_fp, &tmp_A.x0);  //a^2
  fp2_add(&tmp2_fp, &tmp_A.x0, &tmp_A.x1);  //(a+b)
  fp2_sqr(&tmp2_fp,&tmp2_fp);  //(a+b)^2

  fp2_sub_ui(&tmp3_fp,&tmp1_fp,1);//a^2-1
  fp2_add(&ANS->x0,&tmp3_fp,&tmp1_fp);//2a^2-1

  fp2_mul_base_inv(&tmp3_fp,&tmp3_fp);//(a^2-1)/2
  
  fp2_sub(&ANS->x1,&tmp2_fp,&tmp1_fp);//(a+b)^2 - a^2
  fp2_sub(&ANS->x1,&ANS->x1, &tmp3_fp);//(a+b)^2 - a^2 - (a^2-1)/2
}

void fp4_add(fp4_t *ANS,fp4_t *A,fp4_t *B){
  fp2_add(&ANS->x0,&A->x0,&B->x0);
  fp2_add(&ANS->x1,&A->x1,&B->x1);

}

void fp4_add_nonmod_single(fp4_t *ANS,fp4_t *A,fp4_t *B){
  fp2_add_nonmod_single(&ANS->x0,&A->x0,&B->x0);
  fp2_add_nonmod_single(&ANS->x1,&A->x1,&B->x1);
}

void fp4_add_nonmod_double(fpd4_t *ANS,fpd4_t *A,fpd4_t *B){
  fp2_add_nonmod_double(&ANS->x0,&A->x0,&B->x0);
  fp2_add_nonmod_double(&ANS->x1,&A->x1,&B->x1);

}

void fp4_add_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
  fp2_add_ui(&ANS->x0,&A->x0,UI);
  fp2_set(&ANS->x1,&A->x1);

}

void fp4_add_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
  fp2_add_ui(&ANS->x0,&A->x0,UI);
  fp2_add_ui(&ANS->x1,&A->x1,UI);

}

void fp4_add_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
  fp2_add_mpn(&ANS->x0,&A->x0,B);
  fp2_add_mpn(&ANS->x1,&A->x1,B);

}

void fp4_sub(fp4_t *ANS,fp4_t *A,fp4_t *B){
  fp2_sub(&ANS->x0,&A->x0,&B->x0);
  fp2_sub(&ANS->x1,&A->x1,&B->x1);

}

void fp4_sub_nonmod_single(fp4_t *ANS,fp4_t *A,fp4_t *B){
  fp2_sub_nonmod_single(&ANS->x0,&A->x0,&B->x0);
  fp2_sub_nonmod_single(&ANS->x1,&A->x1,&B->x1);

}

void fp4_sub_nonmod_double(fpd4_t *ANS,fpd4_t *A,fpd4_t *B){
  fp2_sub_nonmod_double(&ANS->x0,&A->x0,&B->x0);
  fp2_sub_nonmod_double(&ANS->x1,&A->x1,&B->x1);

;}

void fp4_sub_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
  fp2_sub_ui(&ANS->x0,&A->x0,UI);
  fp2_set(&ANS->x1,&A->x1);

}

void fp4_sub_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
  fp2_sub_ui(&ANS->x0,&A->x0,UI);
  fp2_sub_ui(&ANS->x1,&A->x1,UI);

}

void fp4_sub_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
  fp2_sub_mpn(&ANS->x0,&A->x0,B);
  fp2_sub_mpn(&ANS->x1,&A->x1,B);

}

void fp4_inv(fp4_t *ANS,fp4_t *A){  
	static fp2_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
  fp2_set(&tmp1_fp,&A->x0);  //a
  fp2_set_neg(&tmp2_fp,&A->x1); //b

  fp2_sqr(&tmp3_fp,&tmp1_fp);  //a^2
  fp2_mul(&tmp4_fp,&tmp2_fp,&A->x1);  //b^2*c
  fp2_mul_base(&tmp4_fp,&tmp4_fp);
  fp2_add(&tmp3_fp,&tmp3_fp,&tmp4_fp); //a^2 - b^2c //mabe need mul_basis?

  fp2_inv(&tmp3_fp,&tmp3_fp);  //  (a^2 - b^2c)^-1
  fp2_mul(&ANS->x0,&tmp1_fp,&tmp3_fp); // a*(a^2- b^2c)^-1
  fp2_mul(&ANS->x1,&tmp2_fp,&tmp3_fp); //-b*(a^2 - b^2c)^-1
}

void fp4_inv_lazy(fp4_t *ANS, fp4_t *A) {
	static fp2_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
  fp2_set(&tmp1_fp,&A->x0);  //a
  fp2_set_neg(&tmp2_fp,&A->x1); //b

  fp2_sqr(&tmp3_fp,&tmp1_fp);  //a^2
  fp2_mul(&tmp4_fp,&tmp2_fp,&A->x1);  //b^2*c
  fp2_mul_base(&tmp4_fp,&tmp4_fp);
  fp2_add(&tmp3_fp,&tmp3_fp,&tmp4_fp); //a^2 - b^2c //mabe need mul_basis?

  fp2_inv(&tmp3_fp,&tmp3_fp);  //  (a^2 - b^2c)^-1
  fp2_mul(&ANS->x0,&tmp1_fp,&tmp3_fp); // a*(a^2- b^2c)^-1
  fp2_mul(&ANS->x1,&tmp2_fp,&tmp3_fp); //-b*(a^2 - b^2c)^-1
}

void fp4_inv_lazy_montgomery(fp4_t *ANS, fp4_t *A) {
  static fp2_t tmp1_fp, tmp2_fp;
  static fp2_t tmp3;
  static fp2_t tmp1, tmp2;
  fp2_set(&tmp1_fp, &A->x0);
  fp2_set_neg(&tmp2_fp, &A->x1);

  fp2_sqr_lazy_montgomery(&tmp1, &tmp1_fp); //a^2
  fp2_mul_lazy_montgomery(&tmp2, &tmp2_fp, &A->x1);//b^2
  fp2_mul_base(&tmp2,&tmp2);//b^2*c
  fp2_add(&tmp3, &tmp1, &tmp2);//a^2-b^2*c

  fp2_inv_lazy_montgomery(&tmp3, &tmp3);
  fp2_mul_lazy_montgomery(&ANS->x0, &tmp1_fp, &tmp3);
  fp2_mul_lazy_montgomery(&ANS->x1, &tmp2_fp, &tmp3);
}

int fp4_legendre(fp4_t *A){
  fp4_t tmp;
  mpz_t expo;
  fp4_init(&tmp);
  mpz_init(expo);

  //(p^2 -1)/2 を計算
  mpz_pow_ui(expo,prime_z,4);
  mpz_sub_ui(expo,expo,1);
  mpz_tdiv_q_ui(expo,expo,2);
  fp4_pow(&tmp,A,expo);

  if(fp4_cmp_one(&tmp)==0){
    mpz_clear(expo);
    return 1;
  }
  else{
     mpz_clear(expo);
    return -1;
  }
}

void fp4_sqrt(fp4_t *ANS,fp4_t *A){
  fp4_t x,y,t,k,n,tmp;
  fp4_init(&x);
  fp4_init(&y);
  fp4_init(&t);
  fp4_init(&k);
  fp4_init(&n);
  fp4_init(&tmp);
  unsigned long int e,m;
  mpz_t exp,q,z,result;
  mpz_init(exp);
  mpz_init(q);
  mpz_init(z);
  mpz_init(result);
  //gmp_randstate_t state;
  //gmp_randinit_default(state);
  //gmp_randseed_ui(state,(unsigned long)time(NULL));

  fp4_set_random(&n,state);
  while(fp4_legendre(&n)!=-1){
    fp4_set_random(&n,state);
  }
  mpz_pow_ui(q,prime_z,4);
  mpz_sub_ui(q,q,1);
  mpz_mod_ui(result,q,2);
  e=0;
  while(mpz_cmp_ui(result,0)==0){
    mpz_tdiv_q_ui(q,q,2);
    mpz_mod_ui(result,q,2);
    e++;
  }
  fp4_pow(&y,&n,q);
  mpz_set_ui(z,e);
  mpz_sub_ui(exp,q,1);
  mpz_tdiv_q_ui(exp,exp,2);
  fp4_pow(&x,A,exp);
  fp4_mul(&tmp,&x,&x);
  fp4_mul(&k,&tmp,A);
  fp4_mul(&x,&x,A);
  while(fp4_cmp_one(&k)!=0){
    m=1;
    mpz_ui_pow_ui(exp,2,m);
    fp4_pow(&tmp,&k,exp);
    while(fp4_cmp_one(&tmp)!=0){
      m++;
      mpz_ui_pow_ui(exp,2,m);
      fp4_pow(&tmp,&k,exp);
    }
    mpz_sub_ui(exp,z,m);
    mpz_sub_ui(exp,exp,1);
    mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
    fp4_pow(&t,&y,result);
    fp4_mul(&y,&t,&t);
    mpz_set_ui(z,m);
    fp4_mul(&x,&x,&t);
    fp4_mul(&k,&k,&y);
  }
  fp4_set(ANS,&x);

  mpz_clear(exp);
  mpz_clear(q);
  mpz_clear(z);
  mpz_clear(result);
}

void fp4_pow(fp4_t *ANS,fp4_t *A,mpz_t scalar){
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);
  fp4_t tmp;
  fp4_init(&tmp);
  fp4_set(&tmp,A);

  for(i=1;i<length; i++){
    fp4_sqr(&tmp,&tmp);
    if(binary[i]=='1')  fp4_mul(&tmp,A,&tmp);
  }
  fp4_set(ANS,&tmp);
}

void fp4_pow_montgomery(fp4_t *ANS, fp4_t *A, mpz_t scalar) {
  int length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);
  fp4_t tmp;
  fp4_init(&tmp); // not need?

  fp4_set(&tmp, A);

  for (int i = 1; i < length; i++) {
    fp4_sqr_lazy_montgomery(&tmp, &tmp);
    if (binary[i] == '1') {
      fp4_mul_lazy_montgomery(&tmp, A, &tmp);
    }
  }
  fp4_set(ANS, &tmp);
}

int fp4_cmp(fp4_t *A,fp4_t *B){
  if(fp2_cmp(&A->x0,&B->x0)==0 && fp2_cmp(&A->x1,&B->x1)==0){
    return 0;
  }
  return 1;
}

int fp4_cmp_ui(fp4_t *A,unsigned long int UI){
  if(fp2_cmp_ui(&A->x0,UI)==0 && fp2_cmp_ui(&A->x1,UI)==0){
    return 0;
  }
  return 1;
}

int fp4_cmp_mpn(fp4_t *A,mp_limb_t *B){
  if(fp2_cmp_mpn(&A->x0,B)==0 && fp2_cmp_mpn(&A->x1,B)==0){
    return 0;
  }
  return 1;
}

int fp4_cmp_zero(fp4_t *A){
  if(fp2_cmp_zero(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0 ){
    return 0;
  }
  return 1;
}

int fp4_cmp_one(fp4_t *A){
  if(fp2_cmp_one(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0 ){
    return 0;
  }
  return 1;
}

// void fp4_lshift_ui_nonmod_single(fp4_t *ANS, fp4_t *A, int s) {
//   fp2_lshift_ui_nonmod_single(&ANS->x0, &A->x0, s);
//   fp2_lshift_ui_nonmod_single(&ANS->x1, &A->x1, s);
// }
// void fp4_lshift_ui_nonmod_double(fpd4_t *ANS, fpd4_t *A, int s) {
//   fp2_lshift_ui_nonmod_double(&ANS->x0, &A->x0, s);
//   fp2_lshift_ui_nonmod_double(&ANS->x1, &A->x1, s);
// }

void fp4_frobenius_map_p1(fp4_t *ANS,fp4_t *A){
  fp2_set(&ANS->x0,&A->x0);
  fp2_set_neg(&ANS->x1,&A->x1);
}

void fp4_mul_base(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  fp2_mul_base(&ANS->x0, &A->x1);
  fp2_set(&ANS->x1,&tmp_A.x0);    //@^2 = 2
}

void fp4_mul_base_1(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  fp2_mul_base(&ANS->x0, &A->x1); //2lshift,2sub
  fp2_set(&ANS->x1,&tmp_A.x0);    //@^2 = 2//2set

  fp4_add(ANS,ANS,A); //4 fp_add
}

void fp4_mul_base_inv(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  fp2_mul_base_inv(&ANS->x1, &tmp_A.x0);
  fp2_set(&ANS->x0,&tmp_A.x1);    //@^2 = 2
}

void fp4_mul_base_inv_montgomery(fp4_t *ANS,fp4_t *A){
  static fp4_t tmp_A;
  fp4_set(&tmp_A,A);

  fp2_mul_base_inv_montgomery(&ANS->x1, &tmp_A.x0);
  fp2_set(&ANS->x0,&tmp_A.x1);    //@^2 = 2
}
#include "miller.h"
#include "define.h"
#include "fp2.h"
#include "fp8.h"
#include "mpn.h"
#include <cstdio>


void miller_proj_precomp_Costello(efp2_t *S,efp2_t *R){
  //call this after generate P & Q_dash precomputed value wont change

}

//twist
void efp8_to_Jacefp2(efp2_jacobian_t *ANS,efp8_t *A){ 
  fp2_mul_base(&ANS->x,&A->x.x0.x1);
  fp2_mul_base(&ANS->y,&A->y.x1.x0);
  fp2_set_ui(&ANS->z, 1);
  ANS->infinity = 0;
}

void efp8_to_efp2(efp2_t *ANS,efp8_t *A){ 
  fp2_mul_base(&ANS->x,&A->x.x0.x1);
  fp2_mul_base(&ANS->y,&A->y.x1.x0);
  ANS->infinity = A->infinity;
}

void efp8_to_Jacefp2_montgomery(efp2_jacobian_t *ANS,efp8_t *A){ 
  fp2_mul_base(&ANS->x,&A->x.x0.x1);
  fp2_to_montgomery(&ANS->x,&ANS->x);
  fp2_mul_base(&ANS->y,&A->y.x1.x0);
  fp2_to_montgomery(&ANS->y,&ANS->y);
  fp2_set_ui(&ANS->z, 1);
  fp2_to_montgomery(&ANS->z,&ANS->z);
  ANS->infinity = 0;
}

void efp8_to_efp2_montgomery(efp2_t *ANS,efp8_t *A){ 
  fp2_mul_base(&ANS->x,&A->x.x0.x1);
  fp2_to_montgomery(&ANS->x,&ANS->x);
  fp2_mul_base(&ANS->y,&A->y.x1.x0);
  fp2_to_montgomery(&ANS->y,&ANS->y);
  ANS->infinity = A->infinity;
}

//double line //
void ff_lttp(fp8_t *f, efp2_jacobian_t *S, efp_t *P){
  fp8_sqr(f,f); //update

  static fp2_t tmp1_fp, tmp2_fp,tmp3_fp;
  static fp2_t t1,t2,t3;
  static fp2_t nextX,nextY,nextZ;

  static fp8_t tmp1_fp8;

  fp2_sqr(&t1,&S->y);              //t1 = Y^2
  fp2_l1shift(&tmp2_fp,&t1);      //tmp2 = 2*t1

  fp2_mul(&t2,&tmp2_fp,&S->x);     //t2 = 2*t1*X
  fp2_l1shift(&t2, &t2);          //t2 = 4*t1*X

  fp2_sqr(&t3,&S->x);              //t3 = X^2
  fp2_l1shift(&tmp1_fp, &t3);
  fp2_add(&t3,&t3,&tmp1_fp);       //t3 = 3X^2

  fp2_sqr(&tmp3_fp,&S->z);          //Z2=Z^2
  fp2_sqr(&tmp1_fp,&tmp3_fp);       //Z2^2

  // fp2_mul_mpn(&tmp1_fp, &tmp1_fp, curve_a.x0);
  fp2_mul_base(&tmp1_fp, &tmp1_fp);
  fp2_add(&t3,&t3,&tmp1_fp);    //t3 = 3X^2+aZ^2*alpha

  fp2_sqr(&nextX,&t3);             //nextX = t3^2
  fp2_l1shift(&tmp1_fp, &t2);     //tmp1 = 2*t2
  fp2_sub(&nextX,&nextX,&tmp1_fp); //nextX = t3^2 - 2*t2

  fp2_sub(&nextY,&t2,&nextX);      //nextY = t2-nextX
  fp2_mul(&nextY,&nextY,&t3);      //nextY = (t2-nextX)t3

  fp2_sqr(&tmp1_fp,&tmp2_fp);      //tmp1 = tmp2^2 = 4t1^2
  fp2_l1shift(&tmp1_fp,&tmp1_fp); //tmp1 = 8t1^2

  fp2_sub(&nextY,&nextY,&tmp1_fp); //nextY = (t2-nextX)t3 - 8t1^2

  fp2_l1shift(&nextZ,&S->y);      //nextZ = 2Y
  fp2_mul(&nextZ,&nextZ,&S->z);    //nextZ = 2YZ
//------------------------------------
  fp2_mul(&tmp1_fp8.x1.x1,&nextZ,&tmp3_fp);           // = nextZ*Z^2
  fp2_mul_mpn(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,P->y.x0);     // = nextZ*Z^2Py
//------------------------------------
  fp2_mul_mpn(&tmp1_fp8.x0.x1,&t3,P->x.x0);                 // = t3*Px
  fp2_mul(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,&tmp3_fp);  // = t3*Px*Z^2
  fp2_set_neg(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1);       // = -t3*Px*Z^2
//------------------------------------
  fp2_mul(&tmp3_fp,&t3,&S->x);                        //tmp3 = t3*X
  fp2_sub(&tmp1_fp8.x0.x0,&tmp3_fp,&tmp2_fp);         // = t3*X - 2*t1
//------------------------------------
  fp8_mul_sparse_dbl(f,&tmp1_fp8,f);

  fp2_set(&S->x,&nextX);
  fp2_set(&S->y,&nextY);
  fp2_set(&S->z,&nextZ);
}

//double line //
void ff_lttp_lazy_montgomery(fp8_t *f, efp2_jacobian_t *S, efp_t *P){
  fp8_sqr_lazy_montgomery(f,f); //update

  static fp2_t tmp1_fp, tmp2_fp,tmp3_fp;
  static fp2_t t1,t2,t3;
  static fp2_t nextX,nextY,nextZ;

  static fp8_t tmp1_fp8;

  fp2_sqr_lazy_montgomery(&t1,&S->y);              //t1 = Y^2
  fp2_l1shift_nonmod_single(&tmp2_fp,&t1);      //tmp2 = 2*t1

  fp2_mul_lazy_montgomery(&t2,&tmp2_fp,&S->x);     //t2 = 2*t1*X
  fp2_l1shift_nonmod_single(&t2, &t2);          //t2 = 4*t1*X

  fp2_sqr_lazy_montgomery(&t3,&S->x);              //t3 = X^2
  fp2_l1shift_nonmod_single(&tmp1_fp, &t3);
  fp2_add_nonmod_single(&t3,&t3,&tmp1_fp);       //t3 = 3X^2

  fp2_sqr_lazy_montgomery(&tmp3_fp,&S->z);          //Z2=Z^2
  fp2_sqr_lazy_montgomery(&tmp1_fp,&tmp3_fp);       //Z2^2

  // fp2_mul_mpn_montgomery(&tmp1_fp, &tmp1_fp, curve_a.x0);
  fp2_mul_base(&tmp1_fp, &tmp1_fp);
  fp2_add_nonmod_single(&t3,&t3,&tmp1_fp);    //t3 = 3X^2+aZ^2*alpha

  fp2_sqr_lazy_montgomery(&nextX,&t3);             //nextX = t3^2
  fp2_l1shift_nonmod_single(&tmp1_fp, &t2);     //tmp1 = 2*t2
  fp2_sub_nonmod_single(&nextX,&nextX,&tmp1_fp); //nextX = t3^2 - 2*t2

  fp2_sub_nonmod_single(&nextY,&t2,&nextX);      //nextY = t2-nextX
  fp2_mul_lazy_montgomery(&nextY,&nextY,&t3);      //nextY = (t2-nextX)t3

  fp2_sqr_lazy_montgomery(&tmp1_fp,&tmp2_fp);      //tmp1 = tmp2^2 = 4t1^2
  fp2_l1shift_nonmod_single(&tmp1_fp,&tmp1_fp); //tmp1 = 8t1^2

  fp2_sub_nonmod_single(&nextY,&nextY,&tmp1_fp); //nextY = (t2-nextX)t3 - 8t1^2

  fp2_l1shift_nonmod_single(&nextZ,&S->y);      //nextZ = 2Y
  fp2_mul_lazy_montgomery(&nextZ,&nextZ,&S->z);    //nextZ = 2YZ
//------------------------------------
  fp2_mul_lazy_montgomery(&tmp1_fp8.x1.x1,&nextZ,&tmp3_fp);           // = nextZ*Z^2
  fp2_mul_mpn_montgomery(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,P->y.x0);     // = nextZ*Z^2Py
//------------------------------------
  fp2_mul_mpn_montgomery(&tmp1_fp8.x0.x1,&t3,P->x.x0);                 // = t3*Px
  fp2_mul_lazy_montgomery(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,&tmp3_fp);  // = t3*Px*Z^2
  fp2_set_neg_montgomery(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1);       // = -t3*Px*Z^2
//------------------------------------
  fp2_mul_lazy_montgomery(&tmp3_fp,&t3,&S->x);                        //tmp3 = t3*X
  fp2_sub_nonmod_single(&tmp1_fp8.x0.x0,&tmp3_fp,&tmp2_fp);         // = t3*X - 2*t1
//------------------------------------
  fp8_mul_sparse_dbl_lazy_montgomery(f,&tmp1_fp8,f);

  fp2_set(&S->x,&nextX);
  fp2_set(&S->y,&nextY);
  fp2_set(&S->z,&nextZ);
}


void ff_lttp_Costello(fp8_t *f, efp2_jacobian_t *U, efp_t *S){
  fp8_sqr(f,f); //update

  static fp2_t tmpA_fp2, tmpB_fp2,tmpC_fp2,tmpD_fp2,tmpE_fp2,tmpF_fp2;
  // static fp2_t t1,t2,t3;
  static fp2_t nextX,nextY,nextZ;
  static fp2_t tmp1_fp2,tmp2_fp2;
  static fp8_t tmp1_fp8;

  fp2_sqr(&tmpA_fp2,&U->x);              //A = X^2
  fp2_sqr(&tmpB_fp2,&U->y);              //B = Y^2
  fp2_sqr(&tmpC_fp2,&U->z);              //C = Z^2

  // fp2_mul_mpn(&tmpD_fp2, &tmpC_fp2, curve_a.x0);  //D=aC
  fp2_mul_base(&tmpD_fp2, &tmpD_fp2);

  fp2_sub(&tmp1_fp2,&tmpA_fp2,&tmpD_fp2);  //X3=(A-D)^2
  fp2_sqr(&nextX,&tmp1_fp2);

  fp2_add(&tmp2_fp2,&tmpA_fp2,&tmpD_fp2);  //E=(A+D)^2
  fp2_sqr(&tmpE_fp2,&tmp2_fp2);
  fp2_l1shift(&tmpE_fp2,&tmpE_fp2);
  fp2_sub(&tmpE_fp2,&tmpE_fp2,&nextX);  //E=2(A+D)^2 - X3

  fp2_add(&tmpF_fp2,&tmp1_fp2,&U->y);   //F = (A-D+Y1)
  fp2_sqr(&tmpF_fp2,&tmpF_fp2);         //F = (A-D+Y1)^2
  fp2_sub(&tmpF_fp2,&tmpF_fp2,&tmpB_fp2);  //F=(A-D+Y1)^2 - B
  fp2_sub(&tmpF_fp2,&tmpF_fp2,&nextX);  //F=(A-D+Y1)^2 - B - X3

  fp2_mul(&nextY,&tmpE_fp2 ,&tmpF_fp2); //Y3 = E*F

  fp2_l1shift(&nextZ,&tmpB_fp2);
  fp2_l1shift(&nextZ,&nextZ);           //Z3 = 4B

//------------------------------------
  // fp2_mul(&tmp3_fp,&t3,&S->x);                        //tmp3 = t3*X
  // fp2_sub(&tmp1_fp8.x0.x0,&tmp3_fp,&tmp2_fp);         // = t3*X - 2*t1
  fp2_add(&tmp1_fp8.x0.x0,&U->x,&tmp1_fp2);
  fp2_sqr(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0);
  fp2_sub(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0,&nextX);
  fp2_sub(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0,&tmpA_fp2);

//------------------------------------
  // fp2_mul(&tmp1_fp8.x1.x1,&nextZ,&tmp3_fp);           // = nextZ*Z^2
  // fp2_mul_mpn(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,P->y.x0);     // = nextZ*Z^2Py
  fp2_add(&tmp1_fp8.x1.x1,&U->y,&U->z);
  fp2_sqr(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1);
  fp2_sub(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,&tmpB_fp2);
  fp2_sub(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,&tmpC_fp2);
  fp2_l1shift(&tmp1_fp8.x1.x1, &tmp1_fp8.x1.x1);
  fp2_mul_mpn(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,S->y.x0);     // = nextZ*Z^2Py

//------------------------------------
  // fp2_mul(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,&tmp3_fp);  // = t3*Px*Z^2
  fp2_l1shift(&tmp1_fp2, &tmpA_fp2);
  fp2_add(&tmp1_fp8.x0.x1,&tmp1_fp2,&tmp2_fp2);
  fp2_mul(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,&U->z);
  fp2_l1shift(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1);
  fp2_mul_mpn(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,S->x.x0);                 // = t3*Px
  fp2_set_neg(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1);       // = -t3*Px*Z^2

//------------------------------------
  fp8_mul_sparse_dbl(f,&tmp1_fp8,f);

  fp2_set(&U->x,&nextX);
  fp2_set(&U->y,&nextY);
  fp2_set(&U->z,&nextZ);
  // efp2_proj_w1_2_checkOnCurve_Twist(S);
}

//add line 
void ff_ltqp(fp8_t *f, efp2_jacobian_t *S, efp2_t *Q,efp_t *P){
  static fp2_t tmp1_fp, tmp2_fp,tmp3_fp;
  static fp2_t t1,t2,t3,t4,t5;
  static fp2_t nextX,nextY,nextZ;

  static fp8_t tmp1_fp8;

  fp2_sqr(&tmp1_fp,&S->z);
  fp2_mul(&t1,&Q->x,&tmp1_fp);
  fp2_sub(&t1,&t1,&S->x);
  //t1 = (Z^2)*Qx - X

  fp2_mul(&tmp1_fp,&tmp1_fp,&S->z);
  fp2_mul(&t2,&tmp1_fp,&Q->y);
  fp2_sub(&t2,&t2,&S->y);
  //t2 = Z^3 * Qy - Y

  fp2_sqr(&t3,&t1);
  //t3 = t1^2 = ((Z^2)*Qx - X)^2

  fp2_mul(&t4,&t1,&t3);
  //t4 = ((Z^2)*Qx - X)^3

  fp2_mul(&t5,&S->x,&t3);
  //t5 = X((Z^2)*Qx - X)^2

  fp2_sqr(&tmp1_fp,&t2);
  // fp2_mul_ui(&tmp2_fp,&t5,2);
  fp2_l1shift(&tmp2_fp, &t5);
  fp2_add(&tmp2_fp,&tmp2_fp,&t4);
  fp2_sub(&nextX,&tmp1_fp,&tmp2_fp);
  //X = t2^2 -(2t5+t4)

  fp2_sub(&nextY,&t5,&nextX);
  fp2_mul(&nextY,&nextY,&t2);
  fp2_mul(&tmp1_fp,&t4,&S->y);
  fp2_sub(&nextY,&nextY,&tmp1_fp);
  //Y = (t5 - X)t2 - t4Y

  fp2_mul(&nextZ,&S->z,&t1);
//------------------------------------
  fp2_mul_mpn(&tmp1_fp8.x1.x0,&nextZ,P->y.x0);//Zt1*py
//------------------------------------
  fp2_mul_mpn(&tmp1_fp8.x0.x0,&t2,P->x.x0);     //t2xp
  fp2_set_neg(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0); //-t2xp
//------------------------------------
  fp2_mul(&tmp1_fp,&t2,&Q->x);      //t2Qx
  fp2_mul(&tmp2_fp,&nextZ,&Q->y);   //Zt1Qy
  fp2_sub(&tmp1_fp8.x0.x1,&tmp1_fp,&tmp2_fp);
  fp2_mul_base_inv(&tmp1_fp8.x0.x1, &tmp1_fp8.x0.x1);
//------------------------------------
  fp8_mul_sparse_add(f,&tmp1_fp8,f);

  fp2_set(&S->x,&nextX);
  fp2_set(&S->y,&nextY);
  fp2_set(&S->z,&nextZ);
}

//add line 
void ff_ltqp_lazy_montgomery(fp8_t *f, efp2_jacobian_t *S, efp2_t *Q,efp_t *P){
  static fp2_t tmp1_fp, tmp2_fp,tmp3_fp;
  static fp2_t t1,t2,t3,t4,t5;
  static fp2_t nextX,nextY,nextZ;

  static fp8_t tmp1_fp8;

  fp2_sqr_lazy_montgomery(&tmp1_fp,&S->z);
  fp2_mul_lazy_montgomery(&t1,&Q->x,&tmp1_fp);
  fp2_sub_nonmod_single(&t1,&t1,&S->x);
  //t1 = (Z^2)*Qx - X

  fp2_mul_lazy_montgomery(&tmp1_fp,&tmp1_fp,&S->z);
  fp2_mul_lazy_montgomery(&t2,&tmp1_fp,&Q->y);
  fp2_sub_nonmod_single(&t2,&t2,&S->y);
  //t2 = Z^3 * Qy - Y

  fp2_sqr_lazy_montgomery(&t3,&t1);
  //t3 = t1^2 = ((Z^2)*Qx - X)^2

  fp2_mul_lazy_montgomery(&t4,&t1,&t3);
  //t4 = ((Z^2)*Qx - X)^3

  fp2_mul_lazy_montgomery(&t5,&S->x,&t3);
  //t5 = X((Z^2)*Qx - X)^2

  fp2_sqr_lazy_montgomery(&tmp1_fp,&t2);
  // fp2_mul_ui(&tmp2_fp,&t5,2);
  fp2_l1shift_nonmod_single(&tmp2_fp, &t5);
  fp2_add_nonmod_single(&tmp2_fp,&tmp2_fp,&t4);
  fp2_sub_nonmod_single(&nextX,&tmp1_fp,&tmp2_fp);
  //X = t2^2 -(2t5+t4)

  fp2_sub_nonmod_single(&nextY,&t5,&nextX);
  fp2_mul_lazy_montgomery(&nextY,&nextY,&t2);
  fp2_mul_lazy_montgomery(&tmp1_fp,&t4,&S->y);
  fp2_sub_nonmod_single(&nextY,&nextY,&tmp1_fp);
  //Y = (t5 - X)t2 - t4Y

  fp2_mul_lazy_montgomery(&nextZ,&S->z,&t1);
//------------------------------------
  fp2_mul_mpn_montgomery(&tmp1_fp8.x1.x0,&nextZ,P->y.x0);//Zt1*py
//------------------------------------
  fp2_mul_mpn_montgomery(&tmp1_fp8.x0.x0,&t2,P->x.x0);     //t2xp
  fp2_set_neg_montgomery(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0); //-t2xp
//------------------------------------
  fp2_mul_lazy_montgomery(&tmp1_fp,&t2,&Q->x);      //t2Qx
  fp2_mul_lazy_montgomery(&tmp2_fp,&nextZ,&Q->y);   //Zt1Qy
  fp2_sub_nonmod_single(&tmp1_fp8.x0.x1,&tmp1_fp,&tmp2_fp);
  fp2_mul_base_inv(&tmp1_fp8.x0.x1, &tmp1_fp8.x0.x1);
//------------------------------------
  fp8_mul_sparse_add_lazy_montgomery(f,&tmp1_fp8,f);

  fp2_set(&S->x,&nextX);
  fp2_set(&S->y,&nextY);
  fp2_set(&S->z,&nextZ);
}

// void ff_ltqp_Costello_mixed(fp8_t *f, efp2_jacobian_t *U, efp2_t *R,efp_t *S){

//   static fp2_t tmpA_fp2, tmpB_fp2,tmpC_fp2,tmpD_fp2,tmpE_fp2,tmpF_fp2;
//   static fp2_t tmpG_fp2, tmpH_fp2,tmpI_fp2,tmpI2_fp2,tmpJ_fp2,tmpK_fp2;

//   // static fp2_t t1,t2,t3;
//   static fp2_t nextX,nextY,nextZ;
//   static fp2_t tmp1_fp2,tmp2_fp2;
//   static fp8_t tmp1_fp8;


// //------------------------------------
//   fp2_mul_mpn(&tmp1_fp8.x1.x0,&nextZ,P->y.x0);//Zt1*py
// //------------------------------------
//   fp2_mul_mpn(&tmp1_fp8.x0.x0,&t2,P->x.x0);     //t2xp
//   fp2_set_neg(&tmp1_fp8.x0.x0,&tmp1_fp8.x0.x0); //-t2xp
// //------------------------------------
//   fp2_mul(&tmp1_fp,&t2,&Q->x);      //t2Qx
//   fp2_mul(&tmp2_fp,&nextZ,&Q->y);   //Zt1Qy
//   fp2_sub(&tmp1_fp8.x0.x1,&tmp1_fp,&tmp2_fp);
//   fp2_mul_base_inv(&tmp1_fp8.x0.x1, &tmp1_fp8.x0.x1);
// //------------------------------------
//   fp8_mul_sparse_add(f,&tmp1_fp8,f);

//   fp2_set(&S->x,&nextX);
//   fp2_set(&S->y,&nextY);
//   fp2_set(&S->z,&nextZ);
//   efp2_proj_w1_2_checkOnCurve_Twist(S);
// }

// void ff_lttp_lazy_montgomery(fp8_t *f, efp2_jacobian_t *S, efp_t *P){
//   fp8_sqr_lazy_montgomery(f,f); //update
//   static fp_t tmp1_fp, tmp2_fp,tmp3_fp;
//   static fp_t t1,t2,t3;
//   static fp_t nextX,nextY,nextZ;

//   static fp8_t tmp1_fp8;

//   fp_sqrmod_montgomery(&t1,&S->y);              //t1 = Y^2
//   fp_l1shift_nonmod_single(&tmp2_fp,&t1);      //tmp2 = 2*t1

//   fp_mulmod_montgomery(&t2,&tmp2_fp,&S->x);     //t2 = 2*t1*X
//   fp_l1shift_nonmod_single(&t2, &t2);          //t2 = 4*t1*X

//   fp_sqrmod_montgomery(&t3,&S->x);              //t3 = X^2
//   fp_l1shift_nonmod_single(&tmp1_fp, &t3);
//   fp_add_nonmod_single(&t3,&t3,&tmp1_fp);

//   fp_sqrmod_montgomery(&nextX,&t3);             //nextX = t3^2
//   fp_l1shift_nonmod_single(&tmp1_fp, &t2);             //tmp1 = 2*t2
//   fp_sub_nonmod_single(&nextX,&nextX,&tmp1_fp); //nextX = t3^2 - 2*t2

//   fp_sub_nonmod_single(&nextY,&t2,&nextX);      //nextY = t2-nextX
//   fp_mulmod_montgomery(&nextY,&nextY,&t3);      //nextY = (t2-nextX)t3

//   fp_sqrmod_montgomery(&tmp1_fp,&tmp2_fp);      //tmp1 = tmp2^2 = 4t1^2
//   fp_l1shift_nonmod_single(&tmp1_fp,&tmp1_fp); //tmp1 = 8t1^2

//   fp_sub_nonmod_single(&nextY,&nextY,&tmp1_fp); //nextY = (t2-nextX)t3 - 8t1^2

//   fp_l1shift_nonmod_single(&nextZ,&S->y);      //nextZ = 2Y
//   fp_mulmod_montgomery(&nextZ,&nextZ,&S->z);    //nextZ = 2YZ

//   fp_sqrmod_montgomery(&tmp1_fp,&S->z);         //tmp1 = Z^2
//   fp_mulmod_montgomery(&tmp1_fp8.x0.x1,&nextZ,&tmp1_fp);           // = nextZ*Z^2
//   fp_mul_base_inv_single(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1);  // = nextZ*Z^2*(2^-1)
//   fp_mulmod_montgomery(&tmp1_fp8.x0.x1,&tmp1_fp8.x0.x1,&P->y);     // = nextZ*Z^2*(2^-1)*Py

//   fp_mulmod_montgomery(&tmp3_fp,&t3,&S->x);                        //tmp3 = t3*X
//   fp_sub_nonmod_single(&tmp1_fp8.x0.x0,&tmp3_fp,&tmp2_fp);         // = t3*X - 2*t1

//   fp_mulmod_montgomery(&tmp1_fp8.x1.x1,&t3,&P->x);                 // = t3*Px
//   fp_mulmod_montgomery(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1,&tmp1_fp);  // = t3*Px*Z^2
//   fp_set_neg_montgomery(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1);       // = -t3*Px*Z^2
//   fp_mul_base_inv_single(&tmp1_fp8.x1.x1,&tmp1_fp8.x1.x1);  // = -t3*Px*Z^2*(2^-1)

//   fp8_mul_sparse_dbl_lazy_montgomery(f,&tmp1_fp8,f);        //Capable for further Karatsuba //update

//   fp_set(&S->x,&nextX);
//   fp_set(&S->y,&nextY);
//   fp_set(&S->z,&nextZ);
// }


void miller_opt_ate_jac(fp8_t *f,efp8_t *P,efp8_t *Q){
    static efp_t mapped_P;
    static efp2_t mapped_Q;
    static efp2_jacobian_t S;

    fp8_set_ui_ui(f,0);
    fp_set_ui(&f->x0.x0.x0,1);

    fp_set(&mapped_P.x,&P->x.x0.x0.x0);
    fp_set(&mapped_P.y,&P->y.x0.x0.x0);
    mapped_P.infinity = 0;

    efp8_to_efp2(&mapped_Q,Q);//twist
    efp8_to_Jacefp2(&S,Q);

    mp_bitcnt_t i;
    for(i=(mpz_sizeinbase(miller_loop_s,2)-2);i!=-1;i--){
      ff_lttp(f,&S,&mapped_P);
      if(mpz_tstbit(miller_loop_s,i)==1){
        ff_ltqp(f,&S,&mapped_Q,&mapped_P);
      }
    }
}

void miller_opt_ate_jac_2NAF(fp8_t *f,efp8_t *P,efp8_t *Q){
    static efp_t mapped_P;
    static efp2_t mapped_Q,mapped_Q_neg;
    static efp2_jacobian_t S;
    
    fp8_set_ui_ui(f,0);
    fp_set_ui(&f->x0.x0.x0,1);

    fp_set(&mapped_P.x,&P->x.x0.x0.x0);
    fp_set(&mapped_P.y,&P->y.x0.x0.x0);
    mapped_P.infinity = 0;

    efp8_to_efp2(&mapped_Q,Q);//twist
    efp8_to_Jacefp2(&S,Q);
    efp2_set_neg(&mapped_Q_neg,&mapped_Q);
    mp_bitcnt_t i;
    for(i=(miller_loop_v.size() -2);i!=-1;i--){
      switch(miller_loop_v[i]){
        case 0:
          ff_lttp(f,&S,&mapped_P);
          break;
        case 1:
          ff_lttp(f,&S,&mapped_P);
          ff_ltqp(f,&S,&mapped_Q,&mapped_P);
          break;
        case -1:
          ff_lttp(f,&S,&mapped_P);
          ff_ltqp(f,&S,&mapped_Q_neg,&mapped_P);
          break;
        default:
            break;
      }
    }
}
// void miller_opt_ate_proj_2NAF(fp8_t *f,efp8_t *P,efp8_t *Q){
//     static efp_t mapped_P;
//     static efp2_t mapped_Q,mapped_Q_neg;
//     static efp2_jacobian_t S;
    
//     fp8_set_ui_ui(f,0);
//     fp_set_ui(&f->x0.x0.x0,1);

//     fp_set(&mapped_P.x,&P->x.x0.x0.x0);
//     fp_set(&mapped_P.y,&P->y.x0.x0.x0);
//     mapped_P.infinity = 0;

//     efp8_to_efp2(&mapped_Q,Q);//twist
//     efp8_to_Jacefp2(&S,Q);
//     efp2_set_neg(&mapped_Q_neg,&mapped_Q);

//     mp_bitcnt_t i;
//     for(i=(miller_loop_v.size() -2);i!=-1;i--){
//       switch(miller_loop_v[i]){
//         case 0:
//           ff_lttp(f,&S,&mapped_P);
//           break;
//         case 1:
//           ff_lttp(f,&S,&mapped_P);
//           ff_ltqp(f,&S,&mapped_Q,&mapped_P);
//           break;
//         case -1:
//           ff_lttp(f,&S,&mapped_P);
//           ff_ltqp(f,&S,&mapped_Q_neg,&mapped_P);
//           break;
//         default:
//             break;
//       }
//     }
// }

void miller_opt_ate_jac_2NAF_lazy_montgomery(fp8_t *f,efp8_t *P,efp8_t *Q){
    static efp_t mapped_P;
    static efp2_t mapped_Q,mapped_Q_neg;
    static efp2_jacobian_t S;

    fp8_set_ui_ui(f,0);
    fp_set_ui(&f->x0.x0.x0,1);
    fp_to_montgomery(&f->x0.x0.x0, &f->x0.x0.x0);

    fp_to_montgomery(&mapped_P.x,&P->x.x0.x0.x0);
    fp_to_montgomery(&mapped_P.y,&P->y.x0.x0.x0);
    mapped_P.infinity = 0;

    efp8_to_efp2_montgomery(&mapped_Q,Q);//twist //not yet montogmery
    efp8_to_Jacefp2_montgomery(&S,Q);
    efp2_set_neg_montgomery(&mapped_Q_neg,&mapped_Q);//here
    mp_bitcnt_t i;
    for(i=(miller_loop_v.size() -2);i!=-1;i--){//-1
      switch(miller_loop_v[i]){
        case 0:
          ff_lttp_lazy_montgomery(f,&S,&mapped_P);
          break;
        case 1:
          ff_lttp_lazy_montgomery(f,&S,&mapped_P);
          ff_ltqp_lazy_montgomery(f,&S,&mapped_Q,&mapped_P);
          // fp8_println_montgomery("&f1", f);
          break;
        case -1:
          ff_lttp_lazy_montgomery(f,&S,&mapped_P);
          ff_ltqp_lazy_montgomery(f,&S,&mapped_Q_neg,&mapped_P);

          break;
        default:
            break;
      }
    }
}

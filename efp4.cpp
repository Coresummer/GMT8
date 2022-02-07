#include "efp4.h"
#include "fp4.h"

void efp4_init(efp4_t *P){
  fp4_init(&P->x);
  fp4_init(&P->y);
  P->infinity=0;
}

void efp4_projective_init(efp4_projective_t *P){
  fp4_init(&P->x);
  fp4_init(&P->y);
  fp4_init(&P->z);
  P->infinity=0;
}

void efp4_jacobian_init(efp4_jacobian_t *P){
  fp4_init(&P->x);
  fp4_init(&P->y);
  fp4_init(&P->z);
  P->infinity=0;
}
void efp4_printf(std::string str ,efp4_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp4_printf("",&P->x);
    printf(",");
    fp4_printf("",&P->y);
    printf(")");
  }else{
    printf("Infinity");
  }
}

void efp4_println(std::string str ,efp4_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp4_printf("",&P->x);
    printf(",");
    fp4_printf("",&P->y);
    printf(")\n");
  }else{
    printf("Infinity\n");
  }
}

void efp4_projective_printf(std::string str ,efp4_projective_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp4_printf("",&P->x);
    printf(",");
    fp4_printf("",&P->y);
    printf(",");
    fp4_printf("",&P->z);
    printf(")");
  }else{
    printf("Infinity");
  }
}

void efp4_jacobian_printf(std::string str ,efp4_jacobian_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp4_printf("",&P->x);
    printf(",");
    fp4_printf("",&P->y);
    printf(",");
    fp4_printf("",&P->z);
    printf(")");
  }else{
    printf("Infinity");
  }
}

// void efp4_printf_montgomery(std::string str ,efp4_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp4_printf_montgomery("",&P->x);
//     printf(",");
//     fp4_printf_montgomery("",&P->y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp4_jacobian_printf_montgomery(std::string str ,efp4_jacobian_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp4_printf_montgomery("",&P->x);
//     printf(",");
//     fp4_printf_montgomery("",&P->y);
//     printf(",");
//     fp4_printf_montgomery("",&P->z);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp4_projective_printf_montgomery(std::string str ,efp4_projective_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp4_printf_montgomery("",&P->x);
//     printf(",");
//     fp4_printf_montgomery("",&P->y);
//     printf(",");
//     fp4_printf_montgomery("",&P->z);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp4_projective_printf_affine(std::string str ,efp4_projective_t *P){
//   static efp4_t out;
//   efp4_projective_to_affine(&out,P);
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp4_printf("",&out.x);
//     printf(",");
//     fp4_printf("",&out.y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp4_projective_printf_affine_montgomery(std::string str ,efp4_projective_t *P){
//   static efp4_t out;
//   efp4_projective_to_affine_montgomery(&out,P);
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp4_printf_montgomery("",&out.x);
//     printf(",");
//     fp4_printf_montgomery("",&out.y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

void efp4_set(efp4_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp4_projective_set(efp4_projective_t *ANS,efp4_projective_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp4_jacobian_set(efp4_jacobian_t *ANS,efp4_jacobian_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp4_affine_to_projective(efp4_projective_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

void efp4_affine_to_jacobian(efp4_jacobian_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

void efp4_affine_to_projective_montgomery(efp4_projective_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set_mpn(&ANS->z,RmodP);
  ANS->infinity=A->infinity;
}

void efp4_affine_to_jacobian_montgomery(efp4_jacobian_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set(&ANS->y,&A->y);
  fp4_set_mpn(&ANS->z,RmodP);
  ANS->infinity=A->infinity;
}

void efp4_jacobian_to_affine(efp4_t *ANS,efp4_jacobian_t *A){
  static fp4_t Zi,Zt;
  //TODO:mul->mul_lazy
  fp4_inv(&Zi,&A->z);
  fp4_mul(&Zt,&Zi,&Zi);
  fp4_mul(&ANS->x,&A->x,&Zt);
  fp4_mul(&Zt,&Zt,&Zi);
  fp4_mul(&ANS->y,&A->y,&Zt);
  ANS->infinity=A->infinity;
}

void efp4_projective_to_affine(efp4_t *ANS,efp4_projective_t *A){
  static fp4_t Zi;
  //TODO:mul->mul_lazy
  fp4_inv(&Zi,&A->z);
  fp4_mul(&ANS->x,&A->x,&Zi);
  fp4_mul(&ANS->y,&A->y,&Zi);
  ANS->infinity=A->infinity;
}

// void efp4_jacobian_to_affine_montgomery(efp4_t *ANS,efp4_jacobian_t *A){
//   static fp4_t Zi,Zt;
//   fp4_inv_lazy_montgomery(&Zi,&A->z);
//   fp4_mul_lazy_montgomery(&Zt,&Zi,&Zi);
//   fp4_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
//   fp4_mul_lazy_montgomery(&Zt,&Zt,&Zi);
//   fp4_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
//   ANS->infinity=A->infinity;
// }

// void efp4_projective_to_affine_montgomery(efp4_t *ANS,efp4_projective_t *A){
//   static fp4_t Zi;
//   //TODO:mul->mul_lazy
//   fp4_inv_lazy_montgomery(&Zi,&A->z);
//   fp4_mul_lazy_montgomery(&ANS->x,&A->x,&Zi);
//   fp4_mul_lazy_montgomery(&ANS->y,&A->y,&Zi);
//   ANS->infinity=A->infinity;
// }

void efp4_mix(efp4_jacobian_t *ANS,efp4_jacobian_t *A,fp4_t *Zi){
  static fp4_t Zt;
  //TODO:mul->mul_lazy
  fp4_mul(&Zt,Zi,Zi);
  fp4_mul(&ANS->x,&A->x,&Zt);
  fp4_mul(&Zt,&Zt,Zi);
  fp4_mul(&ANS->y,&A->y,&Zt);
  fp4_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

// void efp4_mix_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *A,fp4_t *Zi){
//   static fp4_t Zt;
//   fp4_mul_lazy_montgomery(&Zt,Zi,Zi);
//   fp4_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
//   fp4_mul_lazy_montgomery(&Zt,&Zt,Zi);
//   fp4_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
//   fp4_set_mpn(&ANS->z,RmodP);
//   ANS->infinity=A->infinity;
// }

void efp4_set_ui(efp4_t *ANS,unsigned long int UI){
  fp4_set_ui(&ANS->x,UI);
  fp4_set_ui(&ANS->y,UI);
  ANS->infinity=0;
}

void efp4_to_montgomery(efp4_t *ANS,efp4_t *A){
  fp4_to_montgomery(&ANS->x,&A->x);
  fp4_to_montgomery(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp4_projective_to_montgomery(efp4_projective_t *ANS,efp4_projective_t *A){
  fp4_to_montgomery(&ANS->x,&A->x);
  fp4_to_montgomery(&ANS->y,&A->y);
  fp4_to_montgomery(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp4_mod_montgomery(efp4_t *ANS,efp4_t *A){
  fp4_mod_montgomery(&ANS->x,&A->x);
  fp4_mod_montgomery(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp4_projective_mod_montgomery(efp4_projective_t *ANS,efp4_projective_t *A){
  fp4_mod_montgomery(&ANS->x,&A->x);
  fp4_mod_montgomery(&ANS->y,&A->y);
  fp4_mod_montgomery(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp4_set_mpn(efp4_t *ANS,mp_limb_t *A){
  fp4_set_mpn(&ANS->x,A);
  fp4_set_mpn(&ANS->y,A);
  ANS->infinity=0;
}

void efp4_set_neg(efp4_t *ANS,efp4_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set_neg(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp4_jacobian_set_neg(efp4_jacobian_t *ANS,efp4_jacobian_t *A){
  fp4_set(&ANS->x,&A->x);
  fp4_set_neg(&ANS->y,&A->y);
  fp4_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

int efp4_cmp(efp4_t *A,efp4_t *B){
  if(A->infinity==1 && B->infinity==1)  return 0;
  else if(A->infinity != B->infinity) return 1;
  else if(fp4_cmp(&A->x,&B->x)==0 && fp4_cmp(&A->y,&B->y)==0) return 0;
  else  return 1;
}

void efp4_rational_point(efp4_t *P){
  fp4_t tmp_ax,tmp_y2;
  fp4_init(&tmp_ax);
  fp4_init(&tmp_y2);

  P->infinity=0;
  while(1){
    fp4_set_random(&P->x,state);
    //y^2 = x^3 + b
    fp4_sqr(&tmp_y2,&P->x);
    fp4_mul(&tmp_y2,&tmp_y2,&P->x);
    fp4_mul_mpn(&tmp_ax,&P->x,curve_a.x0); //ax
    fp4_add(&tmp_y2,&tmp_y2,&tmp_ax);//x^3 + ax
    if(fp4_legendre(&tmp_y2)==1){
      fp4_sqrt(&P->y,&tmp_y2);
      break;
    }
  }
}

void efp4_ecd(efp4_t *ANS,efp4_t *P){
  static efp4_t tmp1_efp4;
  static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
  if(P->infinity==1){
    ANS->infinity=1;
    return;
  }
  if(fp4_cmp_zero(&P->y)==0){
    ANS->infinity=1;
    return;
  }
  ANS->infinity=0;
  efp4_set(&tmp1_efp4,P);

  //tmp1_fp = 1/2yp
  fp4_add(&tmp1_fp4,&tmp1_efp4.y,&tmp1_efp4.y);
  fp4_inv(&tmp1_fp4,&tmp1_fp4);
  //tmp2_fp = 3x^2 +a
  fp4_sqr(&tmp2_fp4,&tmp1_efp4.x);
  fp4_add(&tmp3_fp4,&tmp2_fp4,&tmp2_fp4);
  fp4_add(&tmp2_fp4,&tmp2_fp4,&tmp3_fp4);
  fp2_add_mpn(&tmp2_fp4.x0,&tmp2_fp4.x0,curve_a.x0);
  //tmp3_fp = lambda
  fp4_mul(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);

  //ANS.x
  fp4_sqr(&tmp1_fp4,&tmp3_fp4);
  fp4_add(&tmp2_fp4,&tmp1_efp4.x,&tmp1_efp4.x);
  fp4_sub(&ANS->x,&tmp1_fp4,&tmp2_fp4);
  //ANS.y
  fp4_sub(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
  fp4_mul(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
  fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}

// void efp4_ecd_jacobian_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P){
//   static fp4_t s,m,T;

//   static fp4_t buf,tmp1;
//   static fp4_t tmpY2;
//   static efp4_jacobian_t Pt;
//   if(fp4_cmp_zero(&P->y)==0){
//     ANS->infinity=1;
//     return;
//   }

//   efp4_jacobian_set(&Pt,P);

//   //s
//   fp4_mul_lazy_montgomery(&tmpY2,&Pt.y,&Pt.y);
//   fp4_mul_lazy_montgomery(&tmp1,&tmpY2,&Pt.x);
//   fp4_add(&tmp1,&tmp1,&tmp1);
//   fp4_add(&s,&tmp1,&tmp1);

//   //m
//   fp4_add_nonmod_single(&tmp1,&Pt.x,&Pt.x);
//   fp4_add_nonmod_single(&tmp1,&tmp1,&Pt.x);
//   fp4_mul_lazy_montgomery(&m,&tmp1,&Pt.x);

//   //T
//   fp4_mul_lazy_montgomery(&T,&m,&m);
//   fp4_add(&tmp1,&s,&s);
//   fp4_sub(&T,&T,&tmp1);

//   //ANS->x
//   fp4_set(&ANS->x,&T);

//   //ANS->y
//   fp4_sub_nonmod_single(&tmp1,&s,&T);
//   fp4_mul_lazy_montgomery(&buf,&tmp1,&m);

//   fp4_mul_lazy_montgomery(&tmp1,&tmpY2,&tmpY2);
//   fp4_add(&tmp1,&tmp1,&tmp1);
//   fp4_add(&tmp1,&tmp1,&tmp1);
//   fp4_add(&tmp1,&tmp1,&tmp1);
//   fp4_sub(&ANS->y,&buf,&tmp1);

//   //ANS->z
//   fp4_add_nonmod_single(&tmp1,&Pt.y,&Pt.y);
//   fp4_mul_lazy_montgomery(&ANS->z,&tmp1,&Pt.z);
// }

void efp4_eca(efp4_t *ANS,efp4_t *P1,efp4_t *P2){
  static efp4_t tmp1_efp4,tmp2_efp4;
  static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
  if(P1->infinity==1){
    efp4_set(ANS,P2);
    return;
  }
  else if(P2->infinity==1){
    efp4_set(ANS,P1);
    return;
  }
  else if(fp4_cmp(&P1->x,&P2->x)==0){
    if(fp4_cmp(&P1->y,&P2->y)!=0){
      ANS->infinity=1;
      return;
    }else{
      efp4_ecd(ANS,P1);
      return;
    }
  }
  ANS->infinity=0;
  efp4_set(&tmp1_efp4,P1);
  efp4_set(&tmp2_efp4,P2);

  //tmp3_fp = lambda
  fp4_sub(&tmp1_fp4,&tmp2_efp4.x,&tmp1_efp4.x);
  fp4_inv(&tmp1_fp4,&tmp1_fp4);
  fp4_sub(&tmp2_fp4,&tmp2_efp4.y,&tmp1_efp4.y);
  fp4_mul(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);

  //ANS.x
  fp4_sqr(&tmp1_fp4,&tmp3_fp4);
  fp4_sub(&tmp2_fp4,&tmp1_fp4,&tmp1_efp4.x);
  fp4_sub(&ANS->x,&tmp2_fp4,&tmp2_efp4.x);
  //ANS.y
  fp4_sub(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
  fp4_mul(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
  fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}

// void efp4_eca_jacobian_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P1,efp4_jacobian_t *P2){
//   static efp4_jacobian_t Pt1,Pt2;
//   static fp4_t U1,U2,S1,S2,H,r;

//   static fp4_t buf,tmp1,tmp2;
//   static fp4_t tmpZ1,tmpZ2,tmpH2,tmpH3,tmpU1H2;

//   if(P1->infinity==1){
//     efp4_jacobian_set(ANS,P2);
//     return;
//   }else if(P2->infinity==1){
//     efp4_jacobian_set(ANS,P1);
//     return;
//   }else if(fp4_cmp(&P1->x,&P2->x)==0){
//     if(fp4_cmp(&P1->y,&P2->y)!=0){
//       ANS->infinity=1;
//       return;
//     }else{
//       efp4_ecd_jacobian_lazy_montgomery(ANS,P1);
//       return;
//     }
//   }

//   efp4_jacobian_set(&Pt1,P1);
//   efp4_jacobian_set(&Pt2,P2);

//   //U1
//   fp4_mul_lazy_montgomery(&tmpZ2,&Pt2.z,&Pt2.z);
//   fp4_mul_lazy_montgomery(&U1,&tmpZ2,&Pt1.x);
//   //fp4_printf("U1=",&U1);printf("\n");

//   //U2
//   fp4_mul_lazy_montgomery(&tmpZ1,&Pt1.z,&Pt1.z);
//   fp4_mul_lazy_montgomery(&U2,&tmpZ1,&Pt2.x);
//   //fp4_printf("U2=",&U2);printf("\n");

//   //S1
//   fp4_mul_lazy_montgomery(&tmp1,&tmpZ2,&Pt2.z);
//   fp4_mul_lazy_montgomery(&S1,&tmp1,&Pt1.y);
//   //fp4_printf("S1=",&S1);printf("\n");

//   //S2
//   fp4_mul_lazy_montgomery(&tmp1,&tmpZ1,&Pt1.z);
//   fp4_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);
//   //fp4_printf("S2=",&S2);printf("\n");

//   //H
//   //fp4_printf("U1=",&U1);printf("\n");
//   fp4_sub(&H,&U2,&U1);
//   //fp4_printf("H=",&H);printf("\n");

//   //r
//   fp4_sub(&r,&S2,&S1);
//   //fp4_printf("r=",&r);printf("\n");

//   //ANS->x
//   fp4_mul_lazy_montgomery(&tmp1,&r,&r);

//   fp4_mul_lazy_montgomery(&tmpH2,&H,&H);
//   fp4_mul_lazy_montgomery(&tmpH3,&tmpH2,&H);
//   fp4_sub(&tmp2,&tmp1,&tmpH3);

//   fp4_mul_lazy_montgomery(&tmpU1H2,&tmpH2,&U1);
//   fp4_add(&tmp1,&tmpU1H2,&tmpU1H2);
//   fp4_sub(&ANS->x,&tmp2,&tmp1);

//   //ANS->y
//   fp4_sub_nonmod_single(&tmp1,&tmpU1H2,&ANS->x);
//   fp4_mul_lazy_montgomery(&tmp1,&tmp1,&r);

//   fp4_mul_lazy_montgomery(&tmp2,&tmpH3,&S1);
//   fp4_sub(&ANS->y,&tmp1,&tmp2);

//   //ANS->z
//   fp4_mul_lazy_montgomery(&tmp1,&Pt1.z,&Pt2.z);
//   fp4_mul_lazy_montgomery(&ANS->z,&tmp1,&H);
//   // //getchar();
// }
// void efp4_eca_mixture_lazy_montgomery(efp4_jacobian_t *ANS,efp4_jacobian_t *P1,efp4_jacobian_t *P2){
//   static efp4_jacobian_t Pt1,Pt2;
//   static fp4_t Z1Z1,HH,I,J,V;
//   static fp4_t U1,U2,S1,S2,H,r;
//   static fp4_t buf,tmp1,tmp2;

//   if(P1->infinity==1){
//     efp4_jacobian_set(ANS,P2);
//     return;
//   }else if(P2->infinity==1){
//     efp4_jacobian_set(ANS,P1);
//     return;
//   }else if(fp4_cmp(&P1->x,&P2->x)==0){
//     if(fp4_cmp(&P1->y,&P2->y)!=0){
//       ANS->infinity=1;
//       return;
//     }else{
//       efp4_ecd_jacobian_lazy_montgomery(ANS,P1);
//       return;
//     }
//   }

//   efp4_jacobian_set(&Pt1,P1);
//   efp4_jacobian_set(&Pt2,P2);

//   //Z1Z1
//   fp4_mul_lazy_montgomery(&Z1Z1,&Pt1.z,&Pt1.z);

//   //U2
//   fp4_mul_lazy_montgomery(&U2,&Pt2.x,&Z1Z1);

//   //S2
//   fp4_mul_lazy_montgomery(&tmp1,&Z1Z1,&Pt1.z);
//   fp4_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);

//   //H
//   fp4_sub(&H,&U2,&Pt1.x);

//   //HH
//   fp4_mul_lazy_montgomery(&HH,&H,&H);

//   //I
//   fp4_add(&I,&HH,&HH);
//   fp4_add(&I,&I,&I);

//   //J
//   fp4_mul_lazy_montgomery(&J,&HH,&H);

//   //r
//   fp4_sub(&r,&S2,&Pt1.y);

//   //V
//   fp4_mul_lazy_montgomery(&V,&Pt1.x,&HH);

//   //X3
//   fp4_mul_lazy_montgomery(&tmp1,&r,&r);
//   fp4_add(&tmp2,&V,&V);
//   fp4_sub(&buf,&tmp1,&J);
//   fp4_sub(&ANS->x,&buf,&tmp2);

//   //Y3
//   fp4_sub_nonmod_single(&tmp1,&V,&ANS->x);
//   fp4_mul_lazy_montgomery(&tmp2,&tmp1,&r);
//   fp4_mul_lazy_montgomery(&tmp1,&Pt1.y,&J);
//   fp4_sub(&ANS->y,&tmp2,&tmp1);


//   //ANS->z
//   fp4_mul_lazy_montgomery(&ANS->z,&Pt1.z,&H);

// }

void efp4_scm(efp4_t *ANS,efp4_t *P,mpz_t scalar){
  if(mpz_cmp_ui(scalar,0)==0){
    ANS->infinity=1;
    return;
  }else if(mpz_cmp_ui(scalar,1)==0){
    efp4_set(ANS,P);
    return;
  }

  efp4_t Tmp_P,Next_P;
  efp4_init(&Tmp_P);
  efp4_set(&Tmp_P,P);
  efp4_init(&Next_P);
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);

  efp4_set(&Next_P,&Tmp_P);
  for(i=1;i<length;i++){
    efp4_ecd(&Next_P,&Next_P);
    if(binary[i]=='1')  efp4_eca(&Next_P,&Next_P,&Tmp_P);
  }
  efp4_set(ANS,&Next_P);
}

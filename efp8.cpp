#include "efp8.h"
#include "fp8.h"

void efp8_init(efp8_t *P){
  fp8_init(&P->x);
  fp8_init(&P->y);
  P->infinity=0;
}

void efp8_projective_init(efp8_projective_t *P){
  fp8_init(&P->x);
  fp8_init(&P->y);
  fp8_init(&P->z);
  P->infinity=0;
}

void efp8_jacobian_init(efp8_jacobian_t *P){
  fp8_init(&P->x);
  fp8_init(&P->y);
  fp8_init(&P->z);
  P->infinity=0;
}
void efp8_printf(std::string str ,efp8_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp8_printf("",&P->x);
    printf(",");
    fp8_printf("",&P->y);
    printf(")");
  }else{
    printf("Infinity");
  }
}

void efp8_println(std::string str ,efp8_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp8_printf("",&P->x);
    printf(",");
    fp8_printf("",&P->y);
    printf(")\n");
  }else{
    printf("Infinity\n");
  }
}

void efp8_projective_printf(std::string str ,efp8_projective_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp8_printf("",&P->x);
    printf(",");
    fp8_printf("",&P->y);
    printf(",");
    fp8_printf("",&P->z);
    printf(")");
  }else{
    printf("Infinity");
  }
}

void efp8_jacobian_printf(std::string str ,efp8_jacobian_t *P){
  printf("%s",str.c_str());
  if(P->infinity==0){
    printf("(");
    fp8_printf("",&P->x);
    printf(",");
    fp8_printf("",&P->y);
    printf(",");
    fp8_printf("",&P->z);
    printf(")");
  }else{
    printf("Infinity");
  }
}

// void efp8_printf_montgomery(std::string str ,efp8_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp8_printf_montgomery("",&P->x);
//     printf(",");
//     fp8_printf_montgomery("",&P->y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp8_jacobian_printf_montgomery(std::string str ,efp8_jacobian_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp8_printf_montgomery("",&P->x);
//     printf(",");
//     fp8_printf_montgomery("",&P->y);
//     printf(",");
//     fp8_printf_montgomery("",&P->z);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp8_projective_printf_montgomery(std::string str ,efp8_projective_t *P){
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp8_printf_montgomery("",&P->x);
//     printf(",");
//     fp8_printf_montgomery("",&P->y);
//     printf(",");
//     fp8_printf_montgomery("",&P->z);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp8_projective_printf_affine(std::string str ,efp8_projective_t *P){
//   static efp8_t out;
//   efp8_projective_to_affine(&out,P);
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp8_printf("",&out.x);
//     printf(",");
//     fp8_printf("",&out.y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

// void efp8_projective_printf_affine_montgomery(std::string str ,efp8_projective_t *P){
//   static efp8_t out;
//   efp8_projective_to_affine_montgomery(&out,P);
//   printf("%s",str.c_str());
//   if(P->infinity==0){
//     printf("(");
//     fp8_printf_montgomery("",&out.x);
//     printf(",");
//     fp8_printf_montgomery("",&out.y);
//     printf(")");
//   }else{
//     printf("Infinity");
//   }
// }

void efp8_set(efp8_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp8_projective_set(efp8_projective_t *ANS,efp8_projective_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp8_jacobian_set(efp8_jacobian_t *ANS,efp8_jacobian_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp8_affine_to_projective(efp8_projective_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

void efp8_affine_to_jacobian(efp8_jacobian_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

void efp8_affine_to_projective_montgomery(efp8_projective_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set_mpn(&ANS->z,RmodP);
  ANS->infinity=A->infinity;
}

void efp8_affine_to_jacobian_montgomery(efp8_jacobian_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set(&ANS->y,&A->y);
  fp8_set_mpn(&ANS->z,RmodP);
  ANS->infinity=A->infinity;
}

void efp8_jacobian_to_affine(efp8_t *ANS,efp8_jacobian_t *A){
  static fp8_t Zi,Zt;
  //TODO:mul->mul_lazy
  fp8_inv(&Zi,&A->z);
  fp8_mul(&Zt,&Zi,&Zi);
  fp8_mul(&ANS->x,&A->x,&Zt);
  fp8_mul(&Zt,&Zt,&Zi);
  fp8_mul(&ANS->y,&A->y,&Zt);
  ANS->infinity=A->infinity;
}

void efp8_projective_to_affine(efp8_t *ANS,efp8_projective_t *A){
  static fp8_t Zi;
  //TODO:mul->mul_lazy
  fp8_inv(&Zi,&A->z);
  fp8_mul(&ANS->x,&A->x,&Zi);
  fp8_mul(&ANS->y,&A->y,&Zi);
  ANS->infinity=A->infinity;
}

// void efp8_jacobian_to_affine_montgomery(efp8_t *ANS,efp8_jacobian_t *A){
//   static fp8_t Zi,Zt;
//   fp8_inv_lazy_montgomery(&Zi,&A->z);
//   fp8_mul_lazy_montgomery(&Zt,&Zi,&Zi);
//   fp8_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
//   fp8_mul_lazy_montgomery(&Zt,&Zt,&Zi);
//   fp8_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
//   ANS->infinity=A->infinity;
// }

// void efp8_projective_to_affine_montgomery(efp8_t *ANS,efp8_projective_t *A){
//   static fp8_t Zi;
//   //TODO:mul->mul_lazy
//   fp8_inv_lazy_montgomery(&Zi,&A->z);
//   fp8_mul_lazy_montgomery(&ANS->x,&A->x,&Zi);
//   fp8_mul_lazy_montgomery(&ANS->y,&A->y,&Zi);
//   ANS->infinity=A->infinity;
// }

void efp8_mix(efp8_jacobian_t *ANS,efp8_jacobian_t *A,fp8_t *Zi){
  static fp8_t Zt;
  //TODO:mul->mul_lazy
  fp8_mul(&Zt,Zi,Zi);
  fp8_mul(&ANS->x,&A->x,&Zt);
  fp8_mul(&Zt,&Zt,Zi);
  fp8_mul(&ANS->y,&A->y,&Zt);
  fp8_set_ui(&ANS->z,1);
  ANS->infinity=A->infinity;
}

// void efp8_mix_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *A,fp8_t *Zi){
//   static fp8_t Zt;
//   fp8_mul_lazy_montgomery(&Zt,Zi,Zi);
//   fp8_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
//   fp8_mul_lazy_montgomery(&Zt,&Zt,Zi);
//   fp8_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
//   fp8_set_mpn(&ANS->z,RmodP);
//   ANS->infinity=A->infinity;
// }

void efp8_set_ui(efp8_t *ANS,unsigned long int UI){
  fp8_set_ui(&ANS->x,UI);
  fp8_set_ui(&ANS->y,UI);
  ANS->infinity=0;
}

void efp8_to_montgomery(efp8_t *ANS,efp8_t *A){
  fp8_to_montgomery(&ANS->x,&A->x);
  fp8_to_montgomery(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp8_projective_to_montgomery(efp8_projective_t *ANS,efp8_projective_t *A){
  fp8_to_montgomery(&ANS->x,&A->x);
  fp8_to_montgomery(&ANS->y,&A->y);
  fp8_to_montgomery(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp8_mod_montgomery(efp8_t *ANS,efp8_t *A){
  fp8_mod_montgomery(&ANS->x,&A->x);
  fp8_mod_montgomery(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp8_projective_mod_montgomery(efp8_projective_t *ANS,efp8_projective_t *A){
  fp8_mod_montgomery(&ANS->x,&A->x);
  fp8_mod_montgomery(&ANS->y,&A->y);
  fp8_mod_montgomery(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

void efp8_set_mpn(efp8_t *ANS,mp_limb_t *A){
  fp8_set_mpn(&ANS->x,A);
  fp8_set_mpn(&ANS->y,A);
  ANS->infinity=0;
}

void efp8_set_neg(efp8_t *ANS,efp8_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set_neg(&ANS->y,&A->y);
  ANS->infinity=A->infinity;
}

void efp8_jacobian_set_neg(efp8_jacobian_t *ANS,efp8_jacobian_t *A){
  fp8_set(&ANS->x,&A->x);
  fp8_set_neg(&ANS->y,&A->y);
  fp8_set(&ANS->z,&A->z);
  ANS->infinity=A->infinity;
}

int efp8_cmp(efp8_t *A,efp8_t *B){
  if(A->infinity==1 && B->infinity==1)  return 0;
  else if(A->infinity != B->infinity) return 1;
  else if(fp8_cmp(&A->x,&B->x)==0 && fp8_cmp(&A->y,&B->y)==0) return 0;
  else  return 1;
}

void efp8_rational_point(efp8_t *P){
  fp8_t tmp_ax,tmp_y2;
  fp8_init(&tmp_ax);
  fp8_init(&tmp_y2);

  P->infinity=0;
  while(1){
    fp8_set_random(&P->x,state);
    //y^2 = x^3 + b
    fp8_sqr(&tmp_y2,&P->x);
    fp8_mul(&tmp_y2,&tmp_y2,&P->x);
    fp8_mul_mpn(&tmp_ax,&P->x,curve_a.x0); //ax
    fp8_add(&tmp_y2,&tmp_y2,&tmp_ax);//x^3 + ax
    if(fp8_legendre(&tmp_y2)==1){
      fp8_sqrt(&P->y,&tmp_y2);
      break;
    }
  }
}


void generate_g1(efp8_t *P){
  //g1は定義体が素体である楕円曲線上の位数rの有理点群
  efp_t tmp_P;
  efp_init(&tmp_P);
  mpz_t expo;
  mpz_init(expo);

  efp_rational_point(&tmp_P);
  efp8_set_ui(P,0);
  mpz_tdiv_q(expo,efp_total,order_z);//(#efp)/r
  efp_scm(&tmp_P,&tmp_P,expo);
  fp_set(&P->x.x0.x0.x0,&tmp_P.x);
  fp_set(&P->y.x0.x0.x0,&tmp_P.y);
  P->infinity=tmp_P.infinity;

  mpz_clear(expo);
}

void generate_g2(efp8_t *Q){
  //g2は定義体が6次拡大体である楕円曲線上の位数rの有理点群
  efp8_t tmp_Q,frobenius_Q;
  efp8_init(&tmp_Q);
  efp8_init(&frobenius_Q);

  mpz_t expo;
  mpz_init(expo);
  efp8_rational_point(&tmp_Q);
  mpz_pow_ui(expo,order_z,2);
  mpz_divexact(expo,efp8_total,expo);//(#efp8)/(r^2)
  efp8_scm(&tmp_Q,&tmp_Q,expo);//tmp_Qは位数rの有理点となった
  // efp8_set(Q,&tmp_Q); //Comfirmed on r

  //(Φ-1)tmp_Q = Q となる有理点Qが特別な花びら上の有理点となる
  efp8_frobenius_map_p1(&frobenius_Q,&tmp_Q); //(Q^p) frobenius mapping confirmed
  // efp8_println("in gen G2", &frobenius_Q);
  efp8_set_neg(&tmp_Q,&tmp_Q);//(-Q)
  efp8_eca(Q,&tmp_Q,&frobenius_Q);//Q^p - Q
  mpz_clear(expo);
}

void efp8_ecd(efp8_t *ANS,efp8_t *P){
  static efp8_t tmp1_efp8;
  static fp8_t tmp1_fp8,tmp2_fp8,tmp3_fp8;
  if(P->infinity==1){
    ANS->infinity=1;
    return;
  }
  if(fp8_cmp_zero(&P->y)==0){
    ANS->infinity=1;
    return;
  }
  ANS->infinity=0;
  efp8_set(&tmp1_efp8,P);

  //tmp1_fp = 1/2yp
  fp8_add(&tmp1_fp8,&tmp1_efp8.y,&tmp1_efp8.y);
  fp8_inv(&tmp1_fp8,&tmp1_fp8);
  //tmp2_fp = 3x^2 +a
  fp8_sqr(&tmp2_fp8,&tmp1_efp8.x);
  fp8_add(&tmp3_fp8,&tmp2_fp8,&tmp2_fp8);
  fp8_add(&tmp2_fp8,&tmp2_fp8,&tmp3_fp8);
  fp_add(&tmp2_fp8.x0.x0.x0,&tmp2_fp8.x0.x0.x0,&curve_a);
  //tmp3_fp = lambda
  fp8_mul(&tmp3_fp8,&tmp1_fp8,&tmp2_fp8);

  //ANS.x
  fp8_sqr(&tmp1_fp8,&tmp3_fp8);
  fp8_add(&tmp2_fp8,&tmp1_efp8.x,&tmp1_efp8.x);
  fp8_sub(&ANS->x,&tmp1_fp8,&tmp2_fp8);
  //ANS.y
  fp8_sub(&tmp1_fp8,&tmp1_efp8.x,&ANS->x);
  fp8_mul(&tmp2_fp8,&tmp3_fp8,&tmp1_fp8);
  fp8_sub(&ANS->y,&tmp2_fp8,&tmp1_efp8.y);
}


void efp8_ecd_dash(efp8_t *ANS,efp8_t *P){
  static efp8_t tmp1_efp8;
  static fp8_t tmp1_fp8,tmp2_fp8,tmp3_fp8,tmpa_fp8;
  if(P->infinity==1){
    ANS->infinity=1;
    return;
  }
  if(fp8_cmp_zero(&P->y)==0){
    ANS->infinity=1;
    return;
  }
  ANS->infinity=0;
  efp8_set(&tmp1_efp8,P);
  fp8_set_mpn(&tmpa_fp8,curve_a.x0);
  fp2_mul_base(&tmpa_fp8.x0.x0, &tmpa_fp8.x0.x0);

  //tmp1_fp = 1/2yp
  fp8_add(&tmp1_fp8,&tmp1_efp8.y,&tmp1_efp8.y);
  fp8_inv(&tmp1_fp8,&tmp1_fp8);
  //tmp2_fp = 3x^2 +a
  fp8_sqr(&tmp2_fp8,&tmp1_efp8.x);
  fp8_add(&tmp3_fp8,&tmp2_fp8,&tmp2_fp8);
  fp8_add(&tmp2_fp8,&tmp2_fp8,&tmp3_fp8);
  fp8_add(&tmp2_fp8,&tmp2_fp8,&tmpa_fp8);
  //tmp3_fp = lambda
  fp8_mul(&tmp3_fp8,&tmp1_fp8,&tmp2_fp8);

  //ANS.x
  fp8_sqr(&tmp1_fp8,&tmp3_fp8);
  fp8_add(&tmp2_fp8,&tmp1_efp8.x,&tmp1_efp8.x);
  fp8_sub(&ANS->x,&tmp1_fp8,&tmp2_fp8);
  //ANS.y
  fp8_sub(&tmp1_fp8,&tmp1_efp8.x,&ANS->x);
  fp8_mul(&tmp2_fp8,&tmp3_fp8,&tmp1_fp8);
  fp8_sub(&ANS->y,&tmp2_fp8,&tmp1_efp8.y);
}

// void efp8_ecd_jacobian_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P){
//   static fp8_t s,m,T;

//   static fp8_t buf,tmp1;
//   static fp8_t tmpY2;
//   static efp8_jacobian_t Pt;
//   if(fp8_cmp_zero(&P->y)==0){
//     ANS->infinity=1;
//     return;
//   }

//   efp8_jacobian_set(&Pt,P);

//   //s
//   fp8_mul_lazy_montgomery(&tmpY2,&Pt.y,&Pt.y);
//   fp8_mul_lazy_montgomery(&tmp1,&tmpY2,&Pt.x);
//   fp8_add(&tmp1,&tmp1,&tmp1);
//   fp8_add(&s,&tmp1,&tmp1);

//   //m
//   fp8_add_nonmod_single(&tmp1,&Pt.x,&Pt.x);
//   fp8_add_nonmod_single(&tmp1,&tmp1,&Pt.x);
//   fp8_mul_lazy_montgomery(&m,&tmp1,&Pt.x);

//   //T
//   fp8_mul_lazy_montgomery(&T,&m,&m);
//   fp8_add(&tmp1,&s,&s);
//   fp8_sub(&T,&T,&tmp1);

//   //ANS->x
//   fp8_set(&ANS->x,&T);

//   //ANS->y
//   fp8_sub_nonmod_single(&tmp1,&s,&T);
//   fp8_mul_lazy_montgomery(&buf,&tmp1,&m);

//   fp8_mul_lazy_montgomery(&tmp1,&tmpY2,&tmpY2);
//   fp8_add(&tmp1,&tmp1,&tmp1);
//   fp8_add(&tmp1,&tmp1,&tmp1);
//   fp8_add(&tmp1,&tmp1,&tmp1);
//   fp8_sub(&ANS->y,&buf,&tmp1);

//   //ANS->z
//   fp8_add_nonmod_single(&tmp1,&Pt.y,&Pt.y);
//   fp8_mul_lazy_montgomery(&ANS->z,&tmp1,&Pt.z);
// }

void efp8_eca(efp8_t *ANS,efp8_t *P1,efp8_t *P2){
  static efp8_t tmp1_efp8,tmp2_efp8;
  static fp8_t tmp1_fp8,tmp2_fp8,tmp3_fp8;
  if(P1->infinity==1){
    efp8_set(ANS,P2);
    return;
  }
  else if(P2->infinity==1){
    efp8_set(ANS,P1);
    return;
  }
  else if(fp8_cmp(&P1->x,&P2->x)==0){
    if(fp8_cmp(&P1->y,&P2->y)!=0){
      ANS->infinity=1;
      return;
    }else{
      efp8_ecd(ANS,P1);
      return;
    }
  }
  ANS->infinity=0;
  efp8_set(&tmp1_efp8,P1);
  efp8_set(&tmp2_efp8,P2);

  //tmp3_fp = lambda
  fp8_sub(&tmp1_fp8,&tmp2_efp8.x,&tmp1_efp8.x);
  fp8_inv(&tmp1_fp8,&tmp1_fp8);
  fp8_sub(&tmp2_fp8,&tmp2_efp8.y,&tmp1_efp8.y);
  fp8_mul(&tmp3_fp8,&tmp1_fp8,&tmp2_fp8);

  //ANS.x
  fp8_sqr(&tmp1_fp8,&tmp3_fp8);
  fp8_sub(&tmp2_fp8,&tmp1_fp8,&tmp1_efp8.x);
  fp8_sub(&ANS->x,&tmp2_fp8,&tmp2_efp8.x);
  //ANS.y
  fp8_sub(&tmp1_fp8,&tmp1_efp8.x,&ANS->x);
  fp8_mul(&tmp2_fp8,&tmp3_fp8,&tmp1_fp8);
  fp8_sub(&ANS->y,&tmp2_fp8,&tmp1_efp8.y);
}

// void efp8_eca_jacobian_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P1,efp8_jacobian_t *P2){
//   static efp8_jacobian_t Pt1,Pt2;
//   static fp8_t U1,U2,S1,S2,H,r;

//   static fp8_t buf,tmp1,tmp2;
//   static fp8_t tmpZ1,tmpZ2,tmpH2,tmpH3,tmpU1H2;

//   if(P1->infinity==1){
//     efp8_jacobian_set(ANS,P2);
//     return;
//   }else if(P2->infinity==1){
//     efp8_jacobian_set(ANS,P1);
//     return;
//   }else if(fp8_cmp(&P1->x,&P2->x)==0){
//     if(fp8_cmp(&P1->y,&P2->y)!=0){
//       ANS->infinity=1;
//       return;
//     }else{
//       efp8_ecd_jacobian_lazy_montgomery(ANS,P1);
//       return;
//     }
//   }

//   efp8_jacobian_set(&Pt1,P1);
//   efp8_jacobian_set(&Pt2,P2);

//   //U1
//   fp8_mul_lazy_montgomery(&tmpZ2,&Pt2.z,&Pt2.z);
//   fp8_mul_lazy_montgomery(&U1,&tmpZ2,&Pt1.x);
//   //fp8_printf("U1=",&U1);printf("\n");

//   //U2
//   fp8_mul_lazy_montgomery(&tmpZ1,&Pt1.z,&Pt1.z);
//   fp8_mul_lazy_montgomery(&U2,&tmpZ1,&Pt2.x);
//   //fp8_printf("U2=",&U2);printf("\n");

//   //S1
//   fp8_mul_lazy_montgomery(&tmp1,&tmpZ2,&Pt2.z);
//   fp8_mul_lazy_montgomery(&S1,&tmp1,&Pt1.y);
//   //fp8_printf("S1=",&S1);printf("\n");

//   //S2
//   fp8_mul_lazy_montgomery(&tmp1,&tmpZ1,&Pt1.z);
//   fp8_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);
//   //fp8_printf("S2=",&S2);printf("\n");

//   //H
//   //fp8_printf("U1=",&U1);printf("\n");
//   fp8_sub(&H,&U2,&U1);
//   //fp8_printf("H=",&H);printf("\n");

//   //r
//   fp8_sub(&r,&S2,&S1);
//   //fp8_printf("r=",&r);printf("\n");

//   //ANS->x
//   fp8_mul_lazy_montgomery(&tmp1,&r,&r);

//   fp8_mul_lazy_montgomery(&tmpH2,&H,&H);
//   fp8_mul_lazy_montgomery(&tmpH3,&tmpH2,&H);
//   fp8_sub(&tmp2,&tmp1,&tmpH3);

//   fp8_mul_lazy_montgomery(&tmpU1H2,&tmpH2,&U1);
//   fp8_add(&tmp1,&tmpU1H2,&tmpU1H2);
//   fp8_sub(&ANS->x,&tmp2,&tmp1);

//   //ANS->y
//   fp8_sub_nonmod_single(&tmp1,&tmpU1H2,&ANS->x);
//   fp8_mul_lazy_montgomery(&tmp1,&tmp1,&r);

//   fp8_mul_lazy_montgomery(&tmp2,&tmpH3,&S1);
//   fp8_sub(&ANS->y,&tmp1,&tmp2);

//   //ANS->z
//   fp8_mul_lazy_montgomery(&tmp1,&Pt1.z,&Pt2.z);
//   fp8_mul_lazy_montgomery(&ANS->z,&tmp1,&H);
//   // //getchar();
// }
// void efp8_eca_mixture_lazy_montgomery(efp8_jacobian_t *ANS,efp8_jacobian_t *P1,efp8_jacobian_t *P2){
//   static efp8_jacobian_t Pt1,Pt2;
//   static fp8_t Z1Z1,HH,I,J,V;
//   static fp8_t U1,U2,S1,S2,H,r;
//   static fp8_t buf,tmp1,tmp2;

//   if(P1->infinity==1){
//     efp8_jacobian_set(ANS,P2);
//     return;
//   }else if(P2->infinity==1){
//     efp8_jacobian_set(ANS,P1);
//     return;
//   }else if(fp8_cmp(&P1->x,&P2->x)==0){
//     if(fp8_cmp(&P1->y,&P2->y)!=0){
//       ANS->infinity=1;
//       return;
//     }else{
//       efp8_ecd_jacobian_lazy_montgomery(ANS,P1);
//       return;
//     }
//   }

//   efp8_jacobian_set(&Pt1,P1);
//   efp8_jacobian_set(&Pt2,P2);

//   //Z1Z1
//   fp8_mul_lazy_montgomery(&Z1Z1,&Pt1.z,&Pt1.z);

//   //U2
//   fp8_mul_lazy_montgomery(&U2,&Pt2.x,&Z1Z1);

//   //S2
//   fp8_mul_lazy_montgomery(&tmp1,&Z1Z1,&Pt1.z);
//   fp8_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);

//   //H
//   fp8_sub(&H,&U2,&Pt1.x);

//   //HH
//   fp8_mul_lazy_montgomery(&HH,&H,&H);

//   //I
//   fp8_add(&I,&HH,&HH);
//   fp8_add(&I,&I,&I);

//   //J
//   fp8_mul_lazy_montgomery(&J,&HH,&H);

//   //r
//   fp8_sub(&r,&S2,&Pt1.y);

//   //V
//   fp8_mul_lazy_montgomery(&V,&Pt1.x,&HH);

//   //X3
//   fp8_mul_lazy_montgomery(&tmp1,&r,&r);
//   fp8_add(&tmp2,&V,&V);
//   fp8_sub(&buf,&tmp1,&J);
//   fp8_sub(&ANS->x,&buf,&tmp2);

//   //Y3
//   fp8_sub_nonmod_single(&tmp1,&V,&ANS->x);
//   fp8_mul_lazy_montgomery(&tmp2,&tmp1,&r);
//   fp8_mul_lazy_montgomery(&tmp1,&Pt1.y,&J);
//   fp8_sub(&ANS->y,&tmp2,&tmp1);


//   //ANS->z
//   fp8_mul_lazy_montgomery(&ANS->z,&Pt1.z,&H);

// }

void efp8_scm(efp8_t *ANS,efp8_t *P,mpz_t scalar){
  if(mpz_cmp_ui(scalar,0)==0){
    ANS->infinity=1;
    return;
  }else if(mpz_cmp_ui(scalar,1)==0){
    efp8_set(ANS,P);
    return;
  }

  efp8_t Tmp_P,Next_P;
  efp8_init(&Tmp_P);
  efp8_set(&Tmp_P,P);
  efp8_init(&Next_P);
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);

  efp8_set(&Next_P,&Tmp_P);
  for(i=1;i<length;i++){
    efp8_ecd(&Next_P,&Next_P);
    if(binary[i]=='1')  efp8_eca(&Next_P,&Next_P,&Tmp_P);
  }
  efp8_set(ANS,&Next_P);
}

void efp8_scm_dash(efp8_t *ANS,efp8_t *P,mpz_t scalar){
  if(mpz_cmp_ui(scalar,0)==0){
    ANS->infinity=1;
    return;
  }else if(mpz_cmp_ui(scalar,1)==0){
    efp8_set(ANS,P);
    return;
  }

  efp8_t Tmp_P,Next_P;
  efp8_init(&Tmp_P);
  efp8_set(&Tmp_P,P);
  efp8_init(&Next_P);
  int i,length;
  length=(int)mpz_sizeinbase(scalar,2);
  char binary[length+1];
  mpz_get_str(binary,2,scalar);

  efp8_set(&Next_P,&Tmp_P);
  for(i=1;i<length;i++){
    efp8_ecd_dash(&Next_P,&Next_P);
    if(binary[i]=='1')  efp8_eca(&Next_P,&Next_P,&Tmp_P);
  }
  efp8_set(ANS,&Next_P);
}

void efp8_frobenius_map_p1(efp8_t *ANS,efp8_t *A){
  if(A->infinity==1){
    ANS->infinity=1;
    return;
  }
  ANS->infinity=0;
  fp8_frobenius_map_p1(&ANS->x,&A->x);
  fp8_frobenius_map_p1(&ANS->y,&A->y);
}

void efp8_checkOnCurve(efp8_t *A){
  static fp8_t tmp_left_fp8, tmp_right_fp8,tmp_ax;
  fp8_sqr(&tmp_left_fp8,&A->y);

  fp8_sqr(&tmp_right_fp8,&A->x);
  fp8_mul(&tmp_right_fp8,&tmp_right_fp8,&A->x);
  fp8_mul_mpn(&tmp_ax,&A->x,curve_a.x0);

  fp8_add(&tmp_right_fp8,&tmp_right_fp8,&tmp_ax);

  if(fp8_cmp(&tmp_right_fp8, &tmp_left_fp8)==0){
    printf("efp8 check on curve: On curve\n");
  }else{
    printf("efp8 check on curve: NOT On curve\n");
  }

}

void efp8_checkOnTwsitCurve(efp8_t *A){
  static fp_t twistedb;
  static fp8_t tmp_left_fp8, tmp_right_fp8,tmp_ax;
  fp_set_ui(&twistedb,4);

  fp8_sqr(&tmp_left_fp8,&A->y);

  fp8_sqr(&tmp_right_fp8,&A->x);
  fp8_mul(&tmp_right_fp8,&tmp_right_fp8,&A->x);
  fp8_mul_mpn(&tmp_ax,&A->x,curve_a.x0);
  fp8_add(&tmp_right_fp8,&tmp_right_fp8,&tmp_ax);

  if(fp8_cmp(&tmp_right_fp8, &tmp_left_fp8)==0){
    printf("efp8 check on curve: On curve\n");
  }else{
    printf("efp8 check on curve: NOT On curve\n");
  }

}
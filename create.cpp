#include "create.h"
#include "define.h"
#include "fp.h"
#include "efp.h"
#include "fp2.h"
#include "mpn.h"
#include <gmp.h>


void create_prt(){
  //c=22,HW=6
  mpz_set_str(X_z,"ffc0000004020002",16);
  mpz_set_str(prime_z,"347111bfc75e57d130de7be68437c8d75455d209459d421455023bee14df9fe75aa4734686ca3d08c1fa594100d79421d56c53899ee0f066fad9eb45b0985dbdbba2dcc1",16);
  //p = ((4*(hy**2)+4) + 4*x + ((4*hy)+1)*(x**2) - (4*hy)*(x**3) + (8*(hy**2)+5)*(x**4) + ((4*hy)+1)*(x**6) - (4*hy)*(x**7) + (1+4*(hy**2))*(x**8) )//4
  if(!mpz_probab_prime_p(prime_z,30))printf("Inputed p*(prime_z) is not a prime");

  mpz_set_str(order_z,"ff005ff010fbfd093a41afce5a026f5b55729902a9b26b3179c1806080400011",16);
  mpz_set_str(trace_z,"ff005ff010fbfd093a41afce5a026f5b55729902a9b26b307a0180607c3e000e",16);
  //trace_z is a negative value
  const unsigned char* xai = reinterpret_cast<const unsigned char *>("02");
  mpn_set_str(&X,xai,sizeof(char)*2,16); //ui(&X,1,319014718988379808906617884108577046528);

  mpn_set_mpz(prime,prime_z);
  mpn_mul_n(prime2,prime,prime,FPLIMB);
  mpn_lshift(few_prime,prime,FPLIMB,31);
  gmp_printf("few_prime=%Nx\n",few_prime,FPLIMB);
  // mpn_mul_1(few_prime,prime,1,X);

  gmp_printf("X     (%4dbit length) = %Zd\n",(int)mpz_sizeinbase(X_z,2),X_z);
  gmp_printf("prime (%4dbit length) = %Zd\n",(int)mpz_sizeinbase(prime_z,2),prime_z);
  gmp_printf("order (%4dbit length) = %Zd\n",(int)mpz_sizeinbase(order_z,2),order_z);
  gmp_printf("trace (%4dbit length) = %Zd\n",(int)mpz_sizeinbase(trace_z,2),trace_z);
  printf("X     (HW :%2ld)(binary) = ",mpz_popcount(X_z));
  mpz_out_str(stdout,2,X_z);printf("\n");
  printf("trace (HW :%2ld)(binary) = ",mpz_popcount(trace_z));
  mpz_out_str(stdout,2,trace_z);printf("\n");
  fp_set_ui(&base_c.x0,1);
  fp_set_neg(&base_c.x1,&base_c.x0);

  gmp_printf("\nmodulo polynomial\n");

  gmp_printf("fp2 = fp[alpha] /(alpha^p + alpha + 1 )\n");
  gmp_printf("fp4 = fp2[beta] /(beta^2  - (1,-1))\n");
  gmp_printf("fp8 = fp4[gamma]/(gamma^2 - beta )\n");
  // fp_println("base_c     = ",&base_c);
  // fp_println("base_c_inv = ",&base_c_inv);
  printf("---------------------------------\n");
  miller_loop_v =    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  finalexp_pow_x =   {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  finalexp_pow_x_2 = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  finalexp_pow_hy =  {0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1};
  finalexp_pow_4hy = {0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1};
  mpn_sub_ui(p_1.x0,prime,FPLIMB,1);
  fp_set_ui(&fp2_neg_1.x0,1);
  fp_set_ui(&fp2_neg_1.x1,1);
  fp2_inv(&base_c_inv, &base_c);
  fp2_to_montgomery(&base_c_invMR, &base_c_inv);
  fp2_println("base_c_inv", &base_c_inv);
  fp_set_ui(&three_1, 3);
  fp_inv(&three_1,&three_1);
  fp_set_neg(&three_1,&three_1);
  fp_to_montgomery(&three_1MR,&three_1);

}

void check_base(){
  fp_t tmp;
  fp2_t tmp2;
  fp_init(&tmp);
  fp2_init(&tmp2);
  mpz_t expo;
  mpz_init(expo);
  std::string x = "7864352814187761672954208385054718516586573559536409605479181947936697562602380972863469522882765880267143715565140764573573332968082475950434711738861351634506540";
  std::string y = "7864352814187761672954208385054718516586573559536409605479181947936697562602380972863469522882765880267143715565140764573573332968082475950434711738861351634506538";

  mpn_set_char(tmp2.x0.x0, FPLIMB, x.c_str());
  mpn_set_char(tmp2.x1.x0, FPLIMB, y.c_str());
  fp2_println("base_c", &tmp2);
  fp2_inv(&tmp2, &tmp2);
  fp2_println("base_c_inv", &tmp2);

  //check base_c = QNR
  // fp_set(&tmp,&base_c);

  // mpz_sub_ui(expo,prime_z,1);
  // mpz_divexact_ui(expo,expo,2);
  // fp2_pow(&tmp,&base_c,expo);
  // if(fp_cmp_one(&tmp)==0) printf("error!!! c^((p-1)/2)==1\n\n");

  // mpz_set_ui(expo,2);
  // fp_pow(&tmp,&tmp,expo);
  // // fp_println("c^(p-1) :",&tmp);
  // if(fp_cmp_one(&tmp)!=0) printf("error!!! c^(p-1)!=1\n\n");
  
  // //check base_c = QNR
  // fp_set_ui(&tmp2.x0,1);
  // // fp_set_neg(&tmp2.x0, &tmp2.x0);
  // fp_set_ui(&tmp2.x1,1);
  // fp_set_neg(&tmp2.x1, &tmp2.x1);
  mpz_pow_ui(expo,prime_z,2);
  mpz_sub_ui(expo,expo,1);
  mpz_divexact_ui(expo,expo,2);
  fp2_pow(&tmp2,&tmp2,expo);
  if(fp2_cmp_one(&tmp2)==0) printf("error!!! alpha^((p^2-1)/2)==1\n\n");

  mpz_set_ui(expo,2);
  fp2_pow(&tmp2,&tmp2,expo);
  if(fp2_cmp_one(&tmp2)!=0) printf("error!!! alpha^(p^2-1)!=1\n\n"),fp2_println("fp2_base:",&tmp2);

  // fp_mul(&tmp,&base_c,&base_c_inv);
  //fp_println("base_c * base_c_inv = ",&tmp);
  // if(fp2_cmp_one(&tmp)!=0) printf("error!!! base_c * base_c_inv!=1\n\n");

  mpz_clear(expo);
}

void frobenius_precalculation(){
  fp2_t tmp;
  mpz_t expo;
  fp2_init(&tmp);
  mpz_init(expo);

  // mpz_mul(expo,prime_z,prime_z);
  mpz_sub_ui(expo,prime_z,1);
  mpz_divexact_ui(expo,expo,4);

  fp2_pow(&tmp,&base_c,expo);
  fp2_set(&frobenius_1_4,&tmp);
  fp2_println("c^(p-1)/4", &frobenius_1_4);
  fp2_to_montgomery(&frobenius_1_4MR, &frobenius_1_4);

  mpz_set_ui(expo,2);
  fp2_pow(&frobenius_2_4,&tmp,expo);
  fp2_println("c^2(p-1)/4", &frobenius_2_4);
  fp2_to_montgomery(&frobenius_2_4MR, &frobenius_2_4);


  mpz_set_ui(expo,3);
  fp2_pow(&frobenius_3_4,&tmp,expo);
  fp2_println("c^3(p-1)/4", &frobenius_3_4);
  fp2_to_montgomery(&frobenius_3_4MR, &frobenius_3_4);

  mpz_clear(expo);
  printf("Frobenius precalculation is done\n");
}

void curve_search(){
  //s = p + 1 - t ???????????????
  // efp_t efp_P,efp_sP;
  // efp_init(&efp_P);
  // efp_init(&efp_sP);
  // mpz_t s;
  // mpz_init(s);
  // mpz_add_ui(s,prime_z,1);
  // mpz_add(s,s,trace_z);

  fp_set_ui(&curve_a,1);
  // fp_set_neg(&curve_a,&curve_a);

  // while(1){
  //   fp_add_ui(&curve_a,&curve_a,1);
  //   efp_rational_point(&efp_P);
  //   efp_scm(&efp_sP,&efp_P,s);
  //   if(efp_sP.infinity==1)  break;
  // }

  // mpz_clear(s);
  fp_println("curve_a:", &curve_a);

  printf("Elliptic curve search is done\n");
}

void frobenius_trace(mpz_t *trace,unsigned int m){
  mpz_t t_m[33];
  mpz_t t_2,tmp1,tmp2;
  for(int i=0;i<33;i++) mpz_init(t_m[i]);
  mpz_init(t_2);
  mpz_init(tmp1);
  mpz_init(tmp2);

  mpz_set(t_m[0],trace_z);
  //t_2=(t^2 -2p)?????????
  mpz_mul(t_2,trace_z,trace_z);
  mpz_mul_ui(tmp1,prime_z,2);
  mpz_sub(t_2,t_2,tmp1);
  
  mpz_set(t_m[1],t_2);

  for(int i=2;i<m;i++){
    mpz_mul(tmp1,trace_z,t_m[i-1]);
    mpz_mul(tmp2,prime_z,t_m[i-2]);
    mpz_sub(t_m[i],tmp1,tmp2);
  }
  mpz_set(*trace,t_m[m-1]);

  for(int i=0;i<33;i++) mpz_clear(t_m[i]);
  mpz_clear(t_2);
  mpz_clear(tmp1);
  mpz_clear(tmp2);
}

void efpm_order(mpz_t *order,unsigned int m){
  mpz_t trace;
  mpz_init(trace);

  frobenius_trace(&trace,m);
  mpz_pow_ui(*order,prime_z,m);
  mpz_add_ui(*order,*order,1);
  if(m==1){
    mpz_add(*order,*order,trace);
  }else{
    mpz_sub(*order,*order,trace);
  }

  mpz_clear(trace);
}

void create_weil(){
  efpm_order(&efp_total,1);
  efpm_order(&efp2_total,2);
  efpm_order(&efp4_total,4);
  efpm_order(&efp8_total,8);
  mpz_t temp;
  mpz_init(temp);

  frobenius_trace(&temp,8);
  mpz_add(fp8_total_r,efp8_total,temp);
  mpz_sub_ui(fp8_total_r,fp8_total_r,2);
  mpz_divexact(fp8_total_r,fp8_total_r,order_z);
  // gmp_printf("fp8_total_r:%Zd\n",fp8_total_r);

  mpz_set(miller_loop_s,X_z);
  //X = 1 (mod 2) ?????????
  //(X+1)/2?????????????????????????????????
  // mpz_set_str(hy_neg,"dc04",16);
  // mpz_set(fourhy_neg,hy_neg);
  // mpz_mul_ui(fourhy_neg,fourhy_neg,4);
  // mpz_mul(hy_neg,fourhy_neg,hy_neg);
  // mpz_set_str(X_1,"",16);
  mpz_set(X_2,X_z);
  mpz_mul(X_2,X_2,X_2);

  // gmp_printf("p-1:%Nu",p_1.x0,FPLIMB);
}

void tmp_init(){
  mpz_init(X_z);
  mpz_init(prime_z);
  mpz_init(order_z);
  mpz_init(trace_z);

  mpz_init(efp_total);
  mpz_init(efp2_total);
  mpz_init(efp4_total);
  mpz_init(efp8_total);
  mpz_init(fp8_total_r);

  mpz_init(miller_loop_s);
  mpz_init(X_1_div2);
  mpz_init(X_1);
  mpz_init(X_2);
  mpz_init(X_2_1);
  mpz_init_set_ui(four,4);

  mpz_init(hardpart);
  mpz_init(hy_neg);
  mpz_init(fourhy_neg);
  mpz_init_set_ui(four,4);
  mpz_init_set_ui(three,3);

}
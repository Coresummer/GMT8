#include "gmp.h"
#include <vector>
#include <iostream>
using namespace std;

int main(){
    mpz_t xai,fourhy,x_2,tmp;
    mpz_init_set_str(xai,"ffc00020fffffffc",16);
    mpz_init_set(x_2,xai);
    mpz_mul(x_2,xai,xai);
    
    mpz_init_set_str(fourhy,"dc04",16);
    // mpz_mul_ui(fourhy,fourhy,4);
    vector<int> v,v2,v3;

    mpz_init_set_ui(tmp,0);
    while(mpz_cmp_ui(xai, 0 )!= 0){
        mpz_mod_ui(tmp,xai,2);
        if(mpz_cmp_ui(tmp, 1)==0){
            mpz_mod_ui(tmp,xai,4);
            int res = 2-(mpz_get_si(tmp));
            v.push_back(res);
            mpz_set_si(tmp,res);
            mpz_sub(xai,xai,tmp);
        }else{
            v.push_back(0);
        }
        mpz_div_ui(xai,xai,2);
    }
    cout << "xai:";
    for (int i =0; i < v.size(); i++){
        cout << v[i] << ", " ;
    }
    cout << endl;
    cout << endl;

    // mpz_init_set_ui(tmp,0);
    while(mpz_cmp_ui(x_2, 0 )!= 0){
        mpz_mod_ui(tmp,x_2,2);
        if(mpz_cmp_ui(tmp, 1)==0){
            mpz_mod_ui(tmp,x_2,4);
            int res = 2-(mpz_get_si(tmp));
            v2.push_back(res);
            mpz_set_si(tmp,res);
            mpz_sub(x_2,x_2,tmp);
        }else{
            v2.push_back(0);
        }
        mpz_div_ui(x_2,x_2,2);
    }
    cout << "x_2:";
    for (int i =0; i < v2.size(); i++){
        cout << v2[i] << ", " ;
    }
    cout << endl;
    cout << endl;

        // mpz_init_set_ui(tmp,0);
    while(mpz_cmp_ui(fourhy, 0 )!= 0){
        mpz_mod_ui(tmp,fourhy,2);
        if(mpz_cmp_ui(tmp, 1)==0){
            mpz_mod_ui(tmp,fourhy,4);
            int res = 2-(mpz_get_si(tmp));
            v3.push_back(res);
            mpz_set_si(tmp,res);
            mpz_sub(fourhy,fourhy,tmp);
        }else{
            v3.push_back(0);
        }
        mpz_div_ui(fourhy,fourhy,2);
    }
    cout << "4hy:";
    for (int i =0; i < v3.size(); i++){
        cout << v3[i] << ", " ;
    }
    cout << endl;
    cout << endl;


    mpz_clear(fourhy);
    mpz_clear(xai);
    mpz_clear(x_2);
    mpz_clear(tmp);

    
    return 0;

}
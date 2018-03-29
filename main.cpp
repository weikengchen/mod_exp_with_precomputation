#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "libfastexp.h"
#include <gmp.h>
#include <x86intrin.h>

#include <ratio>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

void crypto_random_make_mpz(mpz_t ret, unsigned int l){
    unsigned char r[l / 8]; unsigned short tmp;
    int i, j;
    for(i = 0; i < l / 8 - 4; i+=4){
        _rdrand32_step((unsigned int *)&r[i]);
    }
    for(j = i; j < l / 8; j++){
        _rdrand16_step(&tmp);
        r[j] = tmp % 256;
    }
    mpz_import(ret, l / 8, 1, 1, 0, 0, r);
}

void crypto_random_make_bytes(unsigned char *ret, unsigned int l){
    unsigned short tmp;
    int i, j;
    for(i = 0; i < l / 8 - 4; i+=4){
        _rdrand32_step((unsigned int *)&ret[i]);
    }
    for(j = i; j < l / 8; j++){
        _rdrand16_step(&tmp);
        ret[j] = tmp % 256;
    }
}


int main(){
    mpz_t base, mod, result, result2;
    mpz_inits(base, mod, result, result2, NULL);
    
    fastexp_state _st;
    fastexp_state *st = &_st;
    gmp_sscanf("68924559023601125404348155141634098250387986354880927979218529916307435006751819342929439585847854697250347846226166345365324771618545449092289365568135252406023077502345907798723942114061110788740849932509820253584778513830907678795428771468857557448308625068838328972299648270406456374247453331154798129347513308192535733101891725052708698800101271604670214572285811151721752881301034667407375466204927595946092363345480217273608589869269843686131697067099117151180164016214542289791609070891445117545843219290238219595566979697203618784248353290929981466727000848038365556543828613008543841755579751414975516321323", "%Zd", mod);
        
    gmp_sscanf("3", "%Zd", base);
    
    fastexp_prepare(base, mod, 224, 14, st);
    
    unsigned char exp[1000][224 / 8];
    mpz_t exp_mpz[1000];
    
    for(int i = 0; i < 1000; i++){
        memset(exp[i], 0, sizeof(unsigned char) * (224 / 8));
        crypto_random_make_bytes(exp[i], 224);
        mpz_init(exp_mpz[i]);
        mpz_import(exp_mpz[i], 224 / 8, 1, 1, 0, 0, exp[i]);
    }
    
    {
    auto t1 = Clock::now();
    for(int i = 0; i < 1000; i++){
        fastexp_compute(result, exp[i], st);
    }
    auto t2 = Clock::now();
    printf("%lf ms per computation.\n", (std::chrono::duration<double, std::milli>(t2 - t1).count())/1000.0);
    }
    
    {
    auto t1 = Clock::now();
    for(int i = 0; i < 1000; i++){
        mpz_powm_sec(result2, base, exp_mpz[i], mod);
    }
    auto t2 = Clock::now();
    printf("%lf ms per computation.\n", (std::chrono::duration<double, std::milli>(t2 - t1).count())/1000.0);
    }
      
    fastexp_release(st);
    
    return 0;
}

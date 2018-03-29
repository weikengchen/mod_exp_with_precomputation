#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <gmp.h>
#include "libfastexp.h"


void fastexp_prepare(mpz_t base, mpz_t mod, int exp_len, int group_slice_len, fastexp_state* st){
    int num_group;
    num_group = exp_len / group_slice_len;
    
    int group_size;
    group_size = 1 << group_slice_len;
    
    st->exp_len = exp_len;
    st->group_slice_len = group_slice_len;
    st->num_group = num_group;
    st->group_size = group_size;
    st->precompute_table = (mpz_t*) malloc(sizeof(mpz_t) * num_group * group_size);
    
    mpz_init(st->mod);
    mpz_set(st->mod, mod);
    
    printf("size of precompute table: %lu bytes.\n", sizeof(mpz_t) * num_group * group_size);
    
    mpz_t cur;
    mpz_t curmultiply;
    mpz_inits(cur, curmultiply, NULL);
    mpz_set_ui(cur, (unsigned long) 1);
    mpz_set(curmultiply, base);
    
    int counter = 0;
    for(int i = 0; i < num_group; i++){
        for(int j = 0; j < group_size; j++){
            mpz_init(st->precompute_table[counter]);
            mpz_set(st->precompute_table[counter], cur);
            mpz_mul(cur, cur, curmultiply);
            mpz_mod(cur, cur, mod);
            counter++;
        }
        mpz_set(curmultiply, cur);
        mpz_set_ui(cur, (unsigned long) 1);
    }
}

void fastexp_compute(mpz_t result, unsigned char *exp, fastexp_state *st){
    /* exp must be an array of exp_len bits, or exp_len/8 bytes */
    /* highest bits first */ 
    
    unsigned char *expanded_exp = malloc(sizeof(unsigned char) * st->exp_len);
    unsigned char *expanded_exp_rev = malloc(sizeof(unsigned char) * st->exp_len);
    for(int i = 0; i < (st->exp_len) / 8; i++){
        expanded_exp[i * 8 + 0] = (exp[i] >> 7) & 1;
        expanded_exp[i * 8 + 1] = (exp[i] >> 6) & 1;
        expanded_exp[i * 8 + 2] = (exp[i] >> 5) & 1;
        expanded_exp[i * 8 + 3] = (exp[i] >> 4) & 1;
        expanded_exp[i * 8 + 4] = (exp[i] >> 3) & 1;
        expanded_exp[i * 8 + 5] = (exp[i] >> 2) & 1;
        expanded_exp[i * 8 + 6] = (exp[i] >> 1) & 1;
        expanded_exp[i * 8 + 7] = exp[i] & 1;
    }
    
    for(int i = 0; i < st->exp_len; i++){
        expanded_exp_rev[i] = expanded_exp[st->exp_len - 1 - i];
    }
    
    mpz_set_ui(result, (unsigned long) 1);
    
    int counter = 0;
    int step;
    for(int i = 0; i < st->num_group; i++){
        step = expanded_exp_rev[counter];
        //printf("i=%d, expanded_exp_rev[%d]=%d\n", i, counter, step);
        counter++;
        for(int j = 1; j < st->group_slice_len; j++){
            step = step + (expanded_exp_rev[counter] << j);
            //printf("expanded_exp_rev[%d]=%d\n", counter, expanded_exp_rev[counter]);
            counter++;
        }
        
        //printf("step = %d\n", step);
        
        mpz_mul(result, result, st->precompute_table[(i * st->group_size) + step]);
        mpz_mod(result, result, st->mod);
    }
    
    free(expanded_exp);
    free(expanded_exp_rev);
}

void fastexp_release(fastexp_state *st){
    mpz_clear(st->mod);
    int counter = 0;
    for(int i = 0; i < st->num_group; i++){
        for(int j = 0; j < st->group_size; j++){
            mpz_clear(st->precompute_table[counter]);
            counter++;
        }
    }
    free(st->precompute_table);
    st->precompute_table = NULL;
    
    st->exp_len = 0;
    st->group_slice_len = 0;
    st->num_group = 0;
    st->group_size = 0;
}
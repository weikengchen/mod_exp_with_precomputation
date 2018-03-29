#ifndef LIBFASTEXP_H
#define LIBFASTEXP_H

#include <gmp.h>

typedef struct{
    int exp_len;
    int group_slice_len;
    int num_group;
    int group_size;
    mpz_t mod;
    mpz_t *precompute_table;
} fastexp_state;
    

#ifdef __cplusplus
extern "C" {
#endif
    
void fastexp_prepare(mpz_t base, mpz_t mod, int exp_len, int group_slice_len, fastexp_state* st);

void fastexp_compute(mpz_t result, unsigned char *exp, fastexp_state *st);

void fastexp_release(fastexp_state *st);
    
#ifdef __cplusplus
}
#endif
#endif //LIBFASTEXP_H

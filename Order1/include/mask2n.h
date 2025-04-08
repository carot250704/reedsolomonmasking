#ifndef MASK16_H
#define MASK16_H

#include <stdint.h>

typedef struct mask2n_t 
{
    uint8_t m[3];
} mask2n_t;

uint8_t mul2n(uint8_t a, uint8_t b);

void mask_2n(mask2n_t *rep, uint8_t elem);
void unmask_2n(uint8_t *elem, mask2n_t *rep);

void affect_2n(mask2n_t *out, mask2n_t *in);
void op_xor_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2);
void op_xor3_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2, mask2n_t *x3);
void op_and_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2);
void op_not_2n(mask2n_t *rep, mask2n_t *x);

void mask_2n_byte(mask2n_t *out, uint8_t *in, uint16_t in_len);
void unmask_2n_byte(uint8_t *out, mask2n_t *in, uint16_t out_len);
void refresh_2n(mask2n_t *out, mask2n_t *in);

#endif

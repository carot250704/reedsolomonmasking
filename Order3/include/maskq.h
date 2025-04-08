#ifndef MASKQ_H
#define MASKQ_H

#include <stdint.h>

typedef struct maskq_t
{
    int16_t m[8];
} maskq_t;


int16_t add_reduce(int16_t x, int16_t y);
int16_t sub_reduce(int16_t x, int16_t y);
int16_t fq_mul_reduce(int16_t x, int16_t y);


void mask_q(maskq_t *rep, int16_t elem);
void unmask_q(int16_t *elem, maskq_t *rep);
void affect_q(maskq_t *out, maskq_t *in);
void op_xor_q(maskq_t *rep, maskq_t *x1, maskq_t *x2);
void op_xor5_q(maskq_t *rep, maskq_t *x1, maskq_t *x2, maskq_t *x3, maskq_t *x4, maskq_t *x5);
void op_not_q(maskq_t *rep, maskq_t *x);

void op_add_q(maskq_t *rep, maskq_t *x1, maskq_t *x2);
void op_sub_q(maskq_t *rep, maskq_t *x1, maskq_t *x2);
void op_and_q(maskq_t *rep, maskq_t *x1, maskq_t *x2);
void op_mul_q(maskq_t *rep, maskq_t *x1, maskq_t *x2);
void op_mul_scal_q(maskq_t *rep, maskq_t *x, int16_t scal);

void refresh_q(maskq_t *out, maskq_t *in);

void mask_q_byte(maskq_t *out, uint8_t *in, uint64_t in_len);
void unmask_q_byte(uint8_t *out, maskq_t *in, uint64_t out_len);

void get_bit(maskq_t *out, maskq_t *in);
#endif
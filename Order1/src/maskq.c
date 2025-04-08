#include "maskq.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>





int16_t barrett_reduce(int32_t a) {
    int16_t t;
    const int16_t v = 20158;

    t  = (int64_t)((int64_t)v* (int64_t)a) >> 26;
    t *= 3329;
    t = (a - t);
    return (t + 3329) % 3329;
}


int16_t fq_mul_reduce(int16_t x, int16_t y)
{
    int32_t temp;
    temp = (x * y);
    return barrett_reduce(temp);
}

int16_t add_reduce(int16_t x, int16_t y)
{
    return (x + y + 3329) % 3329;
}

int16_t sub_reduce(int16_t x, int16_t y)
{
    return (x - y + 3329 + 3329) % 3329;
}


void DFT_mod2(int16_t out[4], int16_t in[4])
{
    out[0] = add_reduce(in[0], in[2]);
    out[1] = add_reduce(in[1], in[3]);

    out[2] = sub_reduce(in[0], in[2]);
    out[3] = sub_reduce(in[1], in[3]);
}


void DFT_mod1(int16_t out[4], int16_t in[4])
{
    out[0] = add_reduce(in[0], in[1]);
    out[1] = sub_reduce(in[0], in[1]);

    out[2] = add_reduce(in[2], fq_mul_reduce(in[3], 1729));
    out[3] = add_reduce(in[2], fq_mul_reduce(in[3], 1600));
}


void DFT(int16_t out[4],int16_t in[4])
{
    int16_t temp[4];
    
    DFT_mod2(out, in);
    DFT_mod1(temp, out);

    out[0] = temp[0];
    out[1] = temp[2];
    out[2] = temp[1];
    out[3] = temp[3];
}


void IDFT(int16_t out[4],int16_t in[4])
{
    int16_t temp[4];

    DFT_mod2(out, in);
    DFT_mod1(temp, out);

    out[0] = fq_mul_reduce(temp[0], 2497);
    out[1] = fq_mul_reduce(temp[3], 2497);
    out[2] = fq_mul_reduce(temp[1], 2497);
    out[3] = fq_mul_reduce(temp[2], 2497);
}

// mask of element in Z/3329Z
void mask_q(maskq_t *m, int16_t in)
{
    int16_t temp1[2];
    int16_t temp4[4];

    temp1[0] = in;
    temp1[1] = barrett_reduce(rand() & 0xFFFF);   

    temp4[0] = add_reduce(fq_mul_reduce(temp1[0],  824), fq_mul_reduce(2506, temp1[1]));
    temp4[1] = add_reduce(fq_mul_reduce(temp1[0], 2608), fq_mul_reduce(721 , temp1[1]));
    temp4[2] = 0;
    temp4[3] = 0;
  
    DFT(m->m,temp4); 
}

void op_mul_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    int16_t temp4[4];
    for(uint8_t i = 0; i < 4; i++) {
        temp4[i] = fq_mul_reduce(x1->m[i], x2->m[i]);
    }

    IDFT(rep->m, temp4);

    uint16_t temp1 = barrett_reduce(rand() & 0xFFFF);
    temp4[0] = fq_mul_reduce(temp1, 2506);
    temp4[1] = fq_mul_reduce(temp1, 721);
    temp1 = 0;

    temp4[0] = add_reduce(temp4[0], add_reduce(rep->m[0], fq_mul_reduce(rep->m[2], 2379)));
    temp4[1] = add_reduce(temp4[1], add_reduce(rep->m[1], fq_mul_reduce(rep->m[3], 2379)));
    temp4[2] = 0;
    temp4[3] = 0;

    DFT(rep->m, temp4);
}

// unmask  of element in Z/3329Z^4
void unmask_q(int16_t *out,maskq_t *in)
{
    int16_t temp4[4];
    IDFT(temp4, in->m);
    *out = add_reduce(temp4[0], fq_mul_reduce(temp4[1], 2970)) ;   
}


void affect_q(maskq_t *out, maskq_t *in)
{
    for(int i = 0; i < 4; i++)
        out->m[i] = in->m[i];
}

void op_add_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    for (int i=0; i < 4; i++)
        rep->m[i] = add_reduce(x1->m[i], x2->m[i]);
}

void op_sub_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    for (int i=0; i < 4; i++)
        rep->m[i] = sub_reduce(x1->m[i], x2->m[i]);
}

void op_xor_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    maskq_t rep1;
    op_mul_q(&rep1, x1, x2);
    for(uint8_t i = 0; i < 4; i++) {
        rep1.m[i] = fq_mul_reduce(rep1.m[i], 3327);
    }

    op_add_q(rep, &rep1, x1);
    op_add_q(rep, rep, x2);
}

void op_and_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    op_mul_q(rep, x1, x2);
}

void op_not_q(maskq_t *rep, maskq_t *x)
{
    maskq_t one;
    mask_q(&one, 1);

    op_sub_q(rep, &one, x);
}

void op_mul_scal_q(maskq_t *rep, maskq_t *x, int16_t scal)
{
    for(uint8_t i = 0; i < 4; i++) {
        rep->m[i] = fq_mul_reduce(x->m[i], scal);
    }
}

void mask_q_byte(maskq_t *out, uint8_t *in, uint64_t in_len)
{
    for(uint64_t i = 0; i < in_len; i++) {
        for(uint64_t j = 0; j < 8; j++) {
            mask_q(&out[8*i + j], (int32_t)((in[i] >> j) & 1));
        }
    }
}

void unmask_q_byte(uint8_t *out, maskq_t *in, uint64_t out_len)
{
    int16_t temp;
    for(uint64_t i = 0; i < out_len; i++) {
        out[i] = 0;
        for(uint8_t j = 0; j < 8; j++) {
            unmask_q(&temp, &in[8*i + j]);
            out[i] |= temp<<j;
        }
    }
}

void refresh_q(maskq_t *out, maskq_t *in)
{
    maskq_t m;
    mask_q(&m, 0);
    op_add_q(out, &m, in);
}

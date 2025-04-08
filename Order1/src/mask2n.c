#include "mask2n.h"

#include <stdlib.h>
#include <stdio.h>


static uint8_t mul[16][16] = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15},
    { 0,  2,  4,  6,  8, 10, 12, 14,  3,  1,  7,  5, 11,  9, 15, 13},
    { 0,  3,  6,  5, 12, 15, 10,  9, 11,  8, 13, 14,  7,  4,  1,  2},
    { 0,  4,  8, 12,  3,  7, 11, 15,  6,  2, 14, 10,  5,  1, 13,  9},
    { 0,  5, 10, 15,  7,  2, 13,  8, 14, 11,  4,  1,  9, 12,  3,  6},
    { 0,  6, 12, 10, 11, 13,  7,  1,  5,  3,  9, 15, 14,  8,  2,  4},
    { 0,  7, 14,  9, 15,  8,  1,  6, 13, 10,  3,  4,  2,  5, 12, 11},
    { 0,  8,  3, 11,  6, 14,  5, 13, 12,  4, 15,  7, 10,  2,  9,  1},
    { 0,  9,  1,  8,  2, 11,  3, 10,  4, 13,  5, 12,  6, 15,  7, 14},
    { 0, 10,  7, 13, 14,  4,  9,  3, 15,  5,  8,  2,  1, 11,  6, 12},
    { 0, 11,  5, 14, 10,  1, 15,  4,  7, 12,  2,  9, 13,  6,  8,  3},
    { 0, 12, 11,  7,  5,  9, 14,  2, 10,  6,  1, 13, 15,  3,  4,  8},
    { 0, 13,  9,  4,  1, 12,  8,  5,  2, 15, 11,  6,  3, 14, 10,  7},
    { 0, 14, 15,  1, 13,  3,  2, 12,  9,  7,  6,  8,  4, 10, 11,  5},
    { 0, 15, 13,  2,  9,  6,  4, 11,  1, 14, 12,  3,  8,  7,  5, 10}

};

uint8_t mul2n(uint8_t a, uint8_t b)
{
    return mul[a][b];
}


static void mul_C(uint8_t out[3], uint8_t in[2])
{
    out[0] = mul[0xe][in[0]] ^ mul[0xf][in[1]];
    out[1] = mul[0xb][in[0]] ^ mul[0xa][in[1]];         
    out[2] = mul[0x6][in[0]] ^ mul[0x7][in[1]];            
}

static void mul_vand(uint8_t out[5], uint8_t in[5])
{
    out[0] = in[0] ^ in[1] ^ in[2];
    out[1] = in[0] ^ mul[in[1]][0x6] ^ mul[in[2]][0x7];
    out[2] = in[0] ^ mul[in[1]][0x7] ^ mul[in[2]][0x6];
}

static void mul_inv_vand(uint8_t out[5], uint8_t in[5])
{
    out[0] = in[0] ^ in[1] ^ in[2];
    out[1] = in[0] ^ mul[in[1]][0x7] ^ mul[in[2]][0x6];
    out[2] = in[0] ^ mul[in[1]][0x6] ^ mul[in[2]][0x7];
}


void mask_2n(mask2n_t *rep, uint8_t elem)
{
    uint8_t temp[2];
    temp[0] = elem   & 0x0f;
    temp[1] = rand() & 0x0f; 
    
    mul_C(rep->m, temp);
}

void unmask_2n(uint8_t *elem, mask2n_t *rep)
{
    *elem = mul[rep->m[0]][0x9] ^ mul[rep->m[1]][0xc] ^ mul[rep->m[2]][0x4];
}


void affect_2n(mask2n_t *out, mask2n_t *in)
{
    for(uint8_t i = 0; i < 3; i++) {
        out->m[i] = in->m[i];
    }
}

void op_xor_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2)
{
    for(uint8_t i = 0; i < 3; i++) {
        rep->m[i] = x1->m[i] ^ x2->m[i];
    }
}

void op_xor3_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2, mask2n_t *x3)
{
    for(uint8_t i = 0; i < 3; i++) {
        rep->m[i] = x1->m[i] ^ x2->m[i] ^ x3->m[i];
    }
}

void op_and_2n(mask2n_t *rep, mask2n_t *x1, mask2n_t *x2)
{
    mask2n_t temp;
    mask_2n(&temp, 0);
    for(uint8_t i = 0; i < 3; i++) {
        temp.m[i] ^= mul[x1->m[i]][x2->m[i]];
    }

    //IDFT
    mul_inv_vand(rep->m, temp.m);

    
    uint8_t r = rand() & 0x0f;
    temp.m[0] = rep->m[0];
    temp.m[0] ^= mul[r][0x2];
    temp.m[0] ^= mul[rep->m[2]][0xc];

    temp.m[1] = rep->m[1] ^ mul[r][0xd];
    temp.m[2] = 0;

    mul_vand(rep->m, temp.m);
}

void op_not_2n(mask2n_t *rep, mask2n_t *x)
{
    mask2n_t m;
    mask_2n(&m, 1);
    op_xor_2n(rep, x, &m);
}

// in = bytes, out = bits masked
void mask_2n_byte(mask2n_t *out, uint8_t *in, uint16_t in_len)
{
    for(uint64_t i = 0; i < in_len; i++) {
        for(uint8_t j = 0; j < 8; j++) {
            mask_2n(&out[8*i + j], ((in[i]) >> j) & 1);
        }
    }
}

// in = bits masked, out = bytes
void unmask_2n_byte(uint8_t *out, mask2n_t *in, uint16_t out_len)
{
    uint8_t temp;
    for(uint16_t i = 0; i < out_len; i++) {
        out[i] = 0;
        for(uint8_t j = 0; j < 8; j++) {
            unmask_2n(&temp, &in[8*i + j]);
            out[i] |= temp<<j;
        }
    }
}


void refresh_2n(mask2n_t *out, mask2n_t *in)
{
    mask2n_t m;
    mask_2n(&m, 0);
    op_xor_2n(out, in, &m);
}
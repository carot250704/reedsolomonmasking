#include "maskq.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/**********************************  FFT for n = 8 *************************************************

q = 3329;
q-1 = 2^8*13;
nu = 3;
root = ω = 3^{(q-1)/8}=749 = -2580;
ω^2 = 1729 = -1600
ω^3 = 40 = -3289
ω^4 = 3328 = -1
ω^5 = 2580 = -749
ω^6 = 1600 = -1729
ω^7 = 3289 = -40

α = 2970
α^4 = 341;
          1   1   1   1       1    1    1    1 
          α   α^2 α^3 α^4    2970 2379 1492  341 
Vand(α) = α^2 α^4 α^6 α^8  = 2379  341 2292 3095
          α^3 α^6 α^9 α^12   1492 2292  781 1812

            [1306 2847 2865 2852]
            [3259  144  125 1694]
Vand(α)^1 = [1372 1731 2033  320]
            [ 722 1936 1635 1792]


Wij = {ω^{-4i}} = {1,3328,1,3328,1,3328,1,3328}


          [ 1    1    1    1    1    1    1    1]
          [ 1  749 1729   40 3328 2580 1600 3289]
          [ 1 1729 3328 1600    1 1729 3328 1600]
Vand(ω) = [ 1   40 1600  749 3328 3289 1729 2580]
          [ 1 3328    1 3328    1 3328    1 3328]
          [ 1 2580 1729 3289 3328  749 1600   40]
          [ 1 1600 3328 1729    1 1600 3328 1729]
          [ 1 3289 1600 2580 3328   40 1729  749]


               [2913 2913 2913 2913 2913 2913 2913 2913]
               [2913 3324  200 1987  416    5 3129 1342]
               [2913  200  416 3129 2913  200  416 3129]
Vand(ω)^(-1) = [2913 1987 3129 3324  416 1342  200    5]
               [2913  416 2913  416 2913  416 2913  416]
               [2913    5  200 1342  416 3324 3129 1987]
               [2913 3129  416  200 2913 3129  416  200]
               [2913 1342 3129    5  416 1987  200 3324]


*************************************************************************************************/


uint64_t rand64(void){ 
    static uint64_t oneway, rotor0, rotor1, feedback, seed; 
    if(!oneway){ 
        seed = rand();
        rotor0 =seed; 
        rotor1 =seed*0xF108E4CD87654321UL; 
        oneway =1UL; 
    } 

    feedback =(rotor1<<1)^(rotor0>>63)^(rotor0); 
    rotor1 =rotor0; 
    rotor0 =feedback;
    
    return((rotor0*0x1525374657E5F50DUL)+0x8B459B95879A07F3UL);
}

int16_t IVA[4][4] = {
  {1306, 2847, 2865, 2852},
  {3259,  144,  125, 1694},
  {1372, 1731, 2033,  320},
  { 722, 1936, 1635, 1792}
};

int16_t VA[4] = {1, 2970, 2379, 1492};

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
    int64_t temp = (x * y);
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

void mod4(int16_t out[8],int16_t in[8])
{
    out[7] = sub_reduce(in[3], in[7]);
    out[6] = sub_reduce(in[2], in[6]);
    out[5] = sub_reduce(in[1], in[5]);
    out[4] = sub_reduce(in[0], in[4]);

    out[3] = add_reduce(in[3], in[7]);
    out[2] = add_reduce(in[2], in[6]);
    out[1] = add_reduce(in[1], in[5]);
    out[0] = add_reduce(in[0], in[4]);
}

//modulo x^2-1: a3X^3+a2X^2+a1X+a0 = (a3X+a2)(x^2-1)+ (a1+a3)X+(a0+a2)
//modulo x^2+1: a3X^3+a2X^2+a1X+a0 = (a3X+a2)(x^2+1)+ (a1-a3)X+(a0-a2)
//modulo x^2+2620∗x+1: a3X^3+a2X^2+a1X+a0 = (a3X+a2-2620*a3)(x^2+2620∗x+1)+(a1-2620*a2+a3)x+a0-a2+2620*a3  2620^2-1=2
//modulo x^2+709∗x+1: a3X^3+a2X^2+a1X+a0 = (a3X+a2-2600*a3)(x^2+709∗x+1)+(a1-709*a2+a3)x+a0-a2+709*a3

void mod2(int16_t out[8],int16_t in[8])
{
    out[0] = add_reduce(in[0], in[2]);
    out[1] = add_reduce(in[1], in[3]);
    out[2] = sub_reduce(in[0], in[2]);
    out[3] = sub_reduce(in[1], in[3]);
    out[4] = add_reduce(sub_reduce(in[4], in[6]), fq_mul_reduce(in[7], 2620));
    out[5] = add_reduce(add_reduce(in[5], in[7]), fq_mul_reduce(in[6], 709));
    out[6] = add_reduce(sub_reduce(in[4], in[6]), fq_mul_reduce(in[7], 709));
    out[7] = add_reduce(add_reduce(in[5], in[7]), fq_mul_reduce(in[6], 2620));
}

//mod x-1      = P(ω^0) = P(ω^0)
//mod x+1      = P(ω^4) = P(ω^(-4)) 
//mod x + 1600 = P(ω^2) = P(ω^(-6))
//mod x + 1729 = P(ω^6) = P(ω^(-2))
//mod x + 40   = P(ω^7) = P(ω^(-1))
//mod x + 2580 = P(ω^1) = P(ω^(-7))
//mod x + 749  = P(ω^5) = P(ω^(-2))
//mod x + 3289 = P(ω^3) = P(ω^(-5))
//a1x+a0 = (x+c)a1 + (a0-a1*c)

void mod1(int16_t out[8],int16_t in[8])
{
    out[0] = add_reduce(in[0], in[1]); 
    out[1] = sub_reduce(in[0], in[1]);
    out[2] = sub_reduce(in[2], fq_mul_reduce(in[3], 1600));
    out[3] = sub_reduce(in[2], fq_mul_reduce(in[3], 1729));
    out[4] = sub_reduce(in[4], fq_mul_reduce(in[5], 40  ));
    out[5] = sub_reduce(in[4], fq_mul_reduce(in[5], 2580));
    out[6] = sub_reduce(in[6], fq_mul_reduce(in[7], 749 ));
    out[7] = sub_reduce(in[6], fq_mul_reduce(in[7], 3289));
}


void FFT(int16_t out[8],int16_t in[8])
{
    int16_t v;
    int16_t ini[8];
    memcpy(ini,in,8*sizeof(int16_t));
    mod4(out,ini);
    mod2(ini,out);
    mod1(out,ini);

    v = out[4];
    out[4] = out[1];
    out[1] = out[5];
    out[5] = out[6];
    out[6] = out[3];
    out[3] = out[7];
    out[7] = v;
    out[2] = out[2];
    out[0] = out[0];
}


void IFFT(int16_t out[8],int16_t in[8])
{
    int16_t v;
    int16_t ini[8];
    memcpy(ini,in,8*sizeof(int16_t));
    mod4(out,ini);
    mod2(ini,out);
    mod1(out,ini);

    out[0] = fq_mul_reduce(out[0], 2913);
    v = fq_mul_reduce(out[1], 2913);
    out[1] = fq_mul_reduce(out[4], 2913);
    out[4] = v;

    v = fq_mul_reduce(out[2], 2913);
    out[2] = fq_mul_reduce(out[3], 2913);
    out[3] = fq_mul_reduce(out[6], 2913);
    out[6] = v;
  
    v = fq_mul_reduce(out[5], 2913);
    out[5] = fq_mul_reduce(out[7], 2913);
    out[7] = v;
}

// multiplication by inverse of vandermonde(α)
void multv(int16_t out[8], int16_t in[4])
{
    for (int i = 0; i < 4; i++) {
        out[i] = 0;
        for(int j  = 0; j < 4; j++) {
            out[i] = add_reduce(out[i], fq_mul_reduce(IVA[j][i], in[j]));
        }
    }
    memset(&out[4],0,4*sizeof(int16_t));
}


void imultv(int16_t *out, int16_t in[4])
{
    *out = 0;
    for(uint8_t i = 0; i < 4; i++) {
        *out = add_reduce(*out, fq_mul_reduce(in[i], VA[i]));
    }
}


// mask of element in Z/3329Z
void mask_q(maskq_t *m, int16_t in)
{
    int16_t vo[8],vi[4];
    uint32_t s1 = rand() & 0xffffffff;
    uint32_t s2 = rand() & 0xffffffff;
    uint16_t r1 = barrett_reduce(s1 & 0xffff);
    uint16_t r2 = barrett_reduce((s1 >> 16) & 0xffff);
    uint16_t r3 = barrett_reduce(s2 & 0xffff);

    vi[0] = in;
    vi[1] = r1;
    vi[2] = r2;
    vi[3] = r3;
  
    multv(vo,vi);
    FFT(m->m,vo); 
}

// unmask  of element in Z/3329Z^8
void unmask_q(int16_t *out,maskq_t *in)
{
    int16_t vo[8];
    IFFT(vo, in->m);
    imultv(out, vo);
}

void affect_q(maskq_t *out, maskq_t *in)
{
    for(int i = 0; i < 8; i++)
        out->m[i] = in->m[i];
}

void op_add_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    for (int i=0; i < 8; i++) {
        rep->m[i] = add_reduce(x1->m[i], x2->m[i]);
    }
}

void op_sub_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    for (int i=0; i < 8; i++)
        rep->m[i] = (3329+x1->m[i] - x2->m[i])%3329;
}

void op_mul_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    maskq_t tmp;
    for(int i = 0; i < 8; i++)
        tmp.m[i] = 0;

    for (uint8_t i = 0; i < 8; i++) {
        tmp.m[i] = add_reduce(tmp.m[i],fq_mul_reduce(x1->m[i], x2->m[i]));
    }
    
    IFFT(rep->m,tmp.m);

    tmp.m[0] = add_reduce(rep->m[0], fq_mul_reduce(rep->m[4], 341));
    tmp.m[1] = add_reduce(rep->m[1], fq_mul_reduce(rep->m[5], 341));
    tmp.m[2] = add_reduce(rep->m[2], fq_mul_reduce(rep->m[6], 341));
    tmp.m[3] = add_reduce(rep->m[3], fq_mul_reduce(rep->m[7], 341));
    tmp.m[4] = 0;
    tmp.m[5] = 0;
    tmp.m[6] = 0;
    tmp.m[7] = 0;

    FFT(rep->m,tmp.m);
}

void op_xor_q(maskq_t *rep, maskq_t *x1, maskq_t *x2)
{
    maskq_t rep1;
    op_mul_q(&rep1, x1, x2);
    for(uint8_t i = 0; i < 8; i++) {
        rep1.m[i] = fq_mul_reduce(3327, rep1.m[i]);
    }

    op_add_q(rep, &rep1, x1);
    op_add_q(rep, rep, x2);
}

void op_xor5_q(maskq_t *rep, maskq_t *x1, maskq_t *x2, maskq_t *x3, maskq_t *x4, maskq_t *x5)
{
    maskq_t rep1, rep2;
    op_xor_q(&rep1, x1, x2);
    op_xor_q(&rep2, x3, x4);
    op_xor_q(rep, &rep1, &rep2);
    op_xor_q(rep, rep, x5);
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
    for(uint8_t i = 0; i < 8; i++) {
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

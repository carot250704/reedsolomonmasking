#include "conversion.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


uint8_t Vec_b2[6][7] = {
    {0,0,1,0,1,0,1},
    {0,1,1,0,0,0,1},
    {0,1,1,0,0,0,1},
    {1,0,0,0,0,1,0},
    {0,1,1,1,0,0,1},
    {0,1,1,1,0,1,1}
};

//convert GFq--->GF(2^6)
void convertq_to_2n(mask2n_t *out, maskq_t *in){
    mask2n_t F;
    int16_t u,i,j;
    maskq_t Rout[6][7];

    mask_2n(&F, 0);
    
    for (i = 0; i < 7; i++)
        for (j = 0; j < 6; j++)
            mask_q(&Rout[j][i],(int16_t)((F.m[i]>>j)&1));

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 7; j++){
            if (Vec_b2[i][j] == 1){
                op_xor_q(&Rout[i][j], &Rout[i][j], in);
            }
        }
    }

    for(i = 0; i < 7; i++) {
        out->m[i] = 0;
    }

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 7; j++){
            unmask_q(&u,&Rout[i][j]);
            out->m[j]^= (uint8_t)(u<<i);
        } 
    }
}


uint8_t Vec_bq[12][8] = {
    {0,1,0,1,1,1,0,1},
    {0,1,0,0,0,1,0,0},
    {1,0,0,0,0,0,1,0},
    {1,1,1,0,1,1,1,0},
    {0,0,0,0,0,0,0,1},
    {0,1,1,0,0,0,1,1},
    {0,1,0,0,0,0,0,1},
    {1,1,0,1,0,1,1,0},
    {0,0,0,1,1,1,1,0},
    {0,1,0,1,1,1,0,1},
    {1,0,1,0,1,1,0,1},
    {1,0,1,0,0,0,0,0}
};


void retenue(mask2n_t *out, mask2n_t *y, mask2n_t *z){
    mask2n_t v, u;
    op_xor_2n(&u,y,z);
    op_and_2n(&v,out,&u);
    op_and_2n(&u,y,z);
    op_xor_2n(out,&u,&v);
}


//convert GF(2^6)--->GFq
void convert2n_to_q(maskq_t *out, mask2n_t *in){
    int32_t i,j;
    maskq_t F;
    mask2n_t Rout[13][8], ret, vet;
    uint8_t u;

    mask_q(&F, 0);   
    for (i = 0; i < 8; i++)
        for (j = 0; j < 12; j++)
            mask_2n(&Rout[j][i], (uint8_t)((F.m[i]>>j)&1));

    for (i = 0; i < 8; i++){
        mask_2n(&ret, 0);
        mask_2n(&vet, 0);
        for (j = 0; j < 12; j++) {
            if (Vec_bq[j][i] == 1){
                retenue(&vet, in, &Rout[j][i]);
                op_xor3_2n(&Rout[j][i], in, &Rout[j][i], &ret);
                affect_2n(&ret, &vet);
            }
            else{
                op_and_2n(&vet, &ret, &Rout[j][i]);
                op_xor_2n(&Rout[j][i], &Rout[j][i], &ret);
                affect_2n(&ret, &vet);
            }
        }
        affect_2n(&Rout[12][i],&ret);
    }

    for(i = 0; i < 8; i++) {
        out->m[i] = 0;
    }

    for (i = 0; i < 13; i++) {
        for (j = 0; j < 8; j++){
            unmask_2n(&u,&Rout[i][j]);
            out->m[j]^= ((int32_t) u)<<i;
        }
    }

    for (j = 0; j < 8; j++)
        out->m[j] = (out->m[j] + 3329) % 3329;  
}



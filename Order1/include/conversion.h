#ifndef CONVERSION_H
#define CONVERSION_H

#include "mask2n.h"
#include "maskq.h"

void convert2n_to_q(maskq_t *out, mask2n_t *in);
void convertq_to_2n(mask2n_t *out, maskq_t *in);
#endif
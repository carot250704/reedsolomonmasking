#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "conversion.h"
#include "mask2n.h"
#include "maskq.h"

#define NROUND 100000

typedef int (*f_test)();

int test_unit_mask2n()
{
    uint8_t a, b, a_xor_b, a_and_b;
    mask2n_t m_a, m_b, m_a_xor_b, m_a_and_b;

    a = rand() & 0x0F;
    b = rand() & 0x0F;

    mask_2n(&m_a, a);
    mask_2n(&m_b, b);

    op_xor_2n(&m_a_xor_b, &m_a, &m_b);
    op_and_2n(&m_a_and_b, &m_a, &m_b);

    unmask_2n(&a_xor_b, &m_a_xor_b);
    unmask_2n(&a_and_b, &m_a_and_b);

    if((a ^ b) != a_xor_b) return 0;
    if((mul2n(a, b) != a_and_b)) return 0;

    return 1;
}

int test_unit_maskq()
{
    int16_t a, b, a_add_b, a_mul_b;
    maskq_t m_a, m_b, m_a_add_b, m_a_mul_b;

    a = rand() % 3329;
    b = rand() % 3329;

    mask_q(&m_a, a);
    mask_q(&m_b, b);

    int16_t aa, bb;
    unmask_q(&aa, &m_a);
    unmask_q(&bb, &m_b);

    op_add_q(&m_a_add_b, &m_a, &m_b);
    op_mul_q(&m_a_mul_b, &m_a, &m_b);

    unmask_q(&a_add_b, &m_a_add_b);
    unmask_q(&a_mul_b, &m_a_mul_b);

    if(add_reduce(a, b) != (a_add_b % 3329) ) return 0;
    if(fq_mul_reduce(a, b) != (a_mul_b % 3329) ) return 0;

    return 1;
}

int test_unit_conversion2n_q()
{
    uint8_t in;
    int16_t out;

    mask2n_t m_in;
    maskq_t m_out;

    in = rand() & 1;
    mask_2n(&m_in, in);

    convert2n_to_q(&m_out, &m_in);

    unmask_q(&out, &m_out);

    if(in != out) return 0;

    return 1;
}

int test_unit_conversionq_2n()
{
    uint16_t in;
    uint8_t out;

    maskq_t m_in;
    mask2n_t m_out;

    in = rand() & 1;
    mask_q(&m_in, in);

    convertq_to_2n(&m_out, &m_in);

    unmask_2n(&out, &m_out);

    if(in != out) return 0;

    return 1;
}

void multiple_test(f_test f, int n_round, char *name)
{
    for(int i = 0; i < n_round; i++) {
        if(!f()) {
            printf("\rTest %s [%02d%%] : ERROR", name, i * 100 / n_round);
            exit(EXIT_FAILURE);
        }
        printf("\rTest %s [%02d%%]", name, i * 100 / n_round);
        fflush(stdout);
    }

    printf("\rTest %s [%02d%%] : SUCCESS\n", name, 100);
}


int main(void)
{
    srand(time(NULL));

    multiple_test(test_unit_mask2n, NROUND, "Gadget Mask F16     ");
    multiple_test(test_unit_maskq, NROUND, "Gadget Mask Fq      ");
    multiple_test(test_unit_conversion2n_q, NROUND, "Conversion F16 -> Fq");
    multiple_test(test_unit_conversion2n_q, NROUND, "Conversion Fq -> F16");

    return 0;
}
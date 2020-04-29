/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_INIT_HPP_
#define BN382_INIT_HPP_
#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>

namespace libff {

const mp_size_t bn382_r_bitcount = 382;
const mp_size_t bn382_q_bitcount = 382;

const mp_size_t bn382_r_limbs = (bn382_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn382_q_limbs = (bn382_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn382_r_limbs> bn382_modulus_r;
extern bigint<bn382_q_limbs> bn382_modulus_q;

extern bn::Fp bn382_coeff_b;
extern size_t bn382_Fq_s;
extern bn::Fp bn382_Fq_nqr_to_t;
extern mie::Vuint bn382_Fq_t_minus_1_over_2;

extern bn::Fp2 bn382_twist_coeff_b;
extern size_t bn382_Fq2_s;
extern bn::Fp2 bn382_Fq2_nqr_to_t;
extern mie::Vuint bn382_Fq2_t_minus_1_over_2;

typedef Fp_model<bn382_r_limbs, bn382_modulus_r> bn382_Fr;
typedef Fp_model<bn382_q_limbs, bn382_modulus_q> bn382_Fq;

void init_bn382_params();

class bn382_G1;
class bn382_G2;
class bn382_GT;
typedef bn382_GT bn382_Fq12;

} // libff
#endif // BN382_INIT_HPP_

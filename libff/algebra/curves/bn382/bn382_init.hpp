/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_INIT_HPP_
#define BN382_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff {

const mp_size_t bn382_r_bitcount = 382;
const mp_size_t bn382_q_bitcount = 382;

const mp_size_t bn382_r_limbs = (bn382_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn382_q_limbs = (bn382_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn382_r_limbs> bn382_modulus_r;
extern bigint<bn382_q_limbs> bn382_modulus_q;

typedef Fp_model<bn382_r_limbs, bn382_modulus_r> bn382_Fr;
typedef Fp_model<bn382_q_limbs, bn382_modulus_q> bn382_Fq;
typedef Fp2_model<bn382_q_limbs, bn382_modulus_q> bn382_Fq2;
typedef Fp6_3over2_model<bn382_q_limbs, bn382_modulus_q> bn382_Fq6;
typedef Fp12_2over3over2_model<bn382_q_limbs, bn382_modulus_q> bn382_Fq12;
typedef bn382_Fq12 bn382_GT;

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern bn382_Fq bn382_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern bn382_Fq2 bn382_twist;
extern bn382_Fq2 bn382_twist_coeff_b;
extern bn382_Fq bn382_twist_mul_by_b_c0;
extern bn382_Fq bn382_twist_mul_by_b_c1;
extern bn382_Fq2 bn382_twist_mul_by_q_X;
extern bn382_Fq2 bn382_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bn382_q_limbs> bn382_ate_loop_count;
extern bool bn382_ate_is_loop_count_neg;
extern bigint<12*bn382_q_limbs> bn382_final_exponent;
extern bigint<bn382_q_limbs> bn382_final_exponent_z;
extern bool bn382_final_exponent_is_z_neg;

void init_bn382_params();

class bn382_G1;
class bn382_G2;

} // libff
#endif // BN382_INIT_HPP_

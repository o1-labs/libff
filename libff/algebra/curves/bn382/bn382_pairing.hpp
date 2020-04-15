/** @file
 ********************************************************************************
 Declares functions for computing Ate pairings over the bn382 curves, split into a
 offline and online stages.
 ********************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *******************************************************************************/

#ifndef BN382_PAIRING_HPP_
#define BN382_PAIRING_HPP_
#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/curves/bn382/bn382_g1.hpp>
#include <libff/algebra/curves/bn382/bn382_g2.hpp>
#include <libff/algebra/curves/bn382/bn382_gt.hpp>

namespace libff {

struct bn382_ate_G1_precomp {
    bn::Fp P[3];

    bool operator==(const bn382_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn382_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bn382_ate_G1_precomp &prec_P);
};

typedef bn::Fp6 bn382_ate_ell_coeffs;

struct bn382_ate_G2_precomp {
    bn::Fp2 Q[3];
    std::vector<bn382_ate_ell_coeffs> coeffs;

    bool operator==(const bn382_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn382_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bn382_ate_G2_precomp &prec_Q);
};

bn382_ate_G1_precomp bn382_ate_precompute_G1(const bn382_G1& P);
bn382_ate_G2_precomp bn382_ate_precompute_G2(const bn382_G2& Q);

bn382_Fq12 bn382_double_ate_miller_loop(const bn382_ate_G1_precomp &prec_P1,
                                        const bn382_ate_G2_precomp &prec_Q1,
                                        const bn382_ate_G1_precomp &prec_P2,
                                        const bn382_ate_G2_precomp &prec_Q2);
bn382_Fq12 bn382_ate_miller_loop(const bn382_ate_G1_precomp &prec_P,
                                 const bn382_ate_G2_precomp &prec_Q);

bn382_GT bn382_final_exponentiation(const bn382_Fq12 &elt);

} // libff
#endif // BN382_PAIRING_HPP_

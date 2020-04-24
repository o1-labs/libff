/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_PAIRING_HPP_
#define BN382_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bn382/bn382_init.hpp>

namespace libff {

/* final exponentiation */

bn382_GT bn382_final_exponentiation(const bn382_Fq12 &elt);

/* ate pairing */

struct bn382_ate_G1_precomp {
    bn382_Fq PX;
    bn382_Fq PY;

    bool operator==(const bn382_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn382_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bn382_ate_G1_precomp &prec_P);
};

struct bn382_ate_ell_coeffs {
    bn382_Fq2 ell_0;
    bn382_Fq2 ell_VW;
    bn382_Fq2 ell_VV;

    bool operator==(const bn382_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn382_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, bn382_ate_ell_coeffs &dc);
};

struct bn382_ate_G2_precomp {
    bn382_Fq2 QX;
    bn382_Fq2 QY;
    std::vector<bn382_ate_ell_coeffs> coeffs;

    bool operator==(const bn382_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn382_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bn382_ate_G2_precomp &prec_Q);
};

bn382_ate_G1_precomp bn382_ate_precompute_G1(const bn382_G1& P);
bn382_ate_G2_precomp bn382_ate_precompute_G2(const bn382_G2& Q);

bn382_Fq12 bn382_ate_miller_loop(const bn382_ate_G1_precomp &prec_P,
                              const bn382_ate_G2_precomp &prec_Q);
bn382_Fq12 bn382_ate_double_miller_loop(const bn382_ate_G1_precomp &prec_P1,
                                     const bn382_ate_G2_precomp &prec_Q1,
                                     const bn382_ate_G1_precomp &prec_P2,
                                     const bn382_ate_G2_precomp &prec_Q2);

bn382_Fq12 bn382_ate_pairing(const bn382_G1& P,
                          const bn382_G2 &Q);
bn382_GT bn382_ate_reduced_pairing(const bn382_G1 &P,
                                 const bn382_G2 &Q);

/* choice of pairing */

typedef bn382_ate_G1_precomp bn382_G1_precomp;
typedef bn382_ate_G2_precomp bn382_G2_precomp;

bn382_G1_precomp bn382_precompute_G1(const bn382_G1& P);

bn382_G2_precomp bn382_precompute_G2(const bn382_G2& Q);

bn382_Fq12 bn382_miller_loop(const bn382_G1_precomp &prec_P,
                          const bn382_G2_precomp &prec_Q);

bn382_Fq12 bn382_double_miller_loop(const bn382_G1_precomp &prec_P1,
                                 const bn382_G2_precomp &prec_Q1,
                                 const bn382_G1_precomp &prec_P2,
                                 const bn382_G2_precomp &prec_Q2);

bn382_Fq12 bn382_pairing(const bn382_G1& P,
                      const bn382_G2 &Q);

bn382_GT bn382_reduced_pairing(const bn382_G1 &P,
                             const bn382_G2 &Q);

bn382_GT bn382_affine_reduced_pairing(const bn382_G1 &P,
                                    const bn382_G2 &Q);

} // libff
#endif // BN382_PAIRING_HPP_

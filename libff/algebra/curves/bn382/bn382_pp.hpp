/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BN382_PP_HPP_
#define BN382_PP_HPP_
#include <libff/algebra/curves/bn382/bn382_g1.hpp>
#include <libff/algebra/curves/bn382/bn382_g2.hpp>
#include <libff/algebra/curves/bn382/bn382_init.hpp>
#include <libff/algebra/curves/bn382/bn382_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bn382_pp {
public:
    typedef bn382_Fr Fp_type;
    typedef bn382_G1 G1_type;
    typedef bn382_G2 G2_type;
    typedef bn382_G1_precomp G1_precomp_type;
    typedef bn382_G2_precomp G2_precomp_type;
    typedef bn382_Fq Fq_type;
    typedef bn382_Fq2 Fqe_type;
    typedef bn382_Fq12 Fqk_type;
    typedef bn382_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bn382_GT final_exponentiation(const bn382_Fq12 &elt);
    static bn382_G1_precomp precompute_G1(const bn382_G1 &P);
    static bn382_G2_precomp precompute_G2(const bn382_G2 &Q);
    static bn382_Fq12 miller_loop(const bn382_G1_precomp &prec_P,
                                      const bn382_G2_precomp &prec_Q);
    static bn382_Fq12 double_miller_loop(const bn382_G1_precomp &prec_P1,
                                             const bn382_G2_precomp &prec_Q1,
                                             const bn382_G1_precomp &prec_P2,
                                             const bn382_G2_precomp &prec_Q2);
    static bn382_Fq12 pairing(const bn382_G1 &P,
                                  const bn382_G2 &Q);
    static bn382_Fq12 reduced_pairing(const bn382_G1 &P,
                                          const bn382_G2 &Q);
};

} // libff

#endif // BN382_PP_HPP_

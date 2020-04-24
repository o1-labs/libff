/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn382/bn382_pp.hpp>

namespace libff {

void bn382_pp::init_public_params()
{
    init_bn382_params();
}

bn382_GT bn382_pp::final_exponentiation(const bn382_Fq12 &elt)
{
    return bn382_final_exponentiation(elt);
}

bn382_G1_precomp bn382_pp::precompute_G1(const bn382_G1 &P)
{
    return bn382_precompute_G1(P);
}

bn382_G2_precomp bn382_pp::precompute_G2(const bn382_G2 &Q)
{
    return bn382_precompute_G2(Q);
}

bn382_Fq12 bn382_pp::miller_loop(const bn382_G1_precomp &prec_P,
                                         const bn382_G2_precomp &prec_Q)
{
    return bn382_miller_loop(prec_P, prec_Q);
}

bn382_Fq12 bn382_pp::double_miller_loop(const bn382_G1_precomp &prec_P1,
                                                const bn382_G2_precomp &prec_Q1,
                                                const bn382_G1_precomp &prec_P2,
                                                const bn382_G2_precomp &prec_Q2)
{
    return bn382_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bn382_Fq12 bn382_pp::pairing(const bn382_G1 &P,
                                     const bn382_G2 &Q)
{
    return bn382_pairing(P, Q);
}

bn382_Fq12 bn382_pp::reduced_pairing(const bn382_G1 &P,
                                             const bn382_G2 &Q)
{
    return bn382_reduced_pairing(P, Q);
}

} // libff

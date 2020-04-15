/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn382/bn382_g1.hpp>
#include <libff/algebra/curves/bn382/bn382_g2.hpp>
#include <libff/algebra/curves/bn382/bn382_gt.hpp>
#include <libff/algebra/curves/bn382/bn382_init.hpp>

namespace libff {

bigint<bn382_r_limbs> bn382_modulus_r;
bigint<bn382_q_limbs> bn382_modulus_q;

bn::Fp bn382_coeff_b;
size_t bn382_Fq_s;
bn::Fp bn382_Fq_nqr_to_t;
mie::Vuint bn382_Fq_t_minus_1_over_2;

bn::Fp2 bn382_twist_coeff_b;
size_t bn382_Fq2_s;
bn::Fp2 bn382_Fq2_nqr_to_t;
mie::Vuint bn382_Fq2_t_minus_1_over_2;

void init_bn382_params()
{
    bn::Param::init(); // init ate-pairing library

    typedef bigint<bn382_r_limbs> bigint_r;
    typedef bigint<bn382_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */
    bn382_modulus_r = bigint_r("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    assert(bn382_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn382_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        bn382_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        bn382_Fr::inv = 0xc2e1f593efffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn382_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        bn382_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        bn382_Fr::inv = 0xefffffff;
    }
    bn382_Fr::num_bits = 254;
    bn382_Fr::euler = bigint_r("10944121435919637611123202872628637544274182200208017171849102093287904247808");
    bn382_Fr::s = 28;
    bn382_Fr::t = bigint_r("81540058820840996586704275553141814055101440848469862132140264610111");
    bn382_Fr::t_minus_1_over_2 = bigint_r("40770029410420498293352137776570907027550720424234931066070132305055");
    bn382_Fr::multiplicative_generator = bn382_Fr("5");
    bn382_Fr::root_of_unity = bn382_Fr("19103219067921713944291392827692070036145651957329286315305642004821462161904");
    bn382_Fr::nqr = bn382_Fr("5");
    bn382_Fr::nqr_to_t = bn382_Fr("19103219067921713944291392827692070036145651957329286315305642004821462161904");

    /* parameters for base field Fq */
    bn382_modulus_q = bigint_q("21888242871839275222246405745257275088696311157297823662689037894645226208583");
    assert(bn382_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn382_Fq::Rsquared = bigint_q("3096616502983703923843567936837374451735540968419076528771170197431451843209");
        bn382_Fq::Rcubed = bigint_q("14921786541159648185948152738563080959093619838510245177710943249661917737183");
        bn382_Fq::inv = 0x87d20782e4866389;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn382_Fq::Rsquared = bigint_q("3096616502983703923843567936837374451735540968419076528771170197431451843209");
        bn382_Fq::Rcubed = bigint_q("14921786541159648185948152738563080959093619838510245177710943249661917737183");
        bn382_Fq::inv = 0xe4866389;
    }
    bn382_Fq::num_bits = 254;
    bn382_Fq::euler = bigint_q("10944121435919637611123202872628637544348155578648911831344518947322613104291");
    bn382_Fq::s = 1;
    bn382_Fq::t = bigint_q("10944121435919637611123202872628637544348155578648911831344518947322613104291");
    bn382_Fq::t_minus_1_over_2 = bigint_q("5472060717959818805561601436314318772174077789324455915672259473661306552145");
    bn382_Fq::multiplicative_generator = bn382_Fq("3");
    bn382_Fq::root_of_unity = bn382_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");
    bn382_Fq::nqr = bn382_Fq("3");
    bn382_Fq::nqr_to_t = bn382_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");

    /* additional parameters for square roots in Fq/Fq2 */
    bn382_coeff_b = bn::Fp(3);
    bn382_Fq_s = 1;
    bn382_Fq_nqr_to_t = bn::Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582");
    bn382_Fq_t_minus_1_over_2 = mie::Vuint("5472060717959818805561601436314318772174077789324455915672259473661306552145");

    bn382_twist_coeff_b = bn::Fp2(bn::Fp("19485874751759354771024239261021720505790618469301721065564631296452457478373"),
                                  bn::Fp("266929791119991161246907387137283842545076965332900288569378510910307636690"));
    bn382_Fq2_s = 4;
    bn382_Fq2_nqr_to_t = bn::Fp2(bn::Fp("5033503716262624267312492558379982687175200734934877598599011485707452665730"),
                                 bn::Fp("314498342015008975724433667930697407966947188435857772134235984660852259084"));
    bn382_Fq2_t_minus_1_over_2 = mie::Vuint("14971724250519463826312126413021210649976634891596900701138993820439690427699319920245032869357433499099632259837909383182382988566862092145199781964621");

    /* choice of group G1 */
    bn382_G1::G1_zero.coord[0] = bn::Fp(1);
    bn382_G1::G1_zero.coord[1] = bn::Fp(1);
    bn382_G1::G1_zero.coord[2] = bn::Fp(0);

    bn382_G1::G1_one.coord[0] = bn::Fp(1);
    bn382_G1::G1_one.coord[1] = bn::Fp(2);
    bn382_G1::G1_one.coord[2] = bn::Fp(1);

    bn382_G1::wnaf_window_table.resize(0);
    bn382_G1::wnaf_window_table.push_back(10);
    bn382_G1::wnaf_window_table.push_back(24);
    bn382_G1::wnaf_window_table.push_back(40);
    bn382_G1::wnaf_window_table.push_back(132);

    bn382_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.24]
    bn382_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.24, 10.43]
    bn382_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [10.43, 24.88]
    bn382_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [24.88, 62.10]
    bn382_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [62.10, 157.80]
    bn382_G1::fixed_base_exp_window_table.push_back(62);
    // window 6 is unbeaten in [157.80, 362.05]
    bn382_G1::fixed_base_exp_window_table.push_back(158);
    // window 7 is unbeaten in [362.05, 806.67]
    bn382_G1::fixed_base_exp_window_table.push_back(362);
    // window 8 is unbeaten in [806.67, 2090.34]
    bn382_G1::fixed_base_exp_window_table.push_back(807);
    // window 9 is unbeaten in [2090.34, 4459.58]
    bn382_G1::fixed_base_exp_window_table.push_back(2090);
    // window 10 is unbeaten in [4459.58, 9280.12]
    bn382_G1::fixed_base_exp_window_table.push_back(4460);
    // window 11 is unbeaten in [9280.12, 43302.64]
    bn382_G1::fixed_base_exp_window_table.push_back(9280);
    // window 12 is unbeaten in [43302.64, 210998.73]
    bn382_G1::fixed_base_exp_window_table.push_back(43303);
    // window 13 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 14 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [210998.73, 506869.47]
    bn382_G1::fixed_base_exp_window_table.push_back(210999);
    // window 16 is unbeaten in [506869.47, 930023.36]
    bn382_G1::fixed_base_exp_window_table.push_back(506869);
    // window 17 is unbeaten in [930023.36, 8350812.20]
    bn382_G1::fixed_base_exp_window_table.push_back(930023);
    // window 18 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [8350812.20, 21708138.87]
    bn382_G1::fixed_base_exp_window_table.push_back(8350812);
    // window 21 is unbeaten in [21708138.87, 29482995.52]
    bn382_G1::fixed_base_exp_window_table.push_back(21708139);
    // window 22 is unbeaten in [29482995.52, inf]
    bn382_G1::fixed_base_exp_window_table.push_back(29482996);

    /* choice of group G2 */
    bn382_G2::G2_zero.coord[0] = bn::Fp2(bn::Fp(1), bn::Fp(0));
    bn382_G2::G2_zero.coord[1] = bn::Fp2(bn::Fp(1), bn::Fp(0));
    bn382_G2::G2_zero.coord[2] = bn::Fp2(bn::Fp(0), bn::Fp(0));

    bn382_G2::G2_one.coord[0] = bn::Fp2(bn::Fp("15267802884793550383558706039165621050290089775961208824303765753922461897946"),
                                        bn::Fp("9034493566019742339402378670461897774509967669562610788113215988055021632533"));
    bn382_G2::G2_one.coord[1] = bn::Fp2(bn::Fp("644888581738283025171396578091639672120333224302184904896215738366765861164"),
                                        bn::Fp("20532875081203448695448744255224543661959516361327385779878476709582931298750"));
    bn382_G2::G2_one.coord[2] = bn::Fp2(bn::Fp(1), bn::Fp(0));

    bn382_G2::wnaf_window_table.resize(0);
    bn382_G2::wnaf_window_table.push_back(7);
    bn382_G2::wnaf_window_table.push_back(18);
    bn382_G2::wnaf_window_table.push_back(35);
    bn382_G2::wnaf_window_table.push_back(116);

    bn382_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.13]
    bn382_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.13, 10.72]
    bn382_G2::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [10.72, 25.60]
    bn382_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.60, 60.99]
    bn382_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.99, 153.66]
    bn382_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [153.66, 353.13]
    bn382_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [353.13, 771.87]
    bn382_G2::fixed_base_exp_window_table.push_back(353);
    // window 8 is unbeaten in [771.87, 2025.85]
    bn382_G2::fixed_base_exp_window_table.push_back(772);
    // window 9 is unbeaten in [2025.85, 4398.65]
    bn382_G2::fixed_base_exp_window_table.push_back(2026);
    // window 10 is unbeaten in [4398.65, 10493.42]
    bn382_G2::fixed_base_exp_window_table.push_back(4399);
    // window 11 is unbeaten in [10493.42, 37054.73]
    bn382_G2::fixed_base_exp_window_table.push_back(10493);
    // window 12 is unbeaten in [37054.73, 49928.78]
    bn382_G2::fixed_base_exp_window_table.push_back(37055);
    // window 13 is unbeaten in [49928.78, 114502.82]
    bn382_G2::fixed_base_exp_window_table.push_back(49929);
    // window 14 is unbeaten in [114502.82, 161445.26]
    bn382_G2::fixed_base_exp_window_table.push_back(114503);
    // window 15 is unbeaten in [161445.26, 470648.01]
    bn382_G2::fixed_base_exp_window_table.push_back(161445);
    // window 16 is unbeaten in [470648.01, 1059821.87]
    bn382_G2::fixed_base_exp_window_table.push_back(470648);
    // window 17 is unbeaten in [1059821.87, 5450848.25]
    bn382_G2::fixed_base_exp_window_table.push_back(1059822);
    // window 18 is never the best
    bn382_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [5450848.25, 5566795.57]
    bn382_G2::fixed_base_exp_window_table.push_back(5450848);
    // window 20 is unbeaten in [5566795.57, 33055217.52]
    bn382_G2::fixed_base_exp_window_table.push_back(5566796);
    // window 21 is never the best
    bn382_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [33055217.52, inf]
    bn382_G2::fixed_base_exp_window_table.push_back(33055218);

    bn382_GT::GT_one.elem = bn::Fp12(1);
}
} // libff

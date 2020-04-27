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
    bn382_modulus_r = bigint_r("5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289");
    assert(bn382_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn382_Fr::Rsquared = bigint_r("3692793417111842983091090013000149450637488604708422183249448572043162112253810843515561570513520991143563448829697");
        bn382_Fr::Rcubed = bigint_r("4913232160216665594556091511630695940892941228853749698866551234284029208211897326645383305263534916462281663165078");
        bn382_Fr::inv = 0xffffffffffffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn382_Fr::Rsquared = bigint_r("3692793417111842983091090013000149450637488604708422183249448572043162112253810843515561570513520991143563448829697");
        bn382_Fr::Rcubed = bigint_r("4913232160216665594556091511630695940892941228853749698866551234284029208211897326645383305263534916462281663165078");
        bn382_Fr::inv = 0xffffffffffffffff;
    }
    bn382_Fr::num_bits = 382;
    bn382_Fr::euler = bigint_r("3894323663264254867678158639472658262036772463439134629391769872854144375326483567660487859589139148692278992014156");
    bn382_Fr::s = 67;
    bn382_Fr::t = bigint_r("4428027988584290732738588662791751599513670173863245544187769873421953922449150291326441323012538034144073627222542");
    bn382_Fr::t_minus_1_over_2 = bigint_r("564703292445634606241957248126406675950231116336936597893086061646769357651770221740930140567000978695373845933138");
    bn382_Fr::multiplicative_generator = bn382_Fr("7");
    bn382_Fr::root_of_unity = bn382_Fr("347738786800231438048736028385710683070191513687154426015823034492371540518418305403341320400316253561516694398836");
    bn382_Fr::nqr = bn382_Fr("7"); // TODO: double check
    bn382_Fr::nqr_to_t = bn382_Fr("2184517633355941102883278064663814533493735786041771354571793524338655829174535545409423981266833306361135553065740");

    /* parameters for base field Fq */
    bn382_modulus_q = bigint_q("5543634365110765627805495722742127385843376434033820803590214255538854698464778703795540858859767700241957783601153");
    assert(bn382_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn382_Fq::Rsquared = bigint_q("1813506259290382479117279305834205666818874691074009175269296895332276467012488971276719893797369804234162283022471");
        bn382_Fq::Rcubed = bigint_q("4553601541650426936969305648385543158864846427139547883987218649667965527929416945749946961797009740634862502205912");
        bn382_Fq::inv = 0xffffffffffffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn382_Fq::Rsquared = bigint_q("1813506259290382479117279305834205666818874691074009175269296895332276467012488971276719893797369804234162283022471");
        bn382_Fq::Rcubed = bigint_q("4553601541650426936969305648385543158864846427139547883987218649667965527929416945749946961797009740634862502205912");
        bn382_Fq::inv = 0xffffffffffffffff;
    }
    bn382_Fq::num_bits = 382
    bn382_Fq::euler = bigint_q("2678596848750647316014813728257168746524174408889950683037738455054638950116074495731045092587204569872971777290803");
    bn382_Fq::s = 67;
    bn382_Fq::t = bigint_q("5465934171272772363783567931435954309364188005437100548925363201014921857064496612009583101810160489789900286844113");
    bn382_Fq::t_minus_1_over_2 = bigint_q("2639746751831650684003849832604082208284580194591590555705312927792672529415933449838066214062400964646943028912283");
    bn382_Fq::multiplicative_generator = bn382_Fq("7");
    bn382_Fq::root_of_unity = bn382_Fq("4297662006852145359943630658140775691253257231821187487981932496787264968852562613474710924795125596632514440450270");
    bn382_Fq::nqr = bn382_Fq("7"); // TODO: double check
    bn382_Fq::nqr_to_t = bn382_Fq("926642806721114476368757157645289729546034110359000734340680542302430854664018239532699015859941333394840267162320"); // TODO: double check

    /* additional parameters for square roots in Fq/Fq2 */
    bn382_coeff_b = bn::Fp(14); // B of G1
    bn382_Fq_s = 67;
    bn382_Fq_nqr_to_t = bn382_Fq("926642806721114476368757157645289729546034110359000734340680542302430854664018239532699015859941333394840267162320"); // TODO: double check
    bn382_Fq_t_minus_1_over_2 = mie::Vuint("2639746751831650684003849832604082208284580194591590555705312927792672529415933449838066214062400964646943028912283");

    bn382_twist_coeff_b = bn::Fp2(bn::Fp("0"),
                                  bn::Fp("671741409037656549287655731709824109253980562797465531047568917158473772953357661607607074171171789249425365013734");
    bn382_Fq2_s = 4;                                                                                                        // TODO
    bn382_Fq2_nqr_to_t = bn::Fp2(bn::Fp("5033503716262624267312492558379982687175200734934877598599011485707452665730"),    // TODO
                                 bn::Fp("314498342015008975724433667930697407966947188435857772134235984660852259084"));    // TODO
    bn382_Fq2_t_minus_1_over_2 = mie::Vuint("14971724250519463826312126413021210649976634891596900701138993820439690427699319920245032869357433499099632259837909383182382988566862092145199781964621"); // TODO

    /* choice of group G1 */
    bn382_G1::G1_zero.coord[0] = bn::Fp(1);
    bn382_G1::G1_zero.coord[1] = bn::Fp(1);
    bn382_G1::G1_zero.coord[2] = bn::Fp(0);

    bn382_G1::G1_one.coord[0] = bn::Fp(1);
    bn382_G1::G1_one.coord[1] = bn::Fp("93360544046129830094757569027791679210844519762232758194920967606984287664392872848607365449491441272860487554919");
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

    bn382_G2::G2_one.coord[0] = bn::Fp2(bn::Fp("5091479006341624589567896397635435258574014748076809289641574502625108749078943401554928186045022840715545119724980"),
                                        bn::Fp("3519382844713541579002133617775236000337302709092053889907196608497211512910083011998063983635946531824025900302318"));
    bn382_G2::G2_one.coord[1] = bn::Fp2(bn::Fp("3934462613855637686263305666415197064493526818650772586512345121679314757894509046665527945441022114959626478116310"),
                                        bn::Fp("4780208203968490754926961189497985186234872265999339695883768193648752722495923801951744797689497641942946290071424"));
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

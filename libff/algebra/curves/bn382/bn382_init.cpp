/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn382/bn382_g1.hpp>
#include <libff/algebra/curves/bn382/bn382_g2.hpp>
#include <libff/algebra/curves/bn382/bn382_init.hpp>

namespace libff {

bigint<bn382_r_limbs> bn382_modulus_r;
bigint<bn382_q_limbs> bn382_modulus_q;

bn382_Fq bn382_coeff_b;
bn382_Fq2 bn382_twist;
bn382_Fq2 bn382_twist_coeff_b;
bn382_Fq bn382_twist_mul_by_b_c0;
bn382_Fq bn382_twist_mul_by_b_c1;
bn382_Fq2 bn382_twist_mul_by_q_X;
bn382_Fq2 bn382_twist_mul_by_q_Y;

bigint<bn382_q_limbs> bn382_ate_loop_count;
bool bn382_ate_is_loop_count_neg;
bigint<12*bn382_q_limbs> bn382_final_exponent;
bigint<bn382_q_limbs> bn382_final_exponent_z;
bool bn382_final_exponent_is_z_neg;

void init_bn382_params()
{
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

    /* parameters for twist field Fq2 */
    bn382_Fq2::euler = bigint<2*bn382_q_limbs>("239547588008311421220994022608339370399626158265550411218223901127035046843189118723920525909718935985594116157406550130918127817069793474323196511433944");
    bn382_Fq2::s = 4;
    bn382_Fq2::t = bigint<2*bn382_q_limbs>("29943448501038927652624252826042421299953269783193801402277987640879380855398639840490065738714866998199264519675818766364765977133724184290399563929243");
    bn382_Fq2::t_minus_1_over_2 = bigint<2*bn382_q_limbs>("14971724250519463826312126413021210649976634891596900701138993820439690427699319920245032869357433499099632259837909383182382988566862092145199781964621");
    bn382_Fq2::non_residue = bn382_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");
    bn382_Fq2::nqr = bn382_Fq2(bn382_Fq("2"),bn382_Fq("1"));
    bn382_Fq2::nqr_to_t = bn382_Fq2(bn382_Fq("5033503716262624267312492558379982687175200734934877598599011485707452665730"),bn382_Fq("314498342015008975724433667930697407966947188435857772134235984660852259084"));
    bn382_Fq2::Frobenius_coeffs_c1[0] = bn382_Fq("1");
    bn382_Fq2::Frobenius_coeffs_c1[1] = bn382_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");

    /* parameters for Fq6 */
    bn382_Fq6::non_residue = bn382_Fq2(bn382_Fq("9"),bn382_Fq("1"));
    bn382_Fq6::Frobenius_coeffs_c1[0] = bn382_Fq2(bn382_Fq("1"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c1[1] = bn382_Fq2(bn382_Fq("21575463638280843010398324269430826099269044274347216827212613867836435027261"),bn382_Fq("10307601595873709700152284273816112264069230130616436755625194854815875713954"));
    bn382_Fq6::Frobenius_coeffs_c1[2] = bn382_Fq2(bn382_Fq("21888242871839275220042445260109153167277707414472061641714758635765020556616"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c1[3] = bn382_Fq2(bn382_Fq("3772000881919853776433695186713858239009073593817195771773381919316419345261"),bn382_Fq("2236595495967245188281701248203181795121068902605861227855261137820944008926"));
    bn382_Fq6::Frobenius_coeffs_c1[4] = bn382_Fq2(bn382_Fq("2203960485148121921418603742825762020974279258880205651966"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c1[5] = bn382_Fq2(bn382_Fq("18429021223477853657660792034369865839114504446431234726392080002137598044644"),bn382_Fq("9344045779998320333812420223237981029506012124075525679208581902008406485703"));
    bn382_Fq6::Frobenius_coeffs_c2[0] = bn382_Fq2(bn382_Fq("1"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c2[1] = bn382_Fq2(bn382_Fq("2581911344467009335267311115468803099551665605076196740867805258568234346338"),bn382_Fq("19937756971775647987995932169929341994314640652964949448313374472400716661030"));
    bn382_Fq6::Frobenius_coeffs_c2[2] = bn382_Fq2(bn382_Fq("2203960485148121921418603742825762020974279258880205651966"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c2[3] = bn382_Fq2(bn382_Fq("5324479202449903542726783395506214481928257762400643279780343368557297135718"),bn382_Fq("16208900380737693084919495127334387981393726419856888799917914180988844123039"));
    bn382_Fq6::Frobenius_coeffs_c2[4] = bn382_Fq2(bn382_Fq("21888242871839275220042445260109153167277707414472061641714758635765020556616"),bn382_Fq("0"));
    bn382_Fq6::Frobenius_coeffs_c2[5] = bn382_Fq2(bn382_Fq("13981852324922362344252311234282257507216387789820983642040889267519694726527"),bn382_Fq("7629828391165209371577384193250820201684255241773809077146787135900891633097"));

    /* parameters for Fq12 */

    bn382_Fq12::non_residue = bn382_Fq2(bn382_Fq("9"),bn382_Fq("1"));
    bn382_Fq12::Frobenius_coeffs_c1[0]  = bn382_Fq2(bn382_Fq("1"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[1]  = bn382_Fq2(bn382_Fq("8376118865763821496583973867626364092589906065868298776909617916018768340080"),bn382_Fq("16469823323077808223889137241176536799009286646108169935659301613961712198316"));
    bn382_Fq12::Frobenius_coeffs_c1[2]  = bn382_Fq2(bn382_Fq("21888242871839275220042445260109153167277707414472061641714758635765020556617"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[3]  = bn382_Fq2(bn382_Fq("11697423496358154304825782922584725312912383441159505038794027105778954184319"),bn382_Fq("303847389135065887422783454877609941456349188919719272345083954437860409601"));
    bn382_Fq12::Frobenius_coeffs_c1[4]  = bn382_Fq2(bn382_Fq("21888242871839275220042445260109153167277707414472061641714758635765020556616"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[5]  = bn382_Fq2(bn382_Fq("3321304630594332808241809054958361220322477375291206261884409189760185844239"),bn382_Fq("5722266937896532885780051958958348231143373700109372999374820235121374419868"));
    bn382_Fq12::Frobenius_coeffs_c1[6]  = bn382_Fq2(bn382_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[7]  = bn382_Fq2(bn382_Fq("13512124006075453725662431877630910996106405091429524885779419978626457868503"),bn382_Fq("5418419548761466998357268504080738289687024511189653727029736280683514010267"));
    bn382_Fq12::Frobenius_coeffs_c1[8]  = bn382_Fq2(bn382_Fq("2203960485148121921418603742825762020974279258880205651966"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[9]  = bn382_Fq2(bn382_Fq("10190819375481120917420622822672549775783927716138318623895010788866272024264"),bn382_Fq("21584395482704209334823622290379665147239961968378104390343953940207365798982"));
    bn382_Fq12::Frobenius_coeffs_c1[10] = bn382_Fq2(bn382_Fq("2203960485148121921418603742825762020974279258880205651967"),bn382_Fq("0"));
    bn382_Fq12::Frobenius_coeffs_c1[11] = bn382_Fq2(bn382_Fq("18566938241244942414004596690298913868373833782006617400804628704885040364344"),bn382_Fq("16165975933942742336466353786298926857552937457188450663314217659523851788715"));

    /* choice of short Weierstrass curve and its twist */

    bn382_coeff_b = bn382_Fq("3");
    bn382_twist = bn382_Fq2(bn382_Fq("9"), bn382_Fq("1"));
    bn382_twist_coeff_b = bn382_coeff_b * bn382_twist.inverse();
    bn382_twist_mul_by_b_c0 = bn382_coeff_b * bn382_Fq2::non_residue;
    bn382_twist_mul_by_b_c1 = bn382_coeff_b * bn382_Fq2::non_residue;
    bn382_twist_mul_by_q_X = bn382_Fq2(bn382_Fq("21575463638280843010398324269430826099269044274347216827212613867836435027261"),
                                           bn382_Fq("10307601595873709700152284273816112264069230130616436755625194854815875713954"));
    bn382_twist_mul_by_q_Y = bn382_Fq2(bn382_Fq("2821565182194536844548159561693502659359617185244120367078079554186484126554"),
                                           bn382_Fq("3505843767911556378687030309984248845540243509899259641013678093033130930403"));

    /* choice of group G1 */
    bn382_G1::G1_zero = bn382_G1(bn382_Fq::zero(),
                                     bn382_Fq::one(),
                                     bn382_Fq::zero());
    bn382_G1::G1_one = bn382_G1(bn382_Fq("1"),
                                    bn382_Fq("2"),
                                    bn382_Fq::one());
    bn382_G1::wnaf_window_table.resize(0);
    bn382_G1::wnaf_window_table.push_back(11);
    bn382_G1::wnaf_window_table.push_back(24);
    bn382_G1::wnaf_window_table.push_back(60);
    bn382_G1::wnaf_window_table.push_back(127);

    bn382_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bn382_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bn382_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bn382_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bn382_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bn382_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bn382_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bn382_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bn382_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bn382_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bn382_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bn382_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bn382_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bn382_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bn382_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bn382_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bn382_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bn382_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bn382_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bn382_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

    bn382_G2::G2_zero = bn382_G2(bn382_Fq2::zero(),
                                     bn382_Fq2::one(),
                                     bn382_Fq2::zero());

    bn382_G2::G2_one = bn382_G2(bn382_Fq2(bn382_Fq("10857046999023057135944570762232829481370756359578518086990519993285655852781"),
                                                bn382_Fq("11559732032986387107991004021392285783925812861821192530917403151452391805634")),
                                    bn382_Fq2(bn382_Fq("8495653923123431417604973247489272438418190587263600148770280649306958101930"),
                                                bn382_Fq("4082367875863433681332203403145435568316851327593401208105741076214120093531")),
                                    bn382_Fq2::one());
    bn382_G2::wnaf_window_table.resize(0);
    bn382_G2::wnaf_window_table.push_back(5);
    bn382_G2::wnaf_window_table.push_back(15);
    bn382_G2::wnaf_window_table.push_back(39);
    bn382_G2::wnaf_window_table.push_back(109);

    bn382_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bn382_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bn382_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bn382_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bn382_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bn382_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bn382_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bn382_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bn382_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bn382_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bn382_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bn382_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bn382_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bn382_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bn382_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bn382_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bn382_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bn382_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bn382_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bn382_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bn382_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bn382_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bn382_G2::fixed_base_exp_window_table.push_back(31673815);

    /* pairing parameters */

    bn382_ate_loop_count = bigint_q("29793968203157093288");
    bn382_ate_is_loop_count_neg = false;
    bn382_final_exponent = bigint<12*bn382_q_limbs>("552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480");
    bn382_final_exponent_z = bigint_q("4965661367192848881");
    bn382_final_exponent_is_z_neg = false;

}
} // libff

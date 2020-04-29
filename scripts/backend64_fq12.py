import math
BN382_p = 5543634365110765627805495722742127385843376434033820803590214255538854698464778703795540858859767700241957783601153
BN382_q = 5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289
# 382 bits = 48 bytes
# field modulus and group order differ only in the 25th - 32nd bytes
MODULUS_NUM = 3 + BN382_q

MODULUS_BITS = 384
LIMB_SIZES = [64, 64, 64, 64, 64, 64]
WORD_SIZE = 64

# Check that MODULUS_BITS is correct
assert(2**MODULUS_BITS > MODULUS_NUM)
# assert(2**(MODULUS_BITS - 1) < MODULUS_NUM)

# Check that limb sizes are correct
tmp = 0
for i in range(0, len(LIMB_SIZES)):
        # assert(LIMB_SIZES[i] < WORD_SIZE)
        tmp += LIMB_SIZES[i]
assert(tmp == MODULUS_BITS)

# Compute the value of the modulus in this representation
def bigint_to_limbs(MODULUS_NUM):
    MODULUS = []
    tmp = MODULUS_NUM
    print("MODULUS = [")
    for i in range(0, len(LIMB_SIZES)):
        this_modulus_num = tmp & ((2**LIMB_SIZES[i]) - 1)
        MODULUS.append(this_modulus_num)
        tmp = tmp >> LIMB_SIZES[i]
        print("\t", hex(this_modulus_num), ",")
    print("]")
    return

inv = 3298621403693021520254674166538938247613207941189372348401597750128415207145609847844581041878536076753325935356266 
def limbs_to_bigint(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    s = s * inv % BN382_q
    return s

MOD_BITS = 382
CAPACITY = MODULUS_BITS - 1
REPR_SHAVE_BITS = 2
INV_p = 18446744073709551615
INV_q = 18446744073709551615
TWO_ADICITY = 67

def print_name(name, num):
    s = str(num)
    print(name + " : " + " " * (140 - (len(s) + len(name))) + s)
    return
    # for Fq12
    # non_residue ( same as fq6?? )

# Frob_coeffs_c1[0]
FROB_COEFFS_C1 = [[0xfffffffffffffff9,
                0xfffffff57fab5757,
                0x7f56ac056aeaf57f,
                0x10388572e3c2c0f5,
                0xe6ce591c2bafc343,
                0x3e03f4104144b1a],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[1]
                [0x51b19a7e3871df8a,
                0xff256c8a6096ca14,
                0x3c5ed207a2e9ac81,
                0xee047eb105d3e89c,
                0x59e5bf1f71597093,
                0x2226c77500bb1b4b],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[2]
                [0x43ac10f69cd0866e,
                0xb67658d4844670fa,
                0x64500aac20e3e056,
                0xe69857d69abfc002,
                0x521ddf42ec5832c5,
                0xee09eba205fe5d8],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[3]
                [0x16b744a7d72fb912,
                0x8db76da14b98776d,
                0xd7d0fda03758326c,
                0x9a05f3af0ce04699,
                0x1c8a66ecb161efb2,
                0x13a9f1d5f1261bfe],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[4]
                [0x43ac10f69cd08675,
                0xb67658df049b19a2,
                0xe4f95ea6b5f8ead6,
                0xd65fd263b6fcff0c,
                0x6b4f8626c0a86f82,
                0xb005f791c4b9abd],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[5]
                [0xc505aa299ebdd989,
                0x8e9201186b0dc570,
                0x1b8a5c2a1771876a,
                0x608bab122fa766ff,
                0x33f4e436f0a63ea7,
                0x1587b3a0cb438841],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[6]
                [0x8,
                0xc0060c0c0,
                0xc1848c18180c00,
                0xa451b0a144d8480c,
                0x8a81e34d84edfc45,
                0x202449fed6c43c73],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[7]
                [0xae4e6581c78e2077,
                0xda93771f754e03,
                0x43b95e89e01954fe,
                0xc685b76322c72065,
                0x176a7d4a3f444ef4,
                0x1ddc1cada1d6c43],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[8]
                [0xbc53ef09632f7993,
                0x4989a72cfbc5a71d,
                0x1bc825e5621f2129,
                0xcdf1de3d8ddb48ff,
                0x1f325d26c4458cc2,
                0x1523ea85ba78a1b6],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[9]
                [0xe948bb5828d046ef,
                0x724892603473a0aa,
                0xa84732f14baacf13,
                0x1a8442651bbac267,
                0x54c5d57cff3bcfd6,
                0x105a9769e9b26b90],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[10]
                [0xbc53ef09632f798c,
                0x4989a7227b70fe75,
                0x9b1ed1eacd0a16a9,
                0xde2a63b0719e09f4,
                0x600b642eff55005,
                0x190429c6be8cecd1],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[11]
                [0x3afa55d661422678,
                0x716dfee914fe52a7,
                0x648dd4676b917a15,
                0x53fe8b01f8f3a202,
                0x3d5b5832bff780e1,
                0xe7cd59f0f94ff4d],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0]]

# choice of group G1

# alt_bn128_G1::G1_zero = alt_bn128_G1(0, 1, 0)
# alt_bn128_G1::G1_one = alt_bn128_G1(alt_bn128_Fq("1"),
                                # alt_bn128_Fq("2"),
                                # alt_bn128_Fq::one());

# choice of group G2

# alt_bn128_G2::G2_zero = alt_bn128_G2(alt_bn128_Fq2::zero(),
                                 # alt_bn128_Fq2::one(),
                                 # alt_bn128_Fq2::zero());

# alt_bn128_G2::G2_one = alt_bn128_G2(alt_bn128_Fq2(alt_bn128_Fq("10857046999023057135944570762232829481370756359578518086990519993285655852781"),
#                                             alt_bn128_Fq("11559732032986387107991004021392285783925812861821192530917403151452391805634")),
#                                 alt_bn128_Fq2(alt_bn128_Fq("8495653923123431417604973247489272438418190587263600148770280649306958101930"),
#                                             alt_bn128_Fq("4082367875863433681332203403145435568316851327593401208105741076214120093531")),
#                                 alt_bn128_Fq2::one());

# pairing parameters 
# alt_bn128_ate_loop_count = bigint_q("29793968203157093288");
# alt_bn128_ate_is_loop_count_neg = false;
# alt_bn128_final_exponent = bigint<12*alt_bn128_q_limbs>("552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480");
# alt_bn128_final_exponent_z = bigint_q("4965661367192848881");
# alt_bn128_final_exponent_is_z_neg = false;
for i, x in enumerate(FROB_COEFFS_C1):
    print_name("frob_coeff " + str(i), limbs_to_bigint(x))

print("_init.hpp needs bit lengths changing")
print("_pairing.cpp needs nothing") 
print("_pairing.hpp needs nothing") 
print("_pp.cpp needs nothing") 
print("_pp.hpp needs nothing") 




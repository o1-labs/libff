import math
# E1/Fp : y^2 = x^3 + 7
BN382_p = 5543634365110765627805495722742127385843376434033820803590214255538854698464778703795540858859767700241957783601153
BN382_q = 5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289
# 382 bits = 48 bytes
# field modulus and group order differ only in the 25th - 32nd bytes

MODULUS_BITS = 384
LIMB_SIZES = [64, 64, 64, 64, 64, 64]
WORD_SIZE = 64

MOD_BITS = 382
CAPACITY = MODULUS_BITS - 1
REPR_SHAVE_BITS = 2
INV_p = 18446744073709551615
INV_q = 18446744073709551615
TWO_ADICITY = 67

ONE = [0xfffffffffffffff9,
       0xfffffff57fab5757,
       0x7f56ac056aeaf57f,
       0x10388572e3c2c0f5,
       0xe6ce591c2bafc343,
       0x3e03f4104144b1a]

one = 596565640619119817640570040948722104176104232228701042800312168817257919202191170334817591186065331324034272460793
inverse = 3298621403693021520254674166538938247613207941189372348401597750128415207145609847844581041878536076753325935356266

R_INV = pow(one, BN382_q - 2, BN382_q)
assert R_INV == 3298621403693021520254674166538938247613207941189372348401597750128415207145609847844581041878536076753325935356266
print(one * R_INV % BN382_q)

def limbs_to_bigint(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    s = R_INV * s % BN382_q
    return s

def print_name(name, num):
    s = str(num)
    print(name + " : " + " " * (140 - (len(s) + len(name))) + s)
    return

one = limbs_to_bigint(ONE)
print_name("R, or ONE", one),
# for Fq6
# non_residue :
"""
// u = sqrt(7)
// 3 * u has no cube nor square nor sixth root
/// NONRESIDUE = (2 * U) = 0, 3
"""
GENERATOR = [0xffffffffffffff9d,
        0xffffff6b7b52aeb7,
        0x76a537ba55d66b7f,
        0x2e8d16344c0b846b,
        0x2df8a320b2feee22,
        0x123eec4e5e4393ea]

print_name("generator", limbs_to_bigint(GENERATOR))

NONRESIDUE_0 = [0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
NONRESIDUE_1 = [0xffffffffffffffeb,
                0xffffffe07f020607,
                0x7e04041040c0e07f,
                0x30a99058ab4842e0,
                0xb46b0b54830f49c9,
                0xba0bdc30c3ce150]
print_name("non residue", limbs_to_bigint(NONRESIDUE_0))
print_name("non residue", limbs_to_bigint(NONRESIDUE_1))

# Frob_coeffs_c1[0]
FROB_COEFFS_C1 = [[0xfffffffffffffff9,
                0xfffffff57fab5757,
                0x7f56ac056aeaf57f,
                0x10388572e3c2c0f5,
                0xe6ce591c2bafc343,
                0x3e03f4104144b1a],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[1]
                [0x43ac10f69cd0866e,
                0xb67658d4844670fa,
                0x64500aac20e3e056,
                0xe69857d69abfc002,
                0x521ddf42ec5832c5,
                0xee09eba205fe5d8],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[2]
                [0x43ac10f69cd08675,
                0xb67658df049b19a2,
                0xe4f95ea6b5f8ead6,
                0xd65fd263b6fcff0c,
                0x6b4f8626c0a86f82,
                0xb005f791c4b9abd],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[3]
                [0x8,
                0xc0060c0c0,
                0xc1848c18180c00,
                0xa451b0a144d8480c,
                0x8a81e34d84edfc45,
                0x202449fed6c43c73],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[4]
                [0xbc53ef09632f7993,
                0x4989a72cfbc5a71d,
                0x1bc825e5621f2129,
                0xcdf1de3d8ddb48ff,
                0x1f325d26c4458cc2,
                0x1523ea85ba78a1b6],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c1[5]
                [0xbc53ef09632f798c,
                0x4989a7227b70fe75,
                0x9b1ed1eacd0a16a9,
                0xde2a63b0719e09f4,
                0x600b642eff55005,
                0x190429c6be8cecd1],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0]]

# Frob_coeffs_c2[0]
FROB_COEFFS_C2 = [[0xfffffffffffffff9,
                0xfffffff57fab5757,
                0x7f56ac056aeaf57f,
                0x10388572e3c2c0f5,
                0xe6ce591c2bafc343,
                0x3e03f4104144b1a],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c2[1]
                [0x43ac10f69cd08675,
                0xb67658df049b19a2,
                0xe4f95ea6b5f8ead6,
                0xd65fd263b6fcff0c,
                0x6b4f8626c0a86f82,
                0xb005f791c4b9abd],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c2[2]
                [0xbc53ef09632f7993,
                0x4989a72cfbc5a71d,
                0x1bc825e5621f2129,
                0xcdf1de3d8ddb48ff,
                0x1f325d26c4458cc2,
                0x1523ea85ba78a1b6],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c2[3]
                [0xfffffffffffffff9,
                0xfffffff57fab5757,
                0x7f56ac056aeaf57f,
                0x10388572e3c2c0f5,
                0xe6ce591c2bafc343,
                0x3e03f4104144b1a],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c2[4]
                [0x43ac10f69cd08675,
                0xb67658df049b19a2,
                0xe4f95ea6b5f8ead6,
                0xd65fd263b6fcff0c,
                0x6b4f8626c0a86f82,
                0xb005f791c4b9abd],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
# Frob_coeffs_c2[5]
                [0xbc53ef09632f7993,
                0x4989a72cfbc5a71d,
                0x1bc825e5621f2129,
                0xcdf1de3d8ddb48ff,
                0x1f325d26c4458cc2,
                0x1523ea85ba78a1b6],
                [0x0, 0x0, 0x0, 0x0, 0x0, 0x0]]
print("*****")
for x in FROB_COEFFS_C1:
    print_name("frobenius_coeffs_c1", limbs_to_bigint(x))
print("*****")
for x in FROB_COEFFS_C2:
    print_name("frobenius_coeffs_c2", limbs_to_bigint(x))

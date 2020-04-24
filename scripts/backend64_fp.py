# zexe/algebra/src/fields/bn_382/fq.rs
import math
BN382_q = 5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289
MODULUS_BITS = 384
LIMB_SIZES = [64, 64, 64, 64, 64, 64]

# R = one
one = 596565640619119817640570040948722104176104232228701042800312168817257919202191170334817591186065331324034272460793
inverse = 3298621403693021520254674166538938247613207941189372348401597750128415207145609847844581041878536076753325935356266
R_INV = pow(one, BN382_q - 2, BN382_q)
assert R_INV == inverse

def limbs_to_bigint_r(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    return s

def limbs_to_bigint(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    s = R_INV * s % BN382_q
    return s

MODULUS = [0x1,
    0x1800c1818,
    0x8018309183030180,
    0xb48a3614289b0901,
    0x71503c69b09dbf88,
    0x2404893fdad8878e]

R = [0xfffffffffffffff9,
    0xfffffff57fab5757,
    0x7f56ac056aeaf57f,
    0x10388572e3c2c0f5,
    0xe6ce591c2bafc343,
    0x3e03f4104144b1a]

R2 = [0xc79c121e98884701,
    0xfd75271b6a2e235d,
    0x1530439e68fe657,
    0xf6b7a72ebfbdbfb,
    0x50c6c2ce8f44951b,
    0x17fe189b54066561]


MODULUS_MINUS_ONE_DIV_TWO = [0x0,
    0xc0060c0c,
    0xc00c1848c18180c0,
    0x5a451b0a144d8480,
    0x38a81e34d84edfc4,
    0x1202449fed6c43c7]

T = [0x30018303,
    0x3003061230606030,
    0x169146c285136120,
    0xce2a078d3613b7f1,
    0x4809127fb5b10f1,
    0x0]

T_MINUS_ONE_DIV_TWO = [0x1800c181,
    0x1801830918303018,
    0x8b48a3614289b090,
    0xe71503c69b09dbf8,
    0x2404893fdad8878,
    0x0]

GENERATOR = 14
GENERATOR = [0xffffffffffffff9d,
    0xffffff6b7b52aeb7,
    0x76a537ba55d66b7f,
    0x2e8d16344c0b846b,
    0x2df8a320b2feee22,
    0x123eec4e5e4393ea]

ROOT_OF_UNITY = [0xe38be9090411d7d0,
        0x579d9745d8f8468b,
        0x4a5514233c9850c5,
        0xa7c5be912557804a,
        0xc69e67da380310d4,
        0x136e8eef9cf4445b]

MODULUS_BITS = 382
REPR_SHAVE_BITS = 2
INV = 18446744073709551615
TWO_ADICITY = 67

def print_name(name, num):
    s = str(num)
    print(name + " : " + " " * (140 - (len(s) + len(name))) + s)
    return

print("For scalar field")
modulus = limbs_to_bigint_r(MODULUS)
print_name("modulus_r", modulus)
RR = limbs_to_bigint_r(R)
print_name("R", RR)
# Rsquared (R = W^k, k = 4, 8, R a power of two)
R_SQUARE = (RR * RR) % modulus
print_name("RSquared", R_SQUARE)
assert R_SQUARE == limbs_to_bigint_r(R2)
# Rcubed
R_CUBED = pow(RR, 3, modulus)
print_name("RCubed", R_CUBED)
print_name("inv", hex(INV))
print_name("num_bits", MODULUS_BITS)
# eulers number (mod - 1 / 2)
euler = limbs_to_bigint(MODULUS_MINUS_ONE_DIV_TWO)
print_name("eulers number", euler)
print_name("s", TWO_ADICITY)
# t (as above, t odd)
t = limbs_to_bigint(T)
print_name("t", t)
# (t - 1) / 2
t_one_two = limbs_to_bigint(T_MINUS_ONE_DIV_TWO)
print_name("t - 1 / 2", t_one_two)
# s (mod = 2^s . t + 1)
MULT_GEN = 7
print_name("multiplicative_generator", MULT_GEN)
# root of unity (generator ^ ((mod - 1)/2s))
rou = limbs_to_bigint(ROOT_OF_UNITY)
print_name("root_of_unity", rou)
# nqr
nqr = 7
print_name("nqr", nqr)
# nqr^t
nqr_to_t = pow(nqr, t, modulus)
print_name("nqr^t", nqr_to_t)

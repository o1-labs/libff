# zexe/algebra/src/fields/bn_382/fp.rs
import math
BN382_p = 5543634365110765627805495722742127385843376434033820803590214255538854698464778703795540858859767700241957783601153
MODULUS_BITS = 384
LIMB_SIZES = [64, 64, 64, 64, 64, 64]

def limbs_to_bigint_r(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    return s

def print_name(name, num):
    s = str(num)
    print(name + " : " + " " * (140 - (len(s) + len(name))) + s)
    return

# from Fp.rs
MODULUS = [0x1,
           0x1800c1818,
           0x2012246d22424120,
           0xb48a3614289b0901,
           0x71503c69b09dbf88,
           0x2404893fdad8878e]

assert limbs_to_bigint_r(MODULUS) == BN382_p

R = [0xfffffffffffffff9,
    0xfffffff57fab5757,
    0x1f8101041030381f,
    0x10388572e3c2c0f8,
    0xe6ce591c2bafc343,
    0x3e03f4104144b1a]

R2 = [0xaa7b14a53b610887,
     0xb22034140d119ca9,
     0xe10d2796937ba75,
     0xe52454bf8b810402,
     0x1b4eec3d89fc0fd3,
     0xbc857aea27171f7]

# R = one
one = 596565640619119817640570040948722104176104232228701042800312168817257919202191170334817591186065331324034272460793
inverse = 186440667609470995775868266227789892795027616253919437514737345429576798232629712333450673685358560496014229019547

R_INV = pow(one, BN382_p - 2, BN382_p)
assert inverse == R_INV

def limbs_to_bigint(NUMBER_LIMBS):
    offset = [0, 64, 128, 192, 256, 320]
    s = 0
    for i in range(0, len(LIMB_SIZES)):
        s += 2**offset[i] * NUMBER_LIMBS[i]
    s = R_INV * s % BN382_p
    return s

MOD_MINUS_ONE_DIV_TWO = [0x0,
                        0xc0060c0c,
                        0x9009123691212090,
                        0x5a451b0a144d8480,
                        0x38a81e34d84edfc4,
                        0x1202449fed6c43c7]

T = [0x30018303,
    0x2402448da4484824,
    0x169146c285136120,
    0xce2a078d3613b7f1,
    0x4809127fb5b10f1,
    0x0]

T_MINUS_ONE_DIV_TWO = [0x1800c181,
                      0x12012246d2242412,
                      0x8b48a3614289b090,
                      0xe71503c69b09dbf8,
                      0x2404893fdad8878,
                      0x0]

GENERATOR = 7
GENERATOR_LIMBS = [0xffffffffffffffcf,
                   0xffffffb67daf6367,
                   0xdc87071c715188df,
                   0x718ba6243a5346c8,
                   0x4fa46fc531ce56d5,
                   0x1b21bac71c8e0dbc]

ROOT_OF_UNITY = [0xdb510d8c5d0d218f,
                0x447119a2f8d5e310,
                0x1373332ba33d5a84,
                0xb830356347b45dbb,
                0x851efb96cb691ec1,
                0x141037c57e9d0173]
MOD_BITS = 382;
CAPACITY = MODULUS_BITS - 1;
REPR_SHAVE_BITS = 2;
INV = 18446744073709551615;
TWO_ADICITY = 67;

def print_name(name, num):
    s = str(num)
    print(name + " : " + " " * (140 - (len(s) + len(name))) + s)
    return

print("For base field")
print_name("modulus", BN382_p)
RR = limbs_to_bigint_r(R)
print_name("R", RR)
# Rsquared (R = W^k, k = 4, 8, R a power of two)
R_SQUARE = (RR * RR) % BN382_p
print_name("RSquared", R_SQUARE)
# Rcubed
R_CUBED = pow(RR, 3, BN382_p)
print_name("RCubed", R_CUBED)
print_name("inv", hex(INV))
print_name("num_bits", MOD_BITS)
# eulers number (mod - 1 / 2)
euler = limbs_to_bigint(MOD_MINUS_ONE_DIV_TWO)
print_name("euler", euler)
print_name("s", TWO_ADICITY)
# t (as above, t odd)
t = limbs_to_bigint(T)
print_name("t", t)
# (t - 1) / 2
t_one_two = limbs_to_bigint(T_MINUS_ONE_DIV_TWO)
print_name("t - 1 / 2", t_one_two)
# s (mod = 2^s . t + 1)
print_name("multiplicative_generator", GENERATOR)
# root of unity (generator ^ ((mod - 1)/2s))
rou = limbs_to_bigint(ROOT_OF_UNITY)
print_name("root_of_unity", rou)
# nqr
# nqr^t
for nqr in range(2000):
    nqr_to_t = pow(nqr, t, BN382_p)
    if nqr_to_t == rou:
        print_name("nqr^t" + str(nqr), nqr_to_t)

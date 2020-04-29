import math
BN382_p = 5543634365110765627805495722742127385843376434033820803590214255538854698464778703795540858859767700241957783601153
BN382_q = 5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289
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


mod = BN382_q ** 2
mod_minus_one = mod - 1 
two_to_the_s = 0b100000000000000000000000000000000000000000000000000000000000000000000
print(math.log2(two_to_the_s))

t = (mod - 1)/ 2**68
print('{:.1f}'.format(t))

# fq2_s : int

# fq2_nqr_to_t : fq2 elt (fq, fq)
nqr = [0, 2]
t = 104123666252559599194082820731495841242513418965476562843279162460041891658890113434602162857661651873535350654374440669491702604974560757949770346152587382894324840566507020516312959205272607438081650000396288

nqr_to_t = pow(2, t, BN382_q)
print_name("nqr_to_t", nqr_to_t)

t_minus_1_over_2 = (t - 1 / 2)
print('{:.1f}'.format(t_minus_1_over_2))


# fq2_t_minus_1_over_2 : "mie::Vuint" (big int)



# for Fq2
NONRESIDUE = 7
NONRESIDUE_LIMBS = [0xffffffffffffffcf,
                    0xffffffb67daf6367,
                    0x7b5eb425ec6cb67f,
                    0x718ba6243a5346b6,
                    0x4fa46fc531ce56d5,
                    0x1b21bac71c8e0dbc]

# U = sqrt(7)
# QUADRATIC_NONRESIDUE = (0 + 2 * U)
QUADRATIC_NONRESIDUE = [0, 2]
QUADRATIC_NONRESIDUE_LIMBS = [[0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                               [0xfffffffffffffff2,
                               0xffffffeaff56aeaf,
                               0xfead580ad5d5eaff,
                               0x20710ae5c78581ea,
                               0xcd9cb238575f8686,
                               0x7c07e8208289635]]

# Coefficients for the Frobenius automorphism.
# [Fq(7)**(((q^0) - 1) / 2), Fq(7)**(((q^1) - 1) / 2)]
FROBENIUS_COEFF_FP2_C1 = [[0xfffffffffffffff9,
                           0xfffffff57fab5757,
                           0x7f56ac056aeaf57f,
                           0x10388572e3c2c0f5,
                           0xe6ce591c2bafc343,
                           0x3e03f4104144b1a],
                           [0x8,
                           0xc0060c0c0,
                           0xc1848c18180c00,
                           0xa451b0a144d8480c,
                           0x8a81e34d84edfc45,
                           0x202449fed6c43c73]]

print_name("nonresidue", limbs_to_bigint(NONRESIDUE_LIMBS))
for x in QUADRATIC_NONRESIDUE_LIMBS:
    print_name("quad nonresidue", limbs_to_bigint(x))
for x in FROBENIUS_COEFF_FP2_C1:
    print_name("frobenius_coeffs_c1", limbs_to_bigint(x))

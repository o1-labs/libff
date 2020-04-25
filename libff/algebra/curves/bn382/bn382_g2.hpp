/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_G2_HPP_
#define BN382_G2_HPP_
#include <iostream>
#include <vector>

#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/curves/bn382/bn382_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bn382_G2;
std::ostream& operator<<(std::ostream &, const bn382_G2&);
std::istream& operator>>(std::istream &, bn382_G2&);

class bn382_G2 {
private:
    static bn::Fp2 sqrt(const bn::Fp2 &el);
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bn382_G2 G2_zero;
    static bn382_G2 G2_one;

    bn::Fp2 coord[3];
    bn382_G2();
    typedef bn382_Fq base_field;
    typedef bn382_Fr scalar_field;

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bn382_G2 &other) const;
    bool operator!=(const bn382_G2 &other) const;

    bn382_G2 operator+(const bn382_G2 &other) const;
    bn382_G2 operator-() const;
    bn382_G2 operator-(const bn382_G2 &other) const;

    bn382_G2 add(const bn382_G2 &other) const;
    bn382_G2 mixed_add(const bn382_G2 &other) const;
    bn382_G2 dbl() const;

    bool is_well_formed() const;

    static bn382_G2 zero();
    static bn382_G2 one();
    static bn382_G2 random_element();

    static size_t size_in_bits() { return 2*base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bn382_G2 &g);
    friend std::istream& operator>>(std::istream &in, bn382_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<bn382_G2> &vec);
};

template<mp_size_t m>
bn382_G2 operator*(const bigint<m> &lhs, const bn382_G2 &rhs)
{
    return scalar_mul<bn382_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bn382_G2 operator*(const Fp_model<m, modulus_p> &lhs, const bn382_G2 &rhs)
{
    return scalar_mul<bn382_G2, m>(rhs, lhs.as_bigint());
}

} // libff
#endif // BN382_G2_HPP_

/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_G1_HPP_
#define BN382_G1_HPP_
#include <vector>

#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/curves/bn382/bn382_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bn382_G1;
std::ostream& operator<<(std::ostream &, const bn382_G1&);
std::istream& operator>>(std::istream &, bn382_G1&);

class bn382_G1 {
private:
    static bn::Fp sqrt(const bn::Fp &el);
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bn382_G1 G1_zero;
    static bn382_G1 G1_one;

    bn::Fp coord[3];
    bn382_G1();
    typedef bn382_Fq base_field;
    typedef bn382_Fr scalar_field;

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bn382_G1 &other) const;
    bool operator!=(const bn382_G1 &other) const;

    bn382_G1 operator+(const bn382_G1 &other) const;
    bn382_G1 operator-() const;
    bn382_G1 operator-(const bn382_G1 &other) const;

    bn382_G1 add(const bn382_G1 &other) const;
    bn382_G1 mixed_add(const bn382_G1 &other) const;
    bn382_G1 dbl() const;

    bool is_well_formed() const;

    static bn382_G1 zero();
    static bn382_G1 one();
    static bn382_G1 random_element();

    static size_t size_in_bits() { return bn382_Fq::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bn382_G1 &g);
    friend std::istream& operator>>(std::istream &in, bn382_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<bn382_G1> &vec);
};

template<mp_size_t m>
bn382_G1 operator*(const bigint<m> &lhs, const bn382_G1 &rhs)
{
    return scalar_mul<bn382_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bn382_G1 operator*(const Fp_model<m,modulus_p> &lhs, const bn382_G1 &rhs)
{
    return scalar_mul<bn382_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<bn382_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<bn382_G1> &v);


} // libff
#endif // BN382_G1_HPP_

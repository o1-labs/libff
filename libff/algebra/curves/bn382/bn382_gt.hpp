/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN382_GT_HPP_
#define BN382_GT_HPP_
#include <iostream>

#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/fields/field_utils.hpp>
#include <libff/algebra/fields/fp.hpp>

namespace libff {

class bn382_GT;
std::ostream& operator<<(std::ostream &, const bn382_GT&);
std::istream& operator>>(std::istream &, bn382_GT&);

class bn382_GT {
public:
    static bn382_GT GT_one;
    bn::Fp12 elem;

    bn382_GT();
    bool operator==(const bn382_GT &other) const;
    bool operator!=(const bn382_GT &other) const;

    bn382_GT operator*(const bn382_GT &other) const;
    bn382_GT unitary_inverse() const;

    static bn382_GT one();

    void print() { std::cout << this->elem << "\n"; };

    friend std::ostream& operator<<(std::ostream &out, const bn382_GT &g);
    friend std::istream& operator>>(std::istream &in, bn382_GT &g);
};

template<mp_size_t m>
bn382_GT operator^(const bn382_GT &rhs, const bigint<m> &lhs)
{
    return power<bn382_GT, m>(rhs, lhs);
}


template<mp_size_t m, const bigint<m>& modulus_p>
bn382_GT operator^(const bn382_GT &rhs, const Fp_model<m,modulus_p> &lhs)
{
    return power<bn382_GT, m>(rhs, lhs.as_bigint());
}

} // libff
#endif // BN382_GT_HPP_

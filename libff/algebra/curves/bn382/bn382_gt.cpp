/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn382/bn382_gt.hpp>

namespace libff {

bn382_GT bn382_GT::GT_one;
bn382_GT::bn382_GT()
{
    this->elem.clear();
}

bool bn382_GT::operator==(const bn382_GT &other) const
{
    return (this->elem == other.elem);
}

bool bn382_GT::operator!=(const bn382_GT& other) const
{
    return !(operator==(other));
}

bn382_GT bn382_GT::operator*(const bn382_GT &other) const
{
    bn382_GT result;
    bn::Fp12::mul(result.elem, this->elem, other.elem);
    return result;
}

bn382_GT bn382_GT::unitary_inverse() const
{
    bn382_GT result(*this);
    bn::Fp6::neg(result.elem.b_, result.elem.b_);
    return result;
}

bn382_GT bn382_GT::one()
{
    return GT_one;
}

std::ostream& operator<<(std::ostream &out, const bn382_GT &g)
{
#ifndef BINARY_OUTPUT
    out << g.elem.a_ << OUTPUT_SEPARATOR << g.elem.b_;
#else
    out.write((char*) &g.elem.a_, sizeof(g.elem.a_));
    out.write((char*) &g.elem.b_, sizeof(g.elem.b_));
#endif
    return out;
}

std::istream& operator>>(std::istream &in, bn382_GT &g)
{
#ifndef BINARY_OUTPUT
    in >> g.elem.a_;
    consume_OUTPUT_SEPARATOR(in);
    in >> g.elem.b_;
#else
    in.read((char*) &g.elem.a_, sizeof(g.elem.a_));
    in.read((char*) &g.elem.b_, sizeof(g.elem.b_));
#endif
    return in;
}
} // libff

#ifndef COMPLEX_H
#define COMPLEX_H

#include <fstream>
#include <cmath>

#define COMPLEX_DEFAULT_BASE double
#define i_unit_d i_unit<COMPLEX_DEFAULT_BASE>

template <typename BASE>
class Complex
{
private:
    BASE real_ = 0;
    BASE imag_ = 0;

public:
    BASE &real = real_;
    BASE &imag = imag_;

    Complex() {}
    ~Complex() {}

    Complex(BASE re) : real_(re) {}
    Complex(BASE re, BASE im) : real_(re), imag_(im) {}

    template <typename X>
    Complex(Complex<X> cm) : real_(cm.real), imag_(cm.imag) {}

    Complex conj() const {
        return Complex(real_, -imag_);
    }
    auto mag() const {
        return std::sqrt(real_*real_ + imag_*imag_);
    }
    auto norm() const {
        return real_*real_ + imag_*imag_;
    }
    auto arg() const {
        return atan2(real_, imag_);
    }


    template <typename X>
    bool operator==(const Complex<X> &rhs) const {
        return real_ == rhs.real && imag_ == rhs.imag;
    }

    template <typename X>
    bool operator!=(const Complex<X> &rhs) const {
        return !operator==(rhs);
    }

    template <typename X>
    auto operator+(const X &rhs) const {
        const Complex fixed(rhs);
        return Complex<decltype(real_+fixed.real)>(real_ + fixed.real, imag_ + fixed.imag);
    }

    template <typename X>
    auto operator-(const X &rhs) const {
        const Complex fixed(rhs);
        return Complex<decltype(real_+fixed.real)>(real_ - fixed.real, imag_ - fixed.imag);
    }

    template <typename X>
    auto operator*(const X &rhs) const {
        const Complex fixed(rhs);
        return Complex<decltype(real_+fixed.real)>(real_ * fixed.real - imag_ * fixed.imag, real_ * fixed.imag + imag_ * fixed.real);
    }

    template <typename X>
    auto operator/(const X &rhs) const {
        const Complex fixed(rhs);
        const auto norm = fixed.real * fixed.real + fixed.imag * fixed.imag;
        return Complex<decltype(real_+fixed.real)>((real_ * fixed.real + imag_ * fixed.imag) / norm, (imag * fixed.real - real_ * fixed.imag) / norm);
    }


    friend auto operator+(const BASE lhs, const Complex &cm) {
        return Complex<decltype(cm.real+lhs)>(lhs + cm.real_, cm.imag_);
    }

    friend auto operator-(const BASE lhs, const Complex &cm) {
        return Complex<decltype(cm.real+lhs)>(lhs - cm.real_, -cm.imag_);
    }

    friend auto operator*(const BASE lhs, const Complex &cm) {
        return Complex<decltype(cm.real_+lhs)>(lhs * cm.real_, lhs * cm.imag_);
    }

    friend auto operator/(const BASE lhs, const Complex &cm) {
        return lhs * cm.conj() / cm.norm();
    }


    template <typename X>
    Complex& operator+=(X &rhs) {
        Complex<BASE> newc(*this + (rhs));
        real_ = newc.real;
        imag_ = newc.imag;
        return *this;
    }
    template <typename X>
    Complex& operator-=(X &rhs) {
        Complex<BASE> newc(*this - (rhs));
        real_ = newc.real;
        imag_ = newc.imag;
        return *this;
    }
    template <typename X>
    Complex& operator*=(X &rhs) {
        Complex<BASE> newc(*this * (rhs));
        real_ = newc.real;
        imag_ = newc.imag;
        return *this;
    }
    template <typename X>
    Complex& operator/=(X &rhs) {
        Complex<BASE> newc(*this / (rhs));
        real_ = newc.real;
        imag_ = newc.imag;
        return *this;
    }


    friend std::ostream &operator<<(std::ostream& os, const Complex<BASE> &cm) {
        return (os << '(' << cm.real << ", " << cm.imag << ')');
    }

};

template <typename BASE = COMPLEX_DEFAULT_BASE>
const Complex<BASE> i_unit = Complex<BASE>(0, 1);

typedef Complex<COMPLEX_DEFAULT_BASE> Complex_t;


#endif

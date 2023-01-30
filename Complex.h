#ifndef COMPLEX_H
#define COMPLEX_H

#include <fstream>
#include <cmath>
#include <type_traits>

#define COMPLEX_DEFAULT_BASE double

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
        return std::atan2(imag_, real_);
    }

    Complex operator-() const {
        return this->operator*(-1);
    }




    
    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator+(const X &rhs) {
        return this->operator+(Complex<X>(rhs, 0));
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator-(const X &rhs) {
        return this->operator-(Complex<X>(rhs, 0));
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator*(const X &rhs) {
        return this->operator*(Complex<X>(rhs, 0));
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator/(const X &rhs) {
        return this->operator/(Complex<X>(rhs, 0));
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
    Complex<decltype(BASE(1) + X(1))> operator+(const Complex<X> rhs) const {
        return Complex<decltype(real_+rhs.real)>(real_ + rhs.real, imag_ + rhs.imag);
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator-(const Complex<X> rhs) const {
        return Complex<decltype(real_+rhs.real)>(real_ - rhs.real, imag_ - rhs.imag);
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator*(const Complex<X> rhs) const {
        return Complex<decltype(real_+rhs.real)>(real_ * rhs.real - imag_ * rhs.imag, real_ * rhs.imag + imag_ * rhs.real);
    }

    template <typename X>
    Complex<decltype(BASE(1) + X(1))> operator/(const Complex<X> rhs) const {
        const auto norm = rhs.real * rhs.real + rhs.imag * rhs.imag;
        return Complex<decltype(real_+rhs.real)>((real_ * rhs.real + imag_ * rhs.imag) / norm, (imag * rhs.real - real_ * rhs.imag) / norm);
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

template <typename>
struct isComplex : public std::false_type
{};

template <typename T>
struct isComplex<Complex<T>> : public std::true_type
{};

template <typename T>
constexpr bool isComplexFunc(T const &) {
    return isComplex<T>::value;
}

const Complex<int> i_unit = Complex<int>(0, 1);


template <typename X>
auto abs(const Complex<X> &cm) {
    return cm.mag();
}

template <typename X>
auto arg(const Complex<X> &cm) {
    return cm.arg();
}

template <typename X>
auto exp(const Complex<X> &cm) {
    return exp(cm.real) * Complex(std::cos(cm.imag), std::sin(cm.imag));
}

template <typename X>
auto log(const Complex<X> &cm) {
    return Complex(std::log(cm.mag()), cm.arg());
}

template <typename X, typename Y>
auto pow(const Complex<X> &lhs, const Complex<Y> &rhs) {
    return exp(log(lhs) * rhs);
}


template <typename X>
auto sin(const Complex<X> &cm) {
    return Complex(std::cosh(cm.imag)*std::sin(cm.real), std::cos(cm.real)*std::sinh(cm.imag));
}

template <typename X>
auto cos(const Complex<X> &cm) {
    return Complex(std::cos(cm.real)*std::cosh(cm.imag), -std::sin(cm.real)*std::sinh(cm.imag));
}

template <typename X>
auto tan(const Complex<X> &cm) {
    return Complex(std::sin(2*cm.real), std::sinh(2*cm.imag)) / (std::cos(2*cm.real) + std::cosh(2*cm.imag));
}

template <typename X>
auto sinh(const Complex<X> &cm) {
    return sin(cm * -i_unit) * i_unit;
}

template <typename X>
auto cosh(const Complex<X> &cm) {
    return cos(cm * i_unit);
}

template <typename X>
auto tanh(const Complex<X> &cm) {
    return sinh(cm)/cosh(cm);
}


#endif

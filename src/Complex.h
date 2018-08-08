#pragma once

struct Complex
{
    double r;
    double i;

    Complex() 
    {
        this->r = 0;
        this->i = 0;
    }

    Complex(double r, double i) 
    {
        this->r = r;
        this->i = i;
    }

    Complex& operator+=(const Complex& rhs)
    {
        this->r += rhs.r;
        this->i += rhs.i;
        return *this;
    }

    Complex& operator-=(const Complex& rhs)
    {
        this->r -= rhs.r;
        this->i -= rhs.i;
        return *this;
    }



};

inline Complex operator+(Complex lhs, const Complex& rhs)
{
    lhs += rhs;
    return lhs;
}

inline Complex operator-(Complex lhs, const Complex& rhs)
{
    lhs -= rhs;
    return lhs;
}

inline Complex operator*(const Complex& lhs, const Complex& rhs)
{
    return Complex(lhs.r * rhs.r - lhs.i * rhs.i, lhs.r * rhs.i + lhs.i * rhs.r);
}

inline Complex operator*(const Complex& lhs, double rhs)
{
    return Complex(lhs.r * rhs, lhs.i * rhs);
}
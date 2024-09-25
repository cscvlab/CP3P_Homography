#pragma once

#ifndef QUARTICASOLVER
#define QUARTICASOLVER
#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

namespace hd{
    template<class T>
    int solveQuartic(const T *factors, T *realRoots) {
        const T &a4 = factors[0];
        const T &a3 = factors[1];
        const T &a2 = factors[2];
        const T &a1 = factors[3];
        const T &a0 = factors[4];

        T a4_2 = a4 * a4;
        T a3_2 = a3 * a3;
        T a4_3 = a4_2 * a4;
        T a2a4 = a2 * a4;

        T p4 = (8 * a2a4 - 3 * a3_2) / (8 * a4_2);
        T q4 = (a3_2 * a3 - 4 * a2a4 * a3 + 8 * a1 * a4_2) / (8 * a4_3);
        T r4 =
                (256 * a0 * a4_3 - 3 * (a3_2 * a3_2) - 64 * a1 * a3 * a4_2 + 16 * a2a4 * a3_2) / (256 * (a4_3 * a4));

        T p3 = ((p4 * p4) / 12 + r4) / 3;                             // /=-3
        T q3 = (72 * r4 * p4 - 2 * p4 * p4 * p4 - 27 * q4 * q4) / 432;// /=2

        T t;// *=2
        complex<T> w;
        if (q3 >= 0)
            w = -sqrt(static_cast<complex<T>>(q3 * q3 - p3 * p3 * p3)) - q3;
        else
            w = sqrt(static_cast<complex<T>>(q3 * q3 - p3 * p3 * p3)) - q3;
        if (w.imag() == 0.0) {
            w.real(cbrt(w.real()));
            t = 2.0 * (w.real() + p3 / w.real());
        } else {
            w = pow(w, 1.0 / 3);
            t = 4.0 * w.real();
        }

        complex<T> sqrt_2m = sqrt(static_cast<complex<T>>(-2 * p4 / 3 + t));
        T B_4A = -a3 / (4 * a4);
        T complex1 = 4 * p4 / 3 + t;
        complex<T> complex2 = 2 * q4 / sqrt_2m;

        T sqrt_2m_rh = sqrt_2m.real() / 2;
        T sqrt_2m_rh_i = sqrt_2m.imag() / 2;

        complex<T> temp1 = sqrt(-(complex1 + complex2));
        complex<T> temp2 = sqrt(-(complex1 - complex2));


        T sqrt1 = temp1.real() / 2;
        T sqrt2 = temp2.real() / 2;
        T sqrt1_i = temp1.imag() / 2;
        T sqrt2_i = temp2.imag() / 2;

        int sol_num = 0;

        if (B_4A + sqrt_2m_rh + sqrt1 > 0 && abs(sqrt_2m_rh_i + sqrt1_i) < 1e-8) {
            realRoots[sol_num++] = B_4A + sqrt_2m_rh + sqrt1;
        }
        if (B_4A + sqrt_2m_rh - sqrt1 > 0 && abs(sqrt_2m_rh_i - sqrt1_i) < 1e-8) {
            realRoots[sol_num++] = B_4A + sqrt_2m_rh - sqrt1;
        }
        if (B_4A - sqrt_2m_rh + sqrt2 > 0 && abs(sqrt_2m_rh_i + sqrt2_i) < 1e-8) {
            realRoots[sol_num++] = B_4A - sqrt_2m_rh + sqrt2;
        }
        if (B_4A - sqrt_2m_rh - sqrt2 > 0 && abs(sqrt_2m_rh_i - sqrt2_i) < 1e-8) {
            realRoots[sol_num++] = B_4A - sqrt_2m_rh - sqrt2;
        }

        return sol_num;
    }

    template<class T>
    T cubic_TC(T k3, T k2, T k1, T k0) {
        T k3_inv = 1.0 / k3;
        T k2_norm = k3_inv * k2;
        T k1_norm = k3_inv * k1;
        T k0_norm = k3_inv * k0;
        T consOneThird = 1.0 / 3.0;

        T k2_norm_13 = k2_norm * consOneThird;
        T alpha_oneThird = k1_norm * consOneThird - k2_norm_13 * k2_norm_13;
        T beta_half = (k0_norm - k1_norm * k2_norm_13) * 0.5 + k2_norm_13 * k2_norm_13 * k2_norm_13;
        T gamma = beta_half * beta_half + alpha_oneThird * alpha_oneThird * alpha_oneThird;

        T s;
        if (gamma < 0) {
            // Handle the case where gamma is less than 0
            T squareRoot_1 = sqrt(-alpha_oneThird);
            T temp1 = beta_half / (alpha_oneThird * squareRoot_1);
            T cosAcosValue = cos(acos(temp1) * consOneThird);
            s = 2.0 * squareRoot_1 * cosAcosValue - k2_norm_13;
        } else if (gamma > 0) {
            // Handle the case where gamma is greater than 0
            T gamma_root = sqrt(gamma);
            T cubicRoot_1 = cbrt(-beta_half + gamma_root);
            T cubicRoot_2 = cbrt(beta_half + gamma_root);
            s = cubicRoot_1 - cubicRoot_2 - k2_norm_13;
        } else {

            s = -k2_norm_13 + (alpha_oneThird != 0 ? (2.0 * beta_half / alpha_oneThird) : 0);
        }

        // do a single newton step on the cubic equation. Flops: 19
        // same to “P3P.cc” in PoseLib, implemented for ECCV18.
        T s_v = s * s * s + k2_norm * s * s + k1_norm * s + k0_norm;
        T ds = 3.0 * s * s + 2.0 * k2_norm * s + k1_norm;
        s = s - s_v / ds;

        return s;
    }
}
#endif 

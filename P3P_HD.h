#pragma once
#include <cmath>
#include <complex>
#include <iostream>
// #include "utils/matrix.h"
#include <Eigen/Core>
#include "utils/quarticSolver.h"

using namespace std;
namespace hd {

    enum class P3P_methods {
        hd_mul3=0,
        hd_mul4=1,
        hd_uni3_1v=2,
        hd_uni3_3v=3
    };

// imagePoints: The 3 bearing measurements (normalized) stored as column vectors.
// worldPoints: The positions of the 3 feature points stored as column vectors.
// solutions: Output of this function. Column i (i=0,4,8,12) is the solution for the camera position, and column i+1 to
//            i+3 h.
//Rs: Multiple R Matrices
//Ts: Multiple t Matrices
//iterMax: number denoting the maximal iteration of polishing.
//method_id: a flag denoting to decide to use 'Mul3','Mul4','Uni_3_1V' or 'Uni_3_3V'.

    template<class T>
    int P3P_HD(T *imagePoints, T *worldPoints,T Intrin[3],std::vector<Eigen::Matrix<T, 3, 3>> &Rs,
        std::vector<Eigen::Vector3<T>> &Ts,int iterMax, P3P_methods method_id) {

        // STAGE 1: Constructing two quadratic eqautions.FLOPs : 21 + 31 + 14 = 66
        // 1 - 1. Constants of three object points M, N, and P.Flops : 3 * 2 + 5 * 3 = 21
        // two 3D vector in world frame 'Lambda'
        T MN[3] = {worldPoints[1] - worldPoints[0], worldPoints[5] - worldPoints[4], worldPoints[9] - worldPoints[8]};
        T MP[3] = {worldPoints[2] - worldPoints[0], worldPoints[6] - worldPoints[4], worldPoints[10] - worldPoints[8]};

        T d_MN = MN[0] * MN[0] + MN[1] * MN[1] + MN[2] * MN[2];
        T d_MP = MP[0] * MP[0] + MP[1] * MP[1] + MP[2] * MP[2];
        T d_MNP = MN[0] * MP[0] + MN[1] * MP[1] + MN[2] * MP[2];

        // 1 - 2. Variables related to three corresponding image points M2, N2, and P2.Flops: 3 * 2 + 1 + 6 * 4 = 31
        // 3D vector from camera center C to 3D image point; also denotes the coordinates in the camera frame 'Gamma'
        T e_M2[3] = {imagePoints[0] - Intrin[1], imagePoints[4] - Intrin[2], Intrin[0]};
        T e_N2[3] = {imagePoints[1] - Intrin[1], imagePoints[5] - Intrin[2], Intrin[0]};
        T e_P2[3] = {imagePoints[2] - Intrin[1], imagePoints[6] - Intrin[2], Intrin[0]};
        //Six related constants, also computed using dot product.
        T f2_A = Intrin[0] * Intrin[0];
        T l_MM = e_M2[0] * e_M2[0] + e_M2[1] * e_M2[1] + f2_A;
        T l_NN = e_N2[0] * e_N2[0] + e_N2[1] * e_N2[1] + f2_A;
        T l_PP = e_P2[0] * e_P2[0] + e_P2[1] * e_P2[1] + f2_A;
        T l_MN = e_M2[0] * e_N2[0] + e_M2[1] * e_N2[1] + f2_A;
        T l_MP = e_M2[0] * e_P2[0] + e_M2[1] * e_P2[1] + f2_A;
        T l_NP = e_P2[0] * e_N2[0] + e_P2[1] * e_N2[1] + f2_A;

        //1-3. Compute the coefficients of two quadratic equations f1-f6 and g1-g6. Flops: 14
        T f1 = d_MNP * l_NN;
        T f2 = -d_MN * l_NP;
        T temp = d_MNP - d_MN;
        T f4 = -l_MN * (d_MNP + temp);
        T f5 = d_MN * l_MP;
        T f6 = l_MM * temp;
        T g1 = d_MP * l_NN;
        T g3 = -d_MN * l_PP;
        T g4_half = -d_MP * l_MN;
        T g4 = 2 * g4_half;
        T g5 = 2 * f5;
        T g6 = l_MM * (d_MP - d_MN);

        int solu_size = 0;
        T lam1_solu[4] = {0};
        T lam2_solu[4] = {0};

        // STAGE 2: Solving a quartic or cubic equaiton.
        if (method_id==P3P_methods::hd_mul4) {
            T f22 = f2 * f2;
            T f25 = f2 * g5;
            T f55 = f5 * f5;
            T f4_2 = 2 * f4;
            T g3f1 = g3 * f1;
            T g3f6 = g3 * f6;
            T gf1 = g1 - f1;
            T gf4 = g4 - f4;
            T gf6 = g6 - f6;

            T c1 = g1 * f22 + g3f1 * f1;
            T c2 = f22 * g4 + g3f1 * f4_2 + f25 * gf1;
            T c3 = f55 * (gf1 - f1) + 2 * g3f1 * f6 + f22 * g6 + g3 * f4 * f4 + f25 * gf4;
            T c4 = f55 * (gf4 - f4) + g3f6 * f4_2 + f25 * gf6;
            T c5 = f55 * (gf6 - f6) + g3f6 * f6;

            T c[5] = {c1, c2, c3, c4, c5};

            T solu[4] = {0};
            int quartic_solu_size = hd::solveQuartic(c, solu);


            for (int i = 0; i < quartic_solu_size; i++) {
                T b_temp = solu[i];
                T v = ((f1 * b_temp + f4) * b_temp + f6) / (-f2 * b_temp - f5);
                if (v > 0) {
                    lam1_solu[solu_size] = b_temp;
                    lam2_solu[solu_size] = v;
                    solu_size++;
                }
            }
        }

        else {

            T f1f55 = f1 * f5 * f5;
            T g1f55 = g1 * f5 * f5;

            T k3 = (g3 * g4_half * g4_half) + g1f55 - g1 * g3 * g6;
            T k2 = g3 * f1 * g6 - f1f55 + f6 * g1 * g3 + (f2 * f5 - f4 * g3) * g4_half - g1f55;
            T k1 = (g6 * f2 * f2 + g1f55 + g3 * f4 * f4 - g5 * f2 * f4 - g4 * f2 * f5) * 0.25 + f1f55 - g3 * f1 * f6;
            T k0 = (f2 * f4 * f5 - f1f55 - f2 * f2 * f6) * 0.25;

            T s = hd::cubic_TC(k3, k2, k1, k0);

            // Calculate CC_33
            T CC_11 = f1 - g1 * s;
            T CC_12 = f2 * 0.5;
            T CC_13 = (f4 - g4 * s) * 0.5;
            T CC_22 = -g3 * s;
            T CC_23 = (0.5 - s) * f5;
            T CC_33 = f6 - g6 * s;
            T CInv_11 = CC_23 * CC_23 - CC_22 * CC_33;
            T CInv_22 = CC_13 * CC_13 - CC_11 * CC_33;
            T CInv_33 = CC_12 * CC_12 - CC_11 * CC_22;

            T v2 = 0, v3 = 0;
            if (CInv_11 > CInv_22) {
                if (CInv_11 > CInv_33) {
                    T v1_inv = 1.0 / sqrt(CInv_11);
                    v2 = (CC_12 * CC_33 - CC_13 * CC_23) * v1_inv;
                    v3 = (CC_13 * CC_22 - CC_12 * CC_23) * v1_inv;
                } else {
                    v3 = sqrt(CInv_33);
                    v2 = (CC_11 * CC_23 - CC_12 * CC_13) / v3;
                }
            } else if (CInv_22 > CInv_33) {
                v2 = sqrt(CInv_22);
                v3 = (CC_11 * CC_23 - CC_12 * CC_13) / v2;
            } else {
                v3 = sqrt(CInv_33);
                v2 = (CC_11 * CC_23 - CC_12 * CC_13) / v3;
            }

            T p1, p2, p3, aaa, bbb, ccc, ddd;
            for (int ii = 0; ii < 2; ++ii) {
                p1 = CC_11;
                p2 = (ii == 0) ? CC_12 - v3 : CC_12 + v3;
                p3 = (ii == 0) ? CC_13 + v2 : CC_13 - v2;
                aaa = f1 * p2 - f2 * p1;
                bbb = f4 * p2 - f5 * p1 - f2 * p3;
                ccc = f6 * p2 - f5 * p3;
                ddd = bbb * bbb - 4 * aaa * ccc;

                if (ddd < 0) {
                    continue;
                }

                T ddd_sqrt = sqrt(ddd);
                T temp1 = -bbb + ddd_sqrt;
                T temp2 = -bbb - ddd_sqrt;
                T temp3 = 0.5 / aaa;
                T temp4 = 1.0 / (-p2);
                T lam1_temp1 = temp1 * temp3;
                T lam2_temp1 = (p3 + p1 * lam1_temp1) * temp4;
                T lam1_temp2 = temp2 * temp3;
                T lam2_temp2 = (p3 + p1 * lam1_temp2) * temp4;

                if (lam2_temp1 > 0 && lam1_temp1 > 0) {
                    lam1_solu[solu_size] = lam1_temp1;
                    lam2_solu[solu_size] = lam2_temp1;
                    solu_size++;
                }


                if (lam2_temp2 > 0 && lam1_temp2 > 0) {
                    lam1_solu[solu_size] = lam1_temp2;
                    lam2_solu[solu_size] = lam2_temp2;
                    solu_size++;
                }
            }

        }
        if (solu_size == 0) {
            return 0;
        }

        // STAGE 3. Use the fourth point Q to remove ambiguous solutions.
        // hd_uni-1v:Flops: 41 + 29 * m, m is the number of real solutions of{ lam_1, lam_2 }
        if (method_id==P3P_methods::hd_uni3_1v) {
            // 3-1. Constants related to the 4th object point Q.Flops: 3+5*2 = 13
            T MQ[3] = {worldPoints[3] - worldPoints[0],
                       worldPoints[7] - worldPoints[4],
                       worldPoints[11] - worldPoints[8]};

            T d_MQ = MQ[0] * MQ[0] + MQ[1] * MQ[1] + MQ[2] * MQ[2];
            T d_MNQ = MN[0] * MQ[0] + MN[1] * MQ[1] + MN[2] * MQ[2];

            // 3-2. Constants related to the 4th image point Q2.Flops: 2+4*3 = 14
            T e_Q2[3] = {imagePoints[3] - Intrin[1], imagePoints[7] - Intrin[2], Intrin[0]};
            T l_QQ = e_Q2[0] * e_Q2[0] + e_Q2[1] * e_Q2[1] + f2_A;
            T l_MQ = e_M2[0] * e_Q2[0] + e_M2[1] * e_Q2[1] + f2_A;
            T l_NQ = e_Q2[0] * e_N2[0] + e_Q2[1] * e_N2[1] + f2_A;

            // 3-3. Compute the coefficients of two quadratic equations formed by three points M, N, and Q.Flops: 14
            T f1_1 = d_MNQ * l_NN;
            T f2_1 = -d_MN * l_NQ;
            T temp_1 = d_MNQ - d_MN;
            T f4_1 = -l_MN * (temp_1 + d_MNQ);
            T f5_1 = d_MN * l_MQ;
            T f6_1 = l_MM * temp_1;

            // 3-4. Remove ambiguous solutions. Flops: 31*m
            T min_value = INT32_MAX;
            for (int ii = 0; ii < solu_size; ++ii) {
                T b_1 = lam1_solu[ii];
                T v_1 = -((f1_1 * b_1 + f4_1) * b_1 + f6_1) / (f2_1 * b_1 + f5_1);
                T value = std::abs(d_MQ*(l_NN*b_1*b_1 - 2*l_MN*b_1+l_MM) / (d_MN*(l_QQ*v_1*v_1 - 2*l_MQ*v_1+l_MM)) - 1);

                if (value < min_value) {
                    lam1_solu[0] = b_1;
                    lam2_solu[0] = lam2_solu[ii];
                    min_value = value;

                }
            }
            solu_size = 1;
        }
        // hd_uni-3v:Flops: 56 + 45 * m, m is the number of real solutions of{ lam_1, lam_2 }
        if (method_id==P3P_methods::hd_uni3_3v) {
            // 3-1. Constants related to the 4th object point Q.Flops: 3*3 + 5*4= 29
            T MQ[3] = {worldPoints[3] - worldPoints[0],
                       worldPoints[7] - worldPoints[4],
                       worldPoints[11] - worldPoints[8]};
            T NQ[3] = {worldPoints[3] - worldPoints[1],
                       worldPoints[7] - worldPoints[5],
                       worldPoints[11] - worldPoints[9]};
            T PQ[3] = {worldPoints[3] - worldPoints[2],
                       worldPoints[7] - worldPoints[6],
                       worldPoints[11] - worldPoints[10]};

            T d_NQ = NQ[0] * NQ[0] + NQ[1] * NQ[1] + NQ[2] * NQ[2];
            T d_MQ = MQ[0] * MQ[0] + MQ[1] * MQ[1] + MQ[2] * MQ[2];
            T d_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];
            T d_MNQ = MN[0] * MQ[0] + MN[1] * MQ[1] + MN[2] * MQ[2];

            // 3-2. Constants related to the 4th image point Q2.Flops: 2+4*4 = 18
            T e_Q2[3] = {imagePoints[3] - Intrin[1], imagePoints[7] - Intrin[2], Intrin[0]};
            T l_QQ = e_Q2[0] * e_Q2[0] + e_Q2[1] * e_Q2[1] + f2_A;
            T l_MQ = e_M2[0] * e_Q2[0] + e_M2[1] * e_Q2[1] + f2_A;
            T l_NQ = e_Q2[0] * e_N2[0] + e_Q2[1] * e_N2[1] + f2_A;
            T l_PQ = e_Q2[0] * e_P2[0] + e_Q2[1] * e_P2[1] + f2_A;

            //3-3. Compute the coefficients of two quadratic equations formed by three points M, N, and Q.Flops: 9
            T f1_1 = d_MNQ * l_NN;
            T f2_1 = -d_MN * l_NQ;
            T temp_1 = d_MNQ - d_MN;
            T f4_1 = -l_MN * (temp_1 + d_MNQ);
            T f5_1 = d_MN * l_MQ;
            T f6_1 = l_MM * temp_1;

            // 3-4. Remove ambiguous solutions.Flops: (11+13+19+2)*m=45*m
            T min_value = INT32_MAX;
            for (int ii = 0; ii < solu_size; ++ii) {
                T b_1 = lam1_solu[ii];
                T b_2 = lam2_solu[ii];
                T v_1 = -((f1_1 * b_1 + f4_1) * b_1 + f6_1) / (f2_1 * b_1 + f5_1);

                T temp_NN=l_NN*b_1*b_1;
                T temp_QQ=l_QQ*v_1*v_1;
                T temp_MN=d_MN/(temp_NN - 2*l_MN*b_1+l_MM);
                T temp_2V=2*v_1;

                // value_dmn_dmq
                T value_dmn_dmq = std::abs((temp_QQ - temp_2V*l_MQ+l_MM)*temp_MN - d_MQ);
                // value_dmn_dnq
                T value_dmn_dnq = std::abs((temp_QQ - temp_2V*l_NQ*b_1+temp_NN)*temp_MN - d_NQ);
                // value_dmn_dpq
                T value_dmn_dpq = std::abs((temp_QQ - temp_2V*l_PQ*b_2+l_PP*b_2*b_2)*temp_MN - d_PQ);

                T theta_3PSE=value_dmn_dmq + value_dmn_dnq + value_dmn_dpq;

                if (theta_3PSE < min_value) {
                    lam1_solu[0] = b_1;
                    lam2_solu[0] = lam2_solu[ii];
                    min_value = theta_3PSE;
                }

            }
            solu_size = 1;

        }
        // STAGE 4 Polishing with Gauss - Newton method.Our two objective functions. Flops:53*n + 23.
        for (int ii = 0; ii < solu_size; ++ii) {
            T lam_1 = lam1_solu[ii];
            T lam_2 = lam2_solu[ii];
            T err_min = INT32_MAX;
            T lam_1_min = lam_1;
            T lam_2_min = lam_2;
            for (int iter = 1; iter <= iterMax; ++iter) {
                T r1 = f1 * lam_1 * lam_1 + f2 * lam_1 * lam_2 + f4 * lam_1 + f5 * lam_2 + f6;
                T r2 = g1 * lam_1 * lam_1 + g3 * lam_2 * lam_2 + g4 * lam_1 + g5 * lam_2 + g6;
                T err_cur = std::abs(r1) + std::abs(r2);

                if (err_cur < err_min) {
                    err_min = err_cur;
                    lam_1_min = lam_1;
                    lam_2_min = lam_2;
                }
                if (err_cur < 1e-7) {
                    break;
                }
                T x1_1 = 2 * f1 * lam_1 + f2 * lam_2 + f4;
                T x1_2 = f2 * lam_1 + f5;
                T x2_1 = 2 * g1 * lam_1 + g4;
                T x2_2 = 2 * g3 * lam_2 + g5;

                T detJ = 1. / (x1_1 * x2_2 - x1_2 * x2_1);// Half minus inverse determinant

                lam_1 = lam_1 - (x2_2 * r1 - x1_2 * r2) * detJ;
                lam_2 = lam_2 - (-x2_1 * r1 + x1_1 * r2) * detJ;

                if (iter == iterMax) {
                    r1 = f1 * lam_1 * lam_1 + f2 * lam_1 * lam_2 + f4 * lam_1 + f5 * lam_2 + f6;
                    r2 = g1 * lam_1 * lam_1 + g3 * lam_2 * lam_2 + g4 * lam_1 + g5 * lam_2 + g6;
                    err_cur = std::abs(r1) + std::abs(r2);

                    if (err_cur < err_min) {
                        err_min = err_cur;
                        lam_1_min = lam_1;
                        lam_2_min = lam_2;
                    }
                }
            }
            lam1_solu[ii] = lam_1_min;
            lam2_solu[ii] = lam_2_min;
        }

        // STAGE 5 Recover the complete pose[R, t]. 34 flops + (118 flops + 1 * sqrt) * m
        // m=1:
        // m>1:
        // 5-1. Compute the normal vector MS in the object frame. Flops: 16
        T MS[3] = {MN[1] * MP[2] - MN[2] * MP[1],
                   MN[2] * MP[0] - MN[0] * MP[2],
                   MN[0] * MP[1] - MN[1] * MP[0]};
        T d_MS_inv = 1 / (d_MN * d_MP - d_MNP * d_MNP);

        Eigen::Matrix<T, 3, 3> R_;
        Eigen::Vector3<T> t_;

        if (solu_size == 1) {

            T lam_1 = lam1_solu[0];
            T lam_2 = lam2_solu[0];
            //5-2. Compute the normal vector M2S4 of the image-related triangle. Flops: 21
            Eigen::Vector3<T> M2N4 = {
                    lam_1 * e_N2[0] - e_M2[0],
                    lam_1 * e_N2[1] - e_M2[1],
                    lam_1 * e_N2[2] - e_M2[2]};

            Eigen::Vector3<T> M2P4 = {
                    lam_2 * e_P2[0] - e_M2[0],
                    lam_2 * e_P2[1] - e_M2[1],
                    lam_2 * e_P2[2] - e_M2[2]};

            Eigen::Vector3<T> M2S4 = {
                    M2N4(1) * M2P4(2) - M2N4(2) * M2P4(1),
                    M2N4(2) * M2P4(0) - M2N4(0) * M2P4(2),
                    M2N4(0) * M2P4(1) - M2N4(1) * M2P4(0)};

            //5-3. Compute lam_M. Flops: 9 + 1*sqrt
            T lam_M_2 = d_MN / (M2N4(0) * M2N4(0) + M2N4(1) * M2N4(1) + M2N4(2) * M2N4(2));
            T lam_M = sqrt(lam_M_2);

            //5-4. Compute an immediate matrix [u1'; u2'; u3'], s.t. R = [M2N4 M2P4 M2S4] * [u1'; u2'; u3'];. Flops: 26
            T ff0 = lam_M_2 * d_MS_inv;
            T ff1 = lam_M * d_MS_inv;
            T ff2 = ff1 * d_MNP;
            T ff3 = ff1 * d_MP;
            T ff4 = ff1 * d_MN;

            T u1[3] = {ff3 * MN[0] - ff2 * MP[0],
                       ff3 * MN[1] - ff2 * MP[1],
                       ff3 * MN[2] - ff2 * MP[2]};
            T u2[3] = {ff4 * MP[0] - ff2 * MN[0],
                       ff4 * MP[1] - ff2 * MN[1],
                       ff4 * MP[2] - ff2 * MN[2]};
            T u3[3] = {ff0 * MS[0],
                       ff0 * MS[1],
                       ff0 * MS[2]};

            //5-5. Compute R, t. Flops: 66
            Eigen::Matrix<T, 3, 3> M1, M2;
            M1 << M2N4[0], M2P4[0], M2S4[0],
                  M2N4[1], M2P4[1], M2S4[1],
                  M2N4[2], M2P4[2], M2S4[2];

            M2 << u1[0], u1[1], u1[2],
                  u2[0], u2[1], u2[2],
                  u3[0], u3[1], u3[2];

            R_ = M1 * M2;

            t_[0] = lam_M * e_M2[0] -
                    (R_(0, 0) * worldPoints[0] + R_(0, 1) * worldPoints[4] + R_(0, 2) * worldPoints[8]);
            t_[1] = lam_M * e_M2[1] -
                    (R_(1, 0) * worldPoints[0] + R_(1, 1) * worldPoints[4] + R_(1, 2) * worldPoints[8]);
            t_[2] = lam_M * e_M2[2] -
                    (R_(2, 0) * worldPoints[0] + R_(2, 1) * worldPoints[4] + R_(2, 2) * worldPoints[8]);

            Ts[0] = t_;
            Rs[0] = R_;

        }
        else {
            T d1 = d_MN * d_MS_inv;
            T d2 = d_MP * d_MS_inv;
            T d3 = d_MNP * d_MS_inv;

            //5-2. compute the matrix [u1'; u2'; u3']. Flops: 21
            T u1[3] = {d2 * MN[0] - d3 * MP[0],
                       d2 * MN[1] - d3 * MP[1],
                       d2 * MN[2] - d3 * MP[2]};
            T u2[3] = {d1 * MP[0] - d3 * MN[0],
                       d1 * MP[1] - d3 * MN[1],
                       d1 * MP[2] - d3 * MN[2]};
            T u3[3] = {d_MS_inv * MS[0], d_MS_inv * MS[1], d_MS_inv * MS[2]};

            for (int ii = 0; ii < solu_size; ++ii) {

                //5-3. Compute lam_M. Flops: 11 + 1*sqrt
                T lam_1 = lam1_solu[ii];
                T lam_2 = lam2_solu[ii];
                T lam_M = sqrt(d_MN / ((l_NN * lam_1 - 2 * l_MN) * lam_1 + l_MM));

                T lam_N = lam_1 * lam_M;
                T lam_P = lam_2 * lam_M;

                T lam_Me_M20 = lam_M * e_M2[0];
                T lam_Me_M21 = lam_M * e_M2[1];
                T lam_Me_M22 = lam_M * e_M2[2];

                //5-4. Compute the first matrix. Flops: 27
                T MN_C[3] = {
                        lam_N * e_N2[0] - lam_Me_M20,
                        lam_N * e_N2[1] - lam_Me_M21,
                        lam_N * e_N2[2] - lam_Me_M22
                };

                T MP_C[3] = {
                        lam_P * e_P2[0] - lam_Me_M20,
                        lam_P * e_P2[1] - lam_Me_M21,
                        lam_P * e_P2[2] - lam_Me_M22
                };
                T MS_C[3] = {
                        MN_C[1] * MP_C[2] - MN_C[2] * MP_C[1],
                        MN_C[2] * MP_C[0] - MN_C[0] * MP_C[2],
                        MN_C[0] * MP_C[1] - MN_C[1] * MP_C[0]
                };

                // 5-5. Compute R, t. Flops: 66
                Eigen::Matrix<T, 3, 3> M1, M2;
                M1 << MN_C[0], MP_C[0], MS_C[0],
                      MN_C[1], MP_C[1], MS_C[1],
                      MN_C[2], MP_C[2],MS_C[2];

                M2 << u1[0], u1[1], u1[2],
                      u2[0], u2[1], u2[2],
                      u3[0], u3[1], u3[2];

                R_ = M1 * M2;
                t_[0] = lam_M * e_M2[0] -
                        (R_(0, 0) * worldPoints[0] + R_(0, 1) * worldPoints[4] + R_(0, 2) * worldPoints[8]);
                t_[1] = lam_M * e_M2[1] -
                        (R_(1, 0) * worldPoints[0] + R_(1, 1) * worldPoints[4] + R_(1, 2) * worldPoints[8]);
                t_[2] = lam_M * e_M2[2] -
                        (R_(2, 0) * worldPoints[0] + R_(2, 1) * worldPoints[4] + R_(2, 2) * worldPoints[8]);

                Ts[ii] = t_;
                Rs[ii] = R_;
            }
        }
        return solu_size;
    }
}
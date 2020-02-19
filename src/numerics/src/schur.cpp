/*********************************************************************
 *
 *  Software License Agreement
 *
 *  Copyright (c) 2020,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Authors: Christoph Rösmann
 *********************************************************************/

#include <corbo-numerics/schur.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>

#include <corbo-core/console.h>
#include <corbo-core/utilities.h>
#include <corbo-numerics/sylvester_continuous.h>  // for swap_blocks_schur_form()

#include <Eigen/Eigenvalues>
#include <Eigen/Householder>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <limits>

namespace corbo {

void schur_decomposition_2d(Eigen::Ref<Eigen::Matrix2d> T, Eigen::Ref<Eigen::Matrix2d> U)
{
    assert(U.rows() == 2 && U.cols() == 2);
    assert(have_equal_size(T, U));

    double c = T(1, 0);

    if (approx_zero(c, 1e-30))
    {
        U.setIdentity();
        return;
    }

    double a = T(0, 0);
    double b = T(0, 1);
    double d = T(1, 1);

    if (approx_zero(b, 1e-30))
    {
        // swap rows and columns
        T(0, 0) = d;
        T(0, 1) = -c;
        T(1, 0) = 0;
        T(1, 1) = a;

        // 2d rotation matrix with cos=0 and sin=1
        U(0, 0) = 0;
        U(0, 1) = -1;
        U(1, 0) = 1;
        U(1, 1) = 0;
        return;
    }

    if (approx_zero(a - d, 1e-30) && util::sign(1, b) != util::sign(1, c))
    {
        U.setIdentity();
        return;
    }

    double a_d    = a - d;
    double p      = 0.5 * a_d;
    double bc_max = std::max(std::abs(b), std::abs(c));
    double bc_mis = std::min(std::abs(b), std::abs(c)) * util::sign(1, b) * util::sign(1, c);
    double scale  = std::max(std::abs(p), bc_max);
    double z      = (p / scale) * p + (bc_max / scale) * bc_mis;

    if (z > 5 * std::numeric_limits<double>::min())
    {
        // real eigenvalues, we need elements a and d
        z = p + util::sign(std::sqrt(scale) * std::sqrt(z), p);
        T(0, 0) = d + z;
        T(1, 1) = d - (bc_max / z) * bc_mis;
        // compute b and U
        T(0, 1) = b - c;
        T(1, 0) = 0;

        double tau = util::l2_norm_2d(c, z);
        U(0, 0) = z / tau;
        U(0, 1) = -c / tau;
        U(1, 0) = c / tau;
        U(1, 1) = z / tau;
    }
    else
    {
        // complex eigenvalues (or real (almost) equal)
        // start with making diagonal elements equal
        double sigma = b + c;
        double tau   = util::l2_norm_2d(sigma, a_d);
        double u_cos = std::sqrt(0.5 * (1 + std::abs(sigma) / tau));
        double u_sin = -(p / (tau * u_cos)) * util::sign(1, sigma);
        U(0, 0) = u_cos;
        U(0, 1) = -u_sin;
        U(1, 0) = u_sin;
        U(1, 1) = u_cos;

        T = U.transpose() * T * U;

        double temp = 0.5 * (T(0, 0) + T(1, 1));  // symmetry
        T(0, 0) = temp;
        T(1, 1) = temp;

        a = T(0, 0);
        b = T(0, 1);
        c = T(1, 0);
        d = T(1, 1);

        if (!approx_zero(c, 1e-30))  // c!= 0
        {
            if (!approx_zero(b, 1e-30))  // b != 0
            {
                if (util::sign(1, b) == util::sign(1, c))  // sign(b) == sign(c)
                {
                    // real eigenvalues -> just reduce to upper triangular form
                    double sab = std::sqrt(std::abs(b));
                    double sac = std::sqrt(std::abs(c));
                    p          = util::sign(sab * sac, c);
                    tau        = 1.0 / std::sqrt(std::abs(b + c));
                    T(0, 0) = temp + p;
                    T(1, 1) = temp - p;
                    T(0, 1) = b - c;
                    T(1, 0) = 0;

                    double cos1 = sab * tau;
                    double sin1 = sac * tau;

                    temp = U(0, 0) * cos1 - U(1, 0) * sin1;
                    U(1, 0) = U(0, 0) * sin1 + U(1, 0) * cos1;
                    U(0, 1) = -U(1, 0);
                    U(0, 0) = temp;
                    U(1, 1) = temp;
                }
            }
            else
            {
                T(0, 1) = -b;
                T(1, 0) = 0;
                temp = U(0, 0);
                U(0, 0) = U(0, 1);
                U(1, 1) = U(0, 1);
                U(0, 1) = -temp;
                U(1, 0) = temp;
            }
        }
    }
}

bool swap_schur_blocks(Eigen::Ref<Eigen::MatrixXd> T, int ra11, int p, int q, Eigen::Ref<Eigen::MatrixXd> Q, bool standardize)
{
    assert(is_square(T));
    assert(have_equal_size(T, Q));
    assert(p > 0 && p < 3);
    assert(q > 0 && q < 3);
    assert(ra11 + p + q <= T.rows());

    bool ret_val = true;

    Eigen::MatrixXd BACKUP = T;

    const int na = p + q;

    // make sure that both eigenvalues differ,
    // otherwise the sylvester equation does not have a unique solution
    if (na == 2)
    {
        if (essentially_equal(T(ra11, ra11), T(ra11 + 1, ra11 + 1), 1e-30)) return true;  // same eigenvalues
    }

    int ra22 = ra11 + p;

    if (na == 4)
    {
        // check if eigen values are equal
        if (essentially_equal(T(ra11, ra11), T(ra22, ra22), 1e-30))  // real part
        {
            double imag1 = std::sqrt(std::abs(T(ra11, ra11 + 1) * T(ra11 + 1, ra11)));
            double imag2 = std::sqrt(std::abs(T(ra22, ra22 + 1) * T(ra22 + 1, ra22)));
            if (essentially_equal(imag1, imag2, 1e-30)) return true;  // same eigenvalues
        }
    }

    // Solve sylvester equation A11 X - X A22 = \gamma A12

    // [1] suggests to use a small value close to machine precision times the norm of the matrix
    // in cases if there is are small diagonal values while solving.
    const double gamma = 1e3 * std::numeric_limits<double>::epsilon() * T.block(ra11, ra11, na, na).norm();

    Eigen::MatrixXd X;
    if (!SylvesterContinuous::solve(T.block(ra11, ra11, p, p), -T.block(ra22, ra22, q, q), -gamma * T.block(ra11, ra22, p, q), X)) return false;

    Eigen::MatrixXd G(p + q, q);
    G.block(0, 0, p, q) = -X;
    G.block(p, 0, q, q).setIdentity();
    G.block(p, 0, q, q) *= gamma;

    // store norm for later stability check
    double A_infnorm = T.block(ra11, ra11, na, na).lpNorm<Eigen::Infinity>();

    Eigen::HouseholderQR<Eigen::MatrixXd> qr_dec(G);  // TODO(roesmann) FullPivHouseholder?

    // perform swap on complete A
    // this code might be used with a different temporary matrix in order to avoid destruction in T
    // if the stability criterion fails:
    // T.block(ra11, ra11, p + q, T.rows() - ra11).applyOnTheLeft(qr_dec.householderQ().transpose());
    // T.block(0, ra11, na, na).applyOnTheRight(qr_dec.householderQ());
    // Q.block(ra11, ra11, na, na) = qr_dec.householderQ();

    // perform swap on complete T
    T.block(ra11, ra11, na, T.cols() - ra11).applyOnTheLeft(qr_dec.householderQ().transpose());
    T.block(0, ra11, ra22 + q, na).applyOnTheRight(qr_dec.householderQ());
    Q.block(0, ra11, Q.rows(), na).applyOnTheRight(qr_dec.householderQ());

    // a11 and a22 are now swapped
    ra22 = ra11 + q;

    constexpr const double threshold = std::numeric_limits<double>::epsilon() * 2;
    const double stability_crit      = 10.0 * threshold * A_infnorm;  // see [1]
    double A21_norm                  = T.block(ra22, ra11, p, q).lpNorm<Eigen::Infinity>();
    if (A21_norm < stability_crit)
    {
        // set lower left part to zero
        T.block(ra22, ra11, p, q).setZero();  // this is now p x q due to swap
    }
    else
    {
        ret_val = false;
        PRINT_WARNING("swap_blocks_schur_form(): numerical instability detected: " << A21_norm << "/" << stability_crit << ",\nA21:\n"
                                                                                   << T.block(ra22, ra11, q, p));
    }

    // standardize complex blocks if desired
    if (standardize)
    {
        if (q == 2)
        {  // this block is now in the top left corner
            Eigen::Matrix2d U;
            schur_decomposition_2d(T.block(ra11, ra11, q, q), U);

            // also update T from left to right starting right of the current block
            T.block(ra11, ra22, q, T.cols() - ra22).applyOnTheLeft(U.transpose());

            // also update T from up to down ending above the current block
            for (int i = 0; i < ra11; ++i)
            {
                // We use a map in following to treat 2d row vector as 2d column vector without transposing etc.
                // With the current formulation, we can only handle a single row at once (hence the for-loop).
                Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<>> map(T.col(ra11).row(i).data(), 2, 1, Eigen::InnerStride<>(T.outerStride()));
                map.applyOnTheLeft(U.transpose());
            }

            // update unitary matrix Q
            for (int i = 0; i < Q.rows(); ++i)
            {
                Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<>> map(Q.col(ra11).row(i).data(), 2, 1, Eigen::InnerStride<>(Q.outerStride()));
                map.applyOnTheLeft(U.transpose());
            }
        }

        // transform new 2 x 2 block to standard-form
        if (p == 2)
        {  // this block is now in the bottom right corner
            Eigen::Matrix2d U;
            schur_decomposition_2d(T.block(ra22, ra22, p, p), U);

            // also update T from left to right starting right of the current block
            if (T.cols() > ra22 + p) T.block(ra22, ra22 + p, p, T.cols() - ra22 - p).applyOnTheLeft(U.transpose());

            // also update T from up to down ending above the current block
            for (int i = 0; i < ra22; ++i)
            {
                // We use a map in following to treat 2d row vector as 2d column vector without transposing etc.
                // With the current formulation, we can only handle a single row at once (hence the for-loop).
                Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<>> map(T.col(ra22).row(i).data(), 2, 1, Eigen::InnerStride<>(T.outerStride()));
                map.applyOnTheLeft(U.transpose());
            }

            // update unitary matrix Q
            for (int i = 0; i < Q.rows(); ++i)
            {
                Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<>> map(Q.col(ra22).row(i).data(), 2, 1, Eigen::InnerStride<>(Q.outerStride()));
                map.applyOnTheLeft(U.transpose());
            }
        }
    }
    return ret_val;
}

}  // namespace corbo

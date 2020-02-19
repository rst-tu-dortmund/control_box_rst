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

#include <corbo-numerics/lyapunov_discrete.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace corbo {

bool LyapunovDiscrete::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& Q, Eigen::MatrixXd& X)
{
    assert(is_square(A));
    assert(have_equal_size(A, Q));

    Eigen::ComplexSchur<Eigen::MatrixXd> A_schur(A, true);

    if (A_schur.info() == Eigen::NoConvergence) return false;

    const Eigen::MatrixXcd& T = A_schur.matrixT();
    const Eigen::MatrixXcd& U = A_schur.matrixU();

    // transform rhs
    Eigen::MatrixXcd F = U.adjoint() * Q * U;

    const int n = A.rows();
    Eigen::MatrixXcd Y(n, n);
    Y.setZero();

    Eigen::JacobiSVD<Eigen::MatrixXcd> lhs_svd(n, n, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXcd lhs(n, n);
    lhs.setIdentity();

    Eigen::MatrixXcd T_adjoint = T.adjoint();

    // solve T Y T^T - Y + F = 0
    // construct transformed solution by backwards substitution
    for (int k = n - 1; k >= 0; --k)
    {
        Eigen::VectorXcd rhs = F.col(k) + T * Y * T_adjoint.col(k);
        // move non-zero diagonal elements to lhs and multiply by T
        // (this is the non-zero part of T Y T^T which depends only on Y.col(k)).
        const std::complex<double> factor = T_adjoint(k, k);
        lhs -= T * factor;
        // compute SVD
        // TODO(roesmann): here is place for improving efficiency
        //                 and robustness by avoiding the SVD (see references).
        lhs_svd.compute(lhs);
        // solve linear system
        Y.col(k) = lhs_svd.solve(rhs);
        // revert rhs diagonal element for subsequent iteration
        if (k > 0) lhs.setIdentity();
    }

    // transform back
    X = (U * Y * U.adjoint()).real();

    return true;
}

bool LyapunovDiscrete::hasUniqueSolution(const Eigen::Ref<const Eigen::MatrixXd>& A)
{
    if (!is_square(A)) return false;

    Eigen::VectorXcd eigvals_A = A.eigenvalues();

    for (int i = 0; i < eigvals_A.size(); ++i)
    {
        for (int j = i; j < eigvals_A.size(); ++j)
        {
            std::complex<double> product = eigvals_A[i] * eigvals_A[j];
            if (approx_equal(product.real(), 1) && approx_equal(product.imag(), 0))  // TODO(roesman) is 1+i0 correct, or unit circle in general?
                return false;
        }
    }
    return true;
}

}  // namespace corbo

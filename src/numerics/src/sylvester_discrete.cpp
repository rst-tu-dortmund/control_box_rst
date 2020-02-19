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

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>
#include <corbo-numerics/sylvester_discrete.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace corbo {

bool SylvesterDiscrete::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                              const Eigen::Ref<const Eigen::MatrixXd>& C, Eigen::MatrixXd& X)
{
    assert(is_square(A));
    assert(is_square(B));
    assert(A.rows() == C.rows());
    assert(B.rows() == C.cols());

    Eigen::ComplexSchur<Eigen::MatrixXd> A_schur(A, true);
    Eigen::ComplexSchur<Eigen::MatrixXd> B_schur(B, true);

    if (A_schur.info() == Eigen::NoConvergence) return false;
    if (B_schur.info() == Eigen::NoConvergence) return false;

    const Eigen::MatrixXcd& T_A = A_schur.matrixT();
    const Eigen::MatrixXcd& U_A = A_schur.matrixU();
    const Eigen::MatrixXcd& T_B = B_schur.matrixT();
    const Eigen::MatrixXcd& U_B = B_schur.matrixU();

    const int n = C.rows();
    const int m = C.cols();

    // transform rhs
    Eigen::MatrixXcd F = U_A.adjoint() * C * U_B;

    Eigen::MatrixXcd Y(n, m);
    Y.setZero();

    Eigen::JacobiSVD<Eigen::MatrixXcd> lhs_svd(n, n, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXcd lhs(n, m);
    lhs.setIdentity();

    // solve T_A Y T_B - Y + F = 0
    // construct transformed solution by forward substitution
    const int mm = m - 1;
    for (int k = 0; k < m; ++k)
    {
        Eigen::VectorXcd rhs = F.col(k) + T_A * Y * T_B.col(k);
        // move non-zero diagonal elements to lhs and multiply by T
        // (this is the non-zero part of T_A Y T_B which depends only on Y.col(k)).
        const std::complex<double> factor = T_B(k, k);
        lhs -= T_A * factor;
        // compute SVD
        // TODO(roesmann): here is place for improving efficiency
        //                 and robustness by avoiding the SVD (see references).
        lhs_svd.compute(lhs);
        // solve linear system
        Y.col(k) = lhs_svd.solve(rhs);
        // revert rhs diagonal element for subsequent iteration
        if (k < mm) lhs.setIdentity();
    }

    // transform back
    X = (U_A * Y * U_B.adjoint()).real();

    return true;
}

bool SylvesterDiscrete::hasUniqueSolution(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B)
{
    if (!is_square(A) || !is_square(B)) return false;

    Eigen::VectorXcd eigvals_A = A.eigenvalues();
    Eigen::VectorXcd eigvals_B = B.eigenvalues();

    for (int i = 0; i < eigvals_A.size(); ++i)
    {
        for (int j = 0; j < eigvals_B.size(); ++j)
        {
            std::complex<double> product = eigvals_A[i] * eigvals_B[j];
            if (approx_equal(product.real(), 1) && approx_equal(product.imag(), 0))  // TODO(roesman) is 1+i0 correct, or unit circle in general?
                return false;
        }
    }
    return true;
}

}  // namespace corbo

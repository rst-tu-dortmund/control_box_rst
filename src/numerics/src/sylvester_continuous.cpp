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

#include <corbo-numerics/sylvester_continuous.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace corbo {

bool SylvesterContinuous::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
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

    Eigen::MatrixXcd T_A        = A_schur.matrixT();  // create copy since we need to modify the diagonal later
    const Eigen::MatrixXcd& U_A = A_schur.matrixU();
    const Eigen::MatrixXcd& T_B = B_schur.matrixT();
    const Eigen::MatrixXcd& U_B = B_schur.matrixU();

    const int n = C.rows();
    const int m = C.cols();

    // transform rhs
    Eigen::MatrixXcd F = U_A.adjoint() * C * U_B;

    Eigen::MatrixXcd Y(n, m);
    Y.setZero();

    Eigen::JacobiSVD<Eigen::MatrixXcd> T_A_svd(n, n, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // solve T_A Y + Y T_B + F = 0
    // construct transformed solution by forward substitution
    const int mm = m - 1;
    for (int k = 0; k < m; ++k)
    {
        Eigen::VectorXcd rhs = F.col(k) + Y * T_B.col(k);
        // move non-zero diagonal element to lhs
        const std::complex<double> offset = T_B(k, k);
        T_A.diagonal().array() += offset;
        // compute SVD
        // TODO(roesmann): here is place for improving efficiency
        //                 and robustness by avoiding the SVD (see references).
        T_A_svd.compute(T_A);
        // solve linear system
        Y.col(k) = T_A_svd.solve(-rhs);
        // revert offset
        if (k < mm) T_A.diagonal().array() -= offset;
    }

    // transform back
    X = (U_A * Y * U_B.adjoint()).real();

    return true;
}

bool SylvesterContinuous::hasUniqueSolution(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B)
{
    if (!is_square(A) || !is_square(B)) return false;

    Eigen::VectorXcd eigvals_A = A.eigenvalues();
    Eigen::VectorXcd eigvals_B = B.eigenvalues();

    for (int i = 0; i < eigvals_A.size(); ++i)
    {
        for (int j = 0; j < eigvals_B.size(); ++j)
        {
            if (approx_equal(eigvals_A[i].real(), -eigvals_B[j].real()) && approx_equal(eigvals_A[i].imag(), -eigvals_B[j].imag())) return false;
        }
    }
    return true;
}

}  // namespace corbo

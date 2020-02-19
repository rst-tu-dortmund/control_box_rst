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

#include <corbo-numerics/lyapunov_continuous.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace corbo {

bool LyapunovContinuous::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& Q, Eigen::MatrixXd& X)
{
    assert(is_square(A));
    assert(have_equal_size(A, Q));

    Eigen::ComplexSchur<Eigen::MatrixXd> A_schur(A, true);

    if (A_schur.info() == Eigen::NoConvergence) return false;

    Eigen::MatrixXcd T        = A_schur.matrixT();  // create copy since we need to modify the diagonal later
    const Eigen::MatrixXcd& U = A_schur.matrixU();

    // transform rhs
    Eigen::MatrixXcd F = U.adjoint() * Q * U;

    const int n = A.rows();
    Eigen::MatrixXcd Y(n, n);
    Y.setZero();

    Eigen::JacobiSVD<Eigen::MatrixXcd> T_svd(n, n, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXcd T_adjoint = T.adjoint();

    // solve T Y + Y T^* + F = 0
    // construct transformed solution by backwards substitution
    for (int k = n - 1; k >= 0; --k)
    {
        Eigen::VectorXcd rhs = F.col(k) + Y * T_adjoint.col(k);
        // move non-zero diagonal element to lhs
        const std::complex<double> offset = T_adjoint(k, k);
        T.diagonal().array() += offset;
        // compute SVD
        // TODO(roesmann): here is place for improving efficiency
        //                 and robustness by avoiding the SVD (see references).
        T_svd.compute(T);
        // solve linear system
        Y.col(k) = T_svd.solve(-rhs);
        // revert rhs diagonal element for subsequent iteration
        if (k > 0) T.diagonal().array() -= offset;
    }

    // transform back
    X = (U * Y * U.adjoint()).real();

    return true;
}

bool LyapunovContinuous::hasUniqueSolution(const Eigen::Ref<const Eigen::MatrixXd>& A)
{
    if (!is_square(A)) return false;

    Eigen::VectorXcd eigvals_A = A.eigenvalues();

    for (int i = 0; i < eigvals_A.size(); ++i)
    {
        for (int j = i; j < eigvals_A.size(); ++j)
        {
            if (approx_equal(eigvals_A[i].real(), -eigvals_A[j].real()) && approx_equal(eigvals_A[i].imag(), -eigvals_A[j].imag())) return false;
        }
    }
    return true;
}

}  // namespace corbo

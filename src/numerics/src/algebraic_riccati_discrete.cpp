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

#include <corbo-numerics/algebraic_riccati_discrete.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>
#include <corbo-numerics/schur.h>

#include <corbo-core/console.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>

namespace corbo {

bool AlgebraicRiccatiDiscrete::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                     const Eigen::Ref<const Eigen::MatrixXd>& Q, const Eigen::Ref<const Eigen::MatrixXd>& R, Eigen::MatrixXd& X,
                                     Eigen::MatrixXd* G)
{
    assert(is_square(A));
    assert(have_equal_size(A, Q));
    assert(A.rows() == B.rows());
    assert(R.rows() == B.cols());
    assert(is_square(R));
    const int p = A.rows();  // dim_states

    Eigen::LDLT<Eigen::MatrixXd> R_cholesky(R);
    Eigen::MatrixXd R_inv_B_t   = R_cholesky.solve(B.transpose());
    Eigen::MatrixXd A_inv_trans = A.inverse().transpose();  // A must be invertible

    // Construct Hamiltonian matrix
    // H = [ A+B R^-1 B^T (A^-1)^T Q      -B R^-1 B^T (A^-1)^T]
    //     [-(A^-1)^T Q                     (A^-1)^T       ]
    Eigen::MatrixXd H(2 * p, 2 * p);
    H << A + B * R_inv_B_t * A_inv_trans * Q, -B * R_inv_B_t * A_inv_trans, -A_inv_trans * Q, A_inv_trans;

    bool success = solveRiccatiHamiltonianSchur(H, X);

    if (G && success)
    {
        Eigen::LDLT<Eigen::MatrixXd> G_cholesky(R + B.transpose() * X * B);
        *G = G_cholesky.solve(B.transpose()) * X * A;
    }
    return success;
}

bool AlgebraicRiccatiDiscrete::solveRiccatiHamiltonianSchur(const Eigen::MatrixXd& H, Eigen::MatrixXd& X)
{
    assert(is_square(H));
    assert(H.rows() % 2 == 0);

    // perform schur to find an orthogonal transformation that reduces H to real Schur form
    Eigen::RealSchur<Eigen::MatrixXd> schur(H);

    Eigen::MatrixXd T = schur.matrixT();
    Eigen::MatrixXd U = schur.matrixU();

    // reorder Schur form such that all negative eigenvalues are located on the upper left blocks.
    int subspace_dim;
    bool success = reorder_schur_blocks(T, U, AlgebraicRiccatiDiscrete::isInsideUnitCircle, &subspace_dim, false);

    const int p = H.rows() / 2;

    if (!success || p != subspace_dim) return false;

    // U = [U1; U2
    //      U2; U3]

    // solution X can be obtained by X = U2 * U1^-1
    // since U1.inverse() is not recomended due to numerical effects
    // we must rewrite the equation in order to facilitate solving with Eigen:
    // X^T = (U2 * U1^-1)^T = (U1^T)^(-1) U2^T
    Eigen::JacobiSVD<Eigen::MatrixXd> U1_svd(U.block(0, 0, p, p).transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    X = U1_svd.solve(U.block(p, 0, p, p).transpose());  // X is symmetric, hence no final transpose
    return true;
}

bool AlgebraicRiccatiDiscrete::isClosedLoopStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& G)
{
    assert(is_square(A));
    assert(A.rows() == B.rows());
    assert(A.rows() == G.cols());
    assert(G.rows() == B.cols());

    Eigen::VectorXcd eig = (A - B * G).eigenvalues();
    for (int i = 0; i < eig.size(); ++i)
    {
        if (std::abs(eig[i]) >= 1) return false;
    }
    return true;
}

}  // namespace corbo

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

#include <corbo-numerics/algebraic_riccati_continuous.h>

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>
#include <corbo-numerics/schur.h>

#include <corbo-core/console.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>

namespace corbo {

bool AlgebraicRiccatiContinuous::solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
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
    Eigen::MatrixXd R_inv_B_t = R_cholesky.solve(B.transpose());

    // Construct Hamiltonian matrix
    // H = [ A    -B R^-1 B^T  ]
    //     [-Q      -A^T       ]
    Eigen::MatrixXd H(2 * p, 2 * p);
    H << A, -B * R_inv_B_t, -Q, -A.transpose();

#ifndef NDEBUG
    PRINT_WARNING_COND_NAMED(hasRealPartsCloseToZero(H),
                             "Eigenvalues of the Hamiltonian are close to the imaginary axis. Numerically unstable results are expected.");
#endif

    bool success = solveRiccatiHamiltonianSchur(H, X);

    if (G && success) *G = R_inv_B_t * X;
    return success;
}

bool AlgebraicRiccatiContinuous::solveRiccatiHamiltonianSchur(const Eigen::MatrixXd& H, Eigen::MatrixXd& X)
{
    assert(is_square(H));
    assert(H.rows() % 2 == 0);

    // perform schur to find an orthogonal transformation that reduces H to real Schur form
    Eigen::RealSchur<Eigen::MatrixXd> schur(H);

    Eigen::MatrixXd T = schur.matrixT();
    Eigen::MatrixXd U = schur.matrixU();

    // reorder Schur form such that all negative eigenvalues are located on the upper left blocks.
    int subspace_dim;
    bool success = reorder_schur_blocks(T, U, AlgebraicRiccatiContinuous::hasNegativeRealPart, &subspace_dim, false);

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

bool AlgebraicRiccatiContinuous::hasRealPartsCloseToZero(const Eigen::Ref<const Eigen::MatrixXd>& matrix)
{
    Eigen::VectorXcd eig = matrix.eigenvalues();
    for (int i = 0; i < eig.size(); ++i)
    {
        if (std::abs(eig[i].real()) <= 1e-12)
        {
            return true;
        }
    }
    return false;
}

/*
void AlgebraicRiccatiContinuous::sovleRiccatiHamiltonianMatrixSign(Eigen::MatrixXd& Z, Eigen::MatrixXd& X)
{
    // TODO(roesmann): WE HAVE CHANGED THE ORDER OF ELEMENTS IN THE HAMILTONIAN, but we haven't
    //                 applied it in this method.

    assert(is_square(Z));
    assert(Z.rows() % 2 == 0);

    const int p = Z.rows() / 2;

    Eigen::MatrixXd Z_old;

    const double tolerance      = 1e-9;
    const double max_iterations = 100;

    double relative_norm  = CORBO_INF_DBL;
    std::size_t iteration = 0;

    double dp = double(2 * p);

    do
    {
        Z_old = Z;
        // R. Byers. Solving the algebraic Riccati equation with the matrix sign function. Linear Algebra Appl., 85:267–279, 1987
        // Added determinant scaling to improve convergence (converges in rough half the iterations with this)
        double ck = std::pow(std::abs(Z.determinant()), -1.0 / dp);   // reciprocal of the geometric mean of the eigenvalues
        // TODO(roesmann): here is place for improving efficiency: compute Z.determinant() as a by-product of Z^-1
        Z *= ck;
        Z             = Z - 0.5 * (Z - Z.inverse());
        relative_norm = (Z - Z_old).norm();
        ++iteration;
    } while (iteration < max_iterations && relative_norm > tolerance);

    //    Eigen::MatrixXd W11 = Z.block(0, 0, p, p); // -> A
    //    Eigen::MatrixXd W12 = Z.block(0, p, p, p); // -> Q
    //    Eigen::MatrixXd W21 = Z.block(p, 0, p, p); // -> B R^-1 B^T
    //    Eigen::MatrixXd W22 = Z.block(p, p, p, p); // -> -A^T

    Eigen::MatrixXd lhs(2 * p, p);
    Eigen::MatrixXd rhs(2 * p, p);
    // Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(p, p);
    //     lhs << W12, W22 + eye;
    //     rhs << W11 + eye, W21;
    // lhs << Z.block(0, p, p, p), Z.block(p, p, p, p) + Eigen::MatrixXd::Identity(p, p);
    // rhs << Z.block(0, 0, p, p) + Eigen::MatrixXd::Identity(p, p), Z.block(p, 0, p, p);

    //     lhs << W11 - eye, -W12;
    //     rhs << W21, eye - W22;
    lhs << Z.block(0, 0, p, p) - Eigen::MatrixXd::Identity(p, p), -Z.block(0, p, p, p);
    rhs << Z.block(p, 0, p, p), Eigen::MatrixXd::Identity(p, p) - Z.block(p, p, p, p);

    // TODO(roesmann) the identity matrix seems to be redundant but we need it to get suitable results (maybe for numerical stability in
    // case A==0?)

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(lhs, Eigen::ComputeThinU | Eigen::ComputeThinV);

    X = svd.solve(rhs);

    // exploit symmetry to reduce numerical error
    X = (X + X.transpose()) * 0.5;
}
*/

bool AlgebraicRiccatiContinuous::isClosedLoopStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                                    const Eigen::Ref<const Eigen::MatrixXd>& G)
{
    assert(is_square(A));
    assert(A.rows() == B.rows());
    assert(A.rows() == G.cols());
    assert(G.rows() == B.cols());

    Eigen::VectorXcd eig = (A - B * G).eigenvalues();
    for (int i = 0; i < eig.size(); ++i)
    {
        if (eig[i].real() >= 0) return false;
    }
    return true;
}

bool AlgebraicRiccatiContinuous::isNumericallyStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                                     const Eigen::Ref<const Eigen::MatrixXd>& Q, const Eigen::Ref<const Eigen::MatrixXd>& R)
{
    assert(is_square(A));
    assert(have_equal_size(A, Q));
    assert(A.rows() == B.rows());
    assert(R.rows() == B.cols());
    assert(is_square(R));

    const int p = A.rows();  // dim_states

    Eigen::LDLT<Eigen::MatrixXd> R_cholesky(R);
    Eigen::MatrixXd R_inv_B_t = R_cholesky.solve(B.transpose());

    Eigen::MatrixXd H(2 * p, 2 * p);
    H << A, -B * R_inv_B_t, -Q, -A.transpose();
    return !hasRealPartsCloseToZero(H);
}

}  // namespace corbo

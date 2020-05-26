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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_CONTINUOUS_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_CONTINUOUS_H_

#include <corbo-core/types.h>

#include <Eigen/Eigenvalues>

namespace corbo {

/**
 * @brief Methods for dealing with continuous-time algebraic Riccati equations
 *
 * @ingroup numerics
 *
 * The continuous-time algebraic Riccati equation is given by
 * \f[
 *       A^T X + X A - X B R^{-1} B^T X + Q = 0
 * \f]
 *
 * \f$ A \f$ and \f$ B \f$ must be stabilizable.
 *
 * The resulting gain matrix is given by
 * \f[ G = R^{-1} B^T X \f]
 *
 * @see AlgebraicRiccatiDiscrete LyapunovDiscrete LyapunovContinuous
 *      SylvesterContinuous SylvesterDiscrete
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Add a method that determines if sufficient conditions hold and that (A,B) is stabilizable.
 * @todo Add support for the more generalized Riccati equation (augmented by matrix S and E)
 *       E.g. Hamiltonian for augmentation with S: http://web.dcsc.tudelft.nl/~cscherer/books/thesis.pdf
 */
class AlgebraicRiccatiContinuous
{
 public:
    /**
     * @brief Solve continuous-time algebraic Riccati equation
     *
     * Solve \f$ A^T X + X A - X B R^{-1} B^T X + Q = 0 \f$
     * w.r.t. \f$ X \f$.
     *
     * \f$ A \f$ and \f$ B \f$ must be stabilizable (all unstable modes are controllable).
     *  In addition, the associated Hamiltonian matrix must have no eigenvalue on the imaginary axis.
     *  Sufficient conditions: R > 0.
     *
     * This method optionally returns the gain matrix \f$ G = R^{-1} B^T X \f$.
     *
     * @warning Input matrix sizes are checked only in debug mode and occuring issues immediately break execution
     *
     * @param[in]  A    [n x n] matrix
     * @param[in]  B    [n x m] matrix
     * @param[in]  Q    [n x n] matrix
     * @param[in]  R    [m x m] matrix
     * @param[out] X    Solution of the algebraic Riccati equation with size [n x n] (must not be preallocated).
     * @param[out] G    Gain matrix G = R^{-1} B^T X with size [m x n] (optional)
     * @returns true if a finite solution is found, false otherwise (but no NaN checking is performend on X, use X.allFinite()).
     */
    static bool solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                      const Eigen::Ref<const Eigen::MatrixXd>& Q, const Eigen::Ref<const Eigen::MatrixXd>& R, Eigen::MatrixXd& X,
                      Eigen::MatrixXd* G = nullptr);

    /**
     * @brief Determine if the closed-loop system is stable
     *
     * The closed-loop system is stable if the closed-loop eigenvalues
     * exhibit only negative real parts.
     * This method returns \c true if \f$ A - B G \f$ has only eigenvalues with negative real parts.
     * Hereby, G denotes the gain matrix obtained from solve().
     *
     * @param[in]  A    Open-loop system matrix [n x n]
     * @param[in]  B    Input matrix [n x m]
     * @param[in]  G    Feedback gain matrix according to the solution of the algebraic Riccati equation [m x n]
     * @return true if the closed-loop system is stable, false otherwise.
     */
    static bool isClosedLoopStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                   const Eigen::Ref<const Eigen::MatrixXd>& G);

    /**
     * @brief Determine if solving the riccati equation seems to be numerically stable
     *
     * However, it is not guaranteed...
     *
     * This method checks the real parts of the egienvalues of the Hamiltonian are close to the imaginary axis.
     * If yet, numerically stability is likely to be not satisfactory (also see Matlabs care solver).
     *
     * @param[in]  A    [n x n] matrix
     * @param[in]  B    [n x m] matrix
     * @param[in]  Q    [n x n] matrix
     * @param[in]  R    [m x m] matrix
     * @return true if the problem seems to be numerically stable
     */
    static bool isNumericallyStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                    const Eigen::Ref<const Eigen::MatrixXd>& Q, const Eigen::Ref<const Eigen::MatrixXd>& R);

 protected:
    /**
     * @brief Solve Hamiltonian via Schur method
     *
     * The Hamiltonian constructed in solve() is
     * solved via the Schur decomposition.
     * The Real Schur form of A is reordered such that blocks
     * on the diagonal with eigenvalues having a negative real part
     * precede other one.
     * The corresponding column in the orthogonal transformation matrix U
     * spans the invariant subspace: [U1 U2]^T.
     * The solution of the Riccati equation is then obtained by
     * \f[
     *    X = U_2 U_1^{-1}
     * \f]
     *
     * The implementation is based on [1].
     *
     * [1] A. J. Laub. "A Schur method for solving algebraic Riccati equations",
     *     IEEE Conference on Decision and Control, 1978.
     *
     * @param[in]    H   Hamiltonian matrix of the Riccati equation [square]
     * @param[out]   X   Resuling solution: symmetric matrix with size(H).
     * @return true if solving was successful (but no NaN checking is performend on X, use X.allFinite()).
     */
    static bool solveRiccatiHamiltonianSchur(const Eigen::MatrixXd& H, Eigen::MatrixXd& X);

    // matrix sign method, currently not in use and needs to be updated for future use
    // [1] R. Byers. "Solving the algebraic Riccati equation with the matrix sign function".
    //     Linear ALgebra and its Applications. Volume 85, 267-278, 1987.
    // static void sovleRiccatiHamiltonianMatrixSign(Eigen::MatrixXd& Z, Eigen::MatrixXd& X);

    /// Predicate for Schur block ordering (see solveRiccatiHamiltonianSchur())
    static bool hasNegativeRealPart(const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; }

    /// Check if the real parts of the provided matrix are close to zero
    static bool hasRealPartsCloseToZero(const Eigen::Ref<const Eigen::MatrixXd>& matrix);
};

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_CONTINUOUS_H_

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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_DISCRETE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_DISCRETE_H_

#include <corbo-core/types.h>

#include <Eigen/Eigenvalues>

namespace corbo {

/**
 * @brief Methods for dealing with discrete-time algebraic Riccati equations
 *
 * @ingroup numerics
 *
 * The discrete-time algebraic Riccati equation is given by
 * \f[
 *       A^T X A - X - A^T X B ( B^T X B + R )^{-1} B^T X A + Q = 0
 * \f]
 *
 * \f$ A \f$ and \f$ B \f$ must be stabilizable.
 *
 * The resulting gain matrix is given by
 * \f[ G = (B^T X B R)^{-1} B^T X A \f]
 *
 * @see AlgebraicRiccatiContinuous LyapunovDiscrete LyapunovContinuous
 *      SylvesterContinuous SylvesterDiscrete
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Add a method that determiens if sufficient conditions hold and that (A,B) is stabilizable.
 * @todo Add support for the more generalized Riccati equation (augmented by matrix S and E)
 *       E.g. Hamiltonian for augmentation with S: http://web.dcsc.tudelft.nl/~cscherer/books/thesis.pdf
 */
class AlgebraicRiccatiDiscrete
{
 public:
    /**
     * @brief Solve discrete-time algebraic Riccati equation
     *
     * Solve \f$ A^T X A - X - A^T X B ( B^T X B + R )^{-1} B^T X A + Q = 0 \f$
     * w.r.t. \f$ X \f$.
     *
     * \f$ A \f$ and \f$ B \f$ must be stabilizable (all eigenvalues of A outside the unit circle must be controllable).
     *  In addition, the associated symplectic pencil must have no eigenvalue on the unit circle.
     *  Sufficient conditions: R > 0.
     *
     * Also, \f$ A \f$ must be  \b invertible for this implementation.
     *
     * This method optionally returns the gain matrix \f$ G = (B^T X B + R)^{-1} B^T X A \f$.
     *
     * @warning Input matrix sizes are checked only in debug mode and occuring issues immediately break execution
     *
     * @param[in]  A    [n x n] matrix
     * @param[in]  B    [n x m] matrix
     * @param[in]  Q    [n x n] matrix
     * @param[in]  R    [m x m] matrix
     * @param[out] X    Solution of the algebraic Riccati equation with size [n x n] (must not be preallocated).
     * @param[out] G    Gain matrix with size [m x n] (optional)
     * @returns true if a finite solution is found, false otherwise (but no NaN checking is performend on X, use X.allFinite()).
     */
    static bool solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                      const Eigen::Ref<const Eigen::MatrixXd>& Q, const Eigen::Ref<const Eigen::MatrixXd>& R, Eigen::MatrixXd& X,
                      Eigen::MatrixXd* G = nullptr);

    /**
     * @brief Determine if the closed-loop system is stable
     *
     * The closed-loop system is stable if the closed-loop eigenvalues
     * are within the unit circle.
     * This method returns \c true if \f$ A - B G \f$ has only eigenvalues within the unit circle.
     * Hereby, G denotes the gain matrix obtained from solve().
     *
     * @param[in]  A    Open-loop system matrix [n x n]
     * @param[in]  B    Input matrix [n x m]
     * @param[in]  G    Feedback gain matrix according to the solution of the algebraic Riccati equation [m x n]
     * @return true if the closed-loop system is stable, false otherwise.
     */
    static bool isClosedLoopStable(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                   const Eigen::Ref<const Eigen::MatrixXd>& G);

 protected:
    /**
     * @brief Solve Hamiltonian via Schur method
     *
     * The Hamiltonian constructed in solve() is
     * solved via the Schur decomposition.
     * The Real Schur form of A is reordered such that blocks
     * on the diagonal with eigenvalues contained in the unit circle.
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

    // Predicate for Schur block ordering (see solveRiccatiHamiltonianSchur())
    static bool isInsideUnitCircle(const Eigen::Ref<const Eigen::MatrixXd>& block) { return std::abs(block.eigenvalues()[0]) < 1; }
};

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_ALGEBRAIC_RICCATI_DISCRETE_H_

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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_SCHUR_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_SCHUR_H_

#include <Eigen/Core>

namespace corbo {

/**
 * @brief Perform the 2D Real Schur decomposition
 * @ingroup numerics
 *
 * In contrast to Eigen::RealSchur this function enforces diagonal blocks
 * of complex eigenvalues to be in the form:
 * \f[T =
 *   \begin{bmatrix}
 *     a & b \\
 *     c & d
 *   \end{bmatrix}
 * \f]
 * with a denoting the real part and \f$ \sqrt{b c} \f$ the imaginary part.
 * It is \f$ a c < 0 \f$.
 *
 * This function is mainly used in swap_blocks_schur_form() to ensure that swapped blocks
 * are in standard form.
 *
 * The implementation is based on LAPACKs SLANV2 routine [1].
 *
 * [1] http://www.netlib.org/lapack/explore-3.1.1-html/slanv2.f.html#SLANV2.1
 *
 * @param[in, out]  T   [2 x 2] matrix which is subject to decomposition, this matrix will be overwritten
 *                      by the resulting upper triangular matrix.
 * @param[out]      U   Unitary matrix which transforms T to Schur form: T = U^T T U
 **/
void schur_decomposition_2d(Eigen::Ref<Eigen::Matrix2d> T, Eigen::Ref<Eigen::Matrix2d> U);

/**
 * @brief Swap two consecutive diagonal blocks in Real Schur form
 * @ingroup numerics
 *
 * Let \f$ T \f$ denote an upper [n x n] triangular matrix in Real Schur form
 * with >= 2 blocks on the diagonal.
 * Each block is either of size 1 x 1 (real eigenvalue)
 * or 2 x 2 (complex eigenvalue).
 * 2 x 2 blocks are assumed to be in standard form:
 * \f[
 *   \begin{bmatrix}
 *     \alpha & \beta \\
 *     \gamma & \alpha
 *   \end{bmatrix}
 * \f]
 * with \f$ \beta \gamma < 0 \f$. Hereby, the eigenvalue is \f$ \alpha \pm \sqrt{\alpha \gamma} \f$
 *
 * This method swaps two diagonal blocks \f$ A_{11} \f$ and \f$ A_{22} \f$ in \f$ T \f$
 *  without changing the actual eigenvalues of T.
 * The submatrix \f$ A \f$ is denoted as follows:
 * \f[
 *   \begin{bmatrix}
 *     A_{11} & A_{12} \\
 *     0      & A_{22}
 *   \end{bmatrix}
 * \f]
 * Hereby, \f$ A_{11} }f$ denotes the [p x p] matrix corresponding to the diagonal block of the first eigenvalue (p=\{1,2\}).
 * \f$ A_{22} \f$ is the [q x q] block matrix for the second eigenvalue (q=\{1,2\}) and \f$ A_{12} \f$ the [p x q] matrix
 * corresponding to the upper offdiagonal block.
 *
 * The swapped matrix is computed according to \f$ T_s = Q^T * T * Q \f$.
 *
 * The implementation is based on [1] with some details of [2].
 *
 * [1] Z. Bai and J. W. Demmel. "On swapping diagonal blocks in real Schur form."
 *     Linear Algebra Appl., 186:73–95, 1993.
 * [2] LAPACK routine SLAEXC: http://www.netlib.org/lapack/explore-3.1.1-html/slaexc.f.html
 *
 * @param[in, out]  T            [n x n] matrix in Real Schur form which is replaced by the one with swapped eigenvalues.
 * @param[in]       ra11         row idx / column idx of the upper left diagonal block (A11) subject to swapping.
 * @param[in]       p            Dimension of the first diagonal block to be swapped (A11).
 * @param[in]       q            Dimension of the second diagonal block to be swapped (A22).
 * @param[in, out]  Q            [n x n] aggregated orthogonal matrix which performs the swapping, in particular Q^T T Q
 *                               (warning: Q must be preallocated, at least as identity matrix!)
 * @param[in]       standardize   if true, swapped 2 x 2 complex blogs are standardized according to the form mentioned above.
 * @returns true if swapping was complete, false otherwise
 **/
bool swap_schur_blocks(Eigen::Ref<Eigen::MatrixXd> T, int ra11, int p, int q, Eigen::Ref<Eigen::MatrixXd> Q, bool standardize = false);

/**
 * @brief Reorder blocks in Real Schur
 * @ingroup numerics
 *
 * This method reorders [1 x 1] and [2 x 2] diagonal blocks of in a Real Schur Matrix
 * according to a given predicate function.
 * Blocks for which the predicate returns \c true are sorted in the upper left part of the diagonal.
 * The columns of the corresponding uniaty matrix span the related invariant subspace which
 * is for example useful for solving the algebraic riccati equation.
 * The signature of the unary predicate function is as follows:
 *
 *    bool predicate(const Eigen::Ref<const Eigen::MatrixXd>& block);
 *
 * The input argument is the current [1 x 1] or [2 x 2] block.
 *
 * Example for sorting a RealSchur Matrix T and its corresponding unitary transformation matrix
 * such that eigenvalues with negative real part are in the upper left part:
 *
 * \code{.cpp}
 *    Eigen::MatrixXd A = Eigen::MatrixXd::Random(6,6); // create random 6 x 6 matrix
 *    Eigen::RealSchur<Eigen::MatrixXd> schur(A); // perform schur decomposition, include <Eigen/Eigenvalues>
 *    Eigen::MatrixXd T = schur.matrixT(); // tridiagonal matrix with 1 x 1 or 2 x 2 blocks on the diagonal containing eigenvalues of T
 *    Eigen::MatrixXd Q = schur.matrixU(); // unitary transformation matrix, such that T = Q* A Q.
 *    // now create predicate function, e.g. here with a lambda function (block is either 1 x 1 or 2 x 2:
 *    auto select = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; };
 *    int subspace_dim; // receive dimension of the invariant subspace if successful
 *    bool success = reorder_schur_blocks(T, Q, select, &subspace_dim);
 * \endcode
 *
 * The implementation is based on Algorithm 1 in [1].
 *
 * [1] D. Kressner. "Block Algorithms for Reordering Standard and Generalized Schur Forms", LAPACK Working Note 171,
 *     http://www.netlib.org/lapack/lawnspdf/lawn171.pdf
 *
 * @tparam          Predicate     Templated type to represent a unary predicate function as specified above;
 *                                Usually, the type is deduced automatically from the function argument list.
 *
 * @param[in, out]  T             [n x n] matrix which is subject to reordering, this matrix will be overwritten
 *                                by the resulting upper triangular matrix.
 * @param[in, out]  Q             Unitary matrix which transforms A to Schur form: T = Q T Q^T
 * @param[in]       predicate     Unary predicate function that determines if the current block should be moved to
 *                                the upper left part of T (if predicate returns \c true). The predicate function
 *                                takes an <tt> const Eigen::Ref<const Eigen::MatrixXd>& </tt> as input argument.
 * @param[out]      subspace_dim  If not \c nullptr: contains the size of the invariant subspace according to the
 *                                selected eigenvalues [optional]
 * @param[in]       standardize   If \c true the underlying swapping method is enforced to standardize swapped
 *                                2 x 2 blocks according to the description in swap_schur_blocks().
 **/
template <class Predicate>
bool reorder_schur_blocks(Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> Q, Predicate predicate, int* subspace_dim = nullptr,
                          bool standardize = false);

}  // namespace corbo

#include <corbo-numerics/schur.hpp>

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_SCHUR_H_

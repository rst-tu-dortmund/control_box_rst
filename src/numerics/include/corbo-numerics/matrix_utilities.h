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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_MATRIX_UTILITIES_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_MATRIX_UTILITIES_H_

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <type_traits>

namespace corbo {

/**
 * @brief Determine if a given matrix is square
 * @ingroup numerics
 *
 * @param matrix generic matrix type to be tested
 * @returns true if square, false otherwise.
 **/
template <typename Derived>
inline bool is_square(const Eigen::MatrixBase<Derived>& matrix)
{
    return matrix.rows() == matrix.cols();
}

/**
 * @brief Determine if a given number of elements defines a square matrix
 * @ingroup numerics
 *
 * @param numel number of elements to be checked for square matrix
 * @returns true if square, false otherwise.
**/
inline bool is_square(int numel)
{
    double q = std::sqrt(numel);
    if (q * q != numel) return false;
    return true;
}

/**
 * @brief Determine if two matrices exhibit equal sizes/dimensions
 * @ingroup numerics
 *
 * @param matrix1 generic matrix type
 * @param matrix2 generic matrix type
 * @returns true if equal dimensions, false otherwise.
 **/
template <typename DerivedA, typename DerivedB>
inline bool have_equal_size(const Eigen::MatrixBase<DerivedA>& matrix1, const Eigen::MatrixBase<DerivedB>& matrix2)
{
    return matrix1.rows() == matrix2.rows() && matrix1.cols() == matrix2.cols();
}

/**
 * @brief Determine if N matrices exhibit equal sizes/dimensions
 * @ingroup numerics
 *
 * @param matrix1 generic matrix type
 * @param matrix2 generic matrix type
 * @returns true if equal dimensions, false otherwise.
 **/
template <typename DerivedA, typename DerivedB, typename... Derived>
inline bool have_equal_size(const Eigen::MatrixBase<DerivedA>& matrix1, const Eigen::MatrixBase<DerivedB>& matrix2,
                            const Eigen::MatrixBase<Derived>&... matrices_other)
{
    return have_equal_size(matrix1, matrix2) ? have_equal_size(matrix2, matrices_other...) : false;
}

/**
 * @brief Determine if a given matrix is positive definite
 * @ingroup numerics
 *
 * The current implementation performs Eigen's Cholesky
 * decomposition in order to test whether the given
 * matrix appears to be positive definite.
 *
 * @param matrix generic matrix type to be tested
 * @returns true if positive definite, false otherwise.
 **/
template <typename Derived>
inline bool is_positive_definite(const Eigen::MatrixBase<Derived>& matrix)
{
    if (!is_square(matrix)) return false;
    Eigen::LLT<Eigen::MatrixXd> cholesky(matrix);
    return cholesky.info() != Eigen::NumericalIssue;
}

/**
 * @brief Compute condition number of a square matrix using the Eigenvalues
 * @ingroup numerics
 *
 * This method computes the eigenvalues and returns abs(eig_max / eig_min).
 * For a version based on the SVD refer to compute_condition_number()
 *
 * @param matrix generic matrix type
 * @returns condition number
 **/
template <typename Derived>
inline double compute_condition_number_square_matrix(const Eigen::MatrixBase<Derived>& matrix)
{
    if (!is_square(matrix))
    {
        return -1.0;
    }
    auto& eigenvalues = matrix.eigenvalues();
    auto eig_min      = eigenvalues.minCoeff();
    auto eig_max      = eigenvalues.maxCoeff();
    // if (std::abs(eig_min) < 1e-15)
    //    return eig_max > 0 ? CORBO_INF_DBL : -CORBO_INF_DBL;

    return std::abs(eig_max / eig_min);
}

/**
 * @brief Compute condition number of matrix using SVD
 * @ingroup numerics
 *
 * @param matrix generic matrix type
 * @returns condition number
 **/
template <typename Derived>
inline double compute_condition_number(const Eigen::MatrixBase<Derived>& matrix)
{
    if (matrix.rows() == 0 || matrix.cols() == 0) return -1.0;

    const typename Eigen::JacobiSVD<typename Derived::PlainObject>::SingularValuesType& sing_vals = matrix.jacobiSvd().singularValues();
    // Singular values are always sorted in decreasing order.
    double sigma_max = sing_vals[0];
    double sigma_min = sing_vals[sing_vals.size() - 1];
    // if (sigma_min < 1e-15)
    //    return sigma_max > 0 ? CORBO_INF_DBL : -CORBO_INF_DBL;

    return sigma_max / sigma_min;
}

/**
 * @brief Compute the pseudo inverse using SVD
 * @ingroup numerics
 *
 **/
template <typename Derived>
void compute_pseudo_inverse(const Eigen::MatrixBase<Derived>& matrix, Eigen::MatrixBase<Derived>& pinvmat, double tolerance = 1e-6)
{
    Eigen::JacobiSVD<typename Derived::PlainObject> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // TODO(roesmann): do we need FullU and FullV, or is Thin sufficient?

    using SingularValuesType    = typename Eigen::JacobiSVD<typename Derived::PlainObject>::SingularValuesType;
    SingularValuesType sing_inv = svd.singularValues();

    for (int i = 0; i < sing_inv.size(); ++i)
    {
        if (sing_inv[i] > tolerance)
            sing_inv[i] = 1.0 / sing_inv[i];
        else
            sing_inv[i] = 0.0;
    }
    pinvmat.noalias() = svd.matrixV() * sing_inv.asDiagonal() * svd.matrixU().transpose();
}

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_MATRIX_UTILITIES_H_

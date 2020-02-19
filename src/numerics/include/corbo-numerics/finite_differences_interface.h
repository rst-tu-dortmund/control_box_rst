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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_FINITE_DIFFERENCES_INTERFACE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_FINITE_DIFFERENCES_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/numerics/finite_differences.pb.h>
#endif

#include <functional>
#include <memory>

namespace corbo {

/**
 * @brief Interface class for finite difference approaches
 *
 * @ingroup numerics
 *
 * This base class provides an interface for computing
 * finite differences as derivative approximation
 * for some relevant functions.
 * For instance for computing the Jacobian matrices
 * of the system dynamics function $\f f(x, u) \$f
 * w.r.t. \f$ x \f$ and/or \f$ u \f$.
 *
 * @remark This interface is provided with factory support (FiniteDifferencesFactory).
 *
 * @see ForwardDifferences CentralDifferences NumericalIntegratorExplicitInterface
 *      NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class FiniteDifferencesInterface
{
 public:
    using Ptr  = std::shared_ptr<FiniteDifferencesInterface>;
    using UPtr = std::unique_ptr<FiniteDifferencesInterface>;

    using StateVector = Eigen::VectorXd;
    using InputVector = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~FiniteDifferencesInterface() {}
    //! Return a newly allocated instances of the inherited class.
    virtual FiniteDifferencesInterface::Ptr getInstance() const = 0;

    /**
     * @brief Compute Jacobian of a desired function
     *
     * Given a function \f$ f(x) \f$ with \f$ f: \mathbb{R}^p \to \mathbb{R}^n \f$.
     * This method computes the Jacobian matrix
     * \f[
     *     Jf_x(x_0) =
     *     \begin{bmatrix}
     *        \partial_{x_1} f_1(x_0) & \partial_{x_2} f_1(x_0) & \cdots & \partial_{x_p} f_1(x_0) \\
     *        \partial_{x_1} f_2(x_0) & \partial_{x_2} f_2(x_0) & \cdots & \partial_{x_p} f_2(x_0) \\
     *        \vdots \\
     *        \partial_{x_1} f_p(x_0) & \partial_{x_2} f_p(x_0) & \cdots & \partial_{x_p} f_p(x_0)
     *     \end{bmatrix}
     * \f]
     * The Jacobian is evaluated at \f$ x_0 \f$.
     *
     * The function \f$ f \f$ is defined by a function of the following prototype:
     * \code
     *     void fun(Eigen::Ref<Eigen::VectorXd> values);
     * \endcode
     * Hereby, \c fun should return the function values via output argument \c values as (n x 1) vector
     * by implicitly taking the current working point \f$ x_0 \f$ into account. Note, you might want to bind other parameters to
     * \c fun by using lambdas or std::bind().
     *
     * Furthermore, a second function needs to be provided, in particular an implempation of the increment operator.
     * Within Jacobian computation, the variable $\f x $\f needs to be varied (for computing finite differences).
     * This can be achieved in a callback fashion by providing the increment function according to the following prototype:
     * \code
     *     void inc(int idx, const double& value);
     * \endcode
     * Herebey, the first parameter determines which index (element in \f$ x \f$) needs to be incremented,
     * and the second contains the actual increment value.
     *
     * Example implementation:
     * \code
     *   Eigen::Vector2d x(1,1); // variable x initialized with [1, 1]^T as current working point (dimension p=2)
     *   auto inc = [&x](int idx, const double& value){x[idx] += value;}; // define increment as simple x[idx]+=value, also bind x.
     *   auto fun = [&x](Eigen::Ref<Eigen::VectorXd> values)
     *   {
     *      values[0] = x[0] * x[0] + 1;  // f1 = x1^2 + 1
     *      values[1] = x[1] * x[0] + 2;  // f2 = x1*x2 + 2
     *   }; // dimension n=2
     *   Eigen::MatrixXd jacobian(2,2); // we must initialize/pre-allocate the jacobian with the correct size [p x n]
     *   computeJacobian(inc, fun, jacobian);
     * \endcode
     *
     * In case the function callback is defined in terms of taking an \c Eigen::VectorXd& into account
     * rather than \c Eigen::Ref<Eigen::VectorXd> refer to overload computeJacobian2()
     *
     * @see computeJacobian2 computeJacobianAndHessian
     *
     * @param[in]   inc_fun   Function callback to the increment operator function
     * @param[in]   eval_fun  Function callback to the function evaluation
     * @param[out]  jacobian  The resulting Jacobian matrix (warning: \c jacobian must be preallocated as n x p matrix)
     */
    virtual void computeJacobian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                 Eigen::Ref<Eigen::MatrixXd> jacobian) = 0;

    /**
     * @brief Compute Jacobian of a desired function (overload which accepts a slightly different callback function)
     *
     * Refer to the documentation of computeJacobian().
     * This method takes type \c Eigen::VectorXd& into account for the function evaluation callback
     * rather than \c Eigen::Ref<Eigen::VectorXd>.
     *
     * @see computeJacobian computeJacobianAndHessian2
     *
     * @param[in]   inc_fun   Function callback to the increment operator function
     * @param[in]   eval_fun  Function callback to the function evaluation
     * @param[out]  jacobian  The resulting Jacobian matrix (warning: \c jacobian must be preallocated as n x p matrix)
     */
    virtual void computeJacobian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                  Eigen::Ref<Eigen::MatrixXd> jacobian) = 0;

    /**
     * @brief Compute Hessian of a desired function
     *
     * Given a function \f$ f(x) \f$ with \f$ f: \mathbb{R}^p \to \mathbb{R}^n \f$.
     * This method computes the \f$ n \f$ Hessian matrices of \f$ f(x) \f$ (element-wise) and accumulates them (optionally with individual scale
     * factors).
     * In case \f$ n=1 \f$, the Hessian is as follows:
     * \f[
     *     Hf_{1x}(x_0) =
     *     \begin{bmatrix}
     *        \partial^2_{x_1 x_1} f_1(x_0) & \partial^2_{x_1 x_2} f_1(x_0) & \cdots & \partial^2_{x_1 x_p} f_1(x_0) \\
     *        \partial^2_{x_2 x_1} f_1(x_0) & \partial^2_{x_2 x_2} f_1(x_0) & \cdots & \partial^2_{x_2 x_p} f_1(x_0) \\
     *        \vdots \\
     *        \partial^2_{x_p x_1} f_1(x_0) & \partial^2_{x_p x_2} f_1(x_0) & \cdots & \partial^2_{x_p x_p} f_1(x_0)
     *     \end{bmatrix}
     * \f]
     * The Hessian is evaluated at \f$ x_0 \f$.
     *
     * Function prototypes \c inc_fun and \c eval_fun are similar to computeJacobian().
     * In case the function callback is defined in terms of taking an \c Eigen::VectorXd& into account
     * rather than \c Eigen::Ref<Eigen::VectorXd> refer to overload computeHessian2()
     *
     * @see computeHessian2 computeJacobianAndHessian
     *
     * @param[in]   inc_fun      Function callback to the increment operator function
     * @param[in]   eval_fun     Function callback to the function evaluation
     * @param[in]   dim_f        Dimension of the function value vector (obtained from eval_fun)
     * @param[out]  hessian      The resulting Hessian matrix (warning: \c hessian must be preallocated as p x p matrix)
     * @param[in]   multipliers  (optional) Vector of multipliers for scaling individual Hessian terms [dim_f x 1].
     */
    virtual void computeHessian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun, int dim_f,
                                Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers = nullptr) = 0;

    /**
     * @brief Compute Hessian of a desired function (overload which accepts a slightly different callback function)
     *
     * Refer to the documentation of computeHessian().
     * This method takes type \c Eigen::VectorXd& into account for the function evaluation callback
     * rather than \c Eigen::Ref<Eigen::VectorXd>.
     *
     * @see computeHessian computeJacobianAndHessian2
     *
     * @param[in]   inc_fun      Function callback to the increment operator function
     * @param[in]   eval_fun     Function callback to the function evaluation
     * @param[in]   dim_f        Dimension of the function value vector (obtained from eval_fun)
     * @param[out]  hessian      The resulting Hessian matrix (warning: \c hessian must be preallocated as p x p matrix)
     * @param[in]   multipliers  (optional) Vector of multipliers for scaling individual Hessian terms [dim_f x 1].
     */
    virtual void computeHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun, int dim_f,
                                 Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers = nullptr) = 0;

    /**
     * @brief Compute Jacobian and Hessian of a desired function
     *
     * This method is often faster than invoking computeJacobian() and computeHessian() separately.
     *
     * Refer to the individual methods to obtain details regarding arguments and usage.
     *
     * EvalFun should be of type void(Eigen::VectorXd& values) or similar and should return current function values.
     * IncFun should be of type void(int idx, double value) and should increment component \c idx by \c value.
     * The dimension of the function provided by \c eval_f is obtained by checking the number of allocated rows
     * of the \c jacobian variable.
     * In case dim_f > 1, the hessian is computed for each dimension and summed up.
     * Optionally, individual hessians might be scaled by a multiplier which is common in solving nonlinear programs
     * (e.g., Lagrangian). If provided, the multiplier vector should be of size [dim_f x 1].
     *
     * Function prototypes \c inc_fun and \c eval_fun are similar to computeJacobian().
     * In case the function callback is defined in terms of taking an \c Eigen::VectorXd& into account
     * rather than \c Eigen::Ref<Eigen::VectorXd> refer to overload computeJacobianAndHessian2()
     *
     * @note The Jacobian is also scaled by the multipliers if provided.
     *
     * @see computeJacobian computeHessian computeJacobianAndHessian2
     *
     * @param[in]   inc_fun      Function callback to the increment operator function
     * @param[in]   eval_fun     Function callback to the function evaluation
     * @param[out]  jacobian     The resulting Jacobian matrix (warning: \c jacobian must be preallocated as n x p matrix)
     * @param[out]  hessian      The resulting Hessian matrix (warning: \c hessian must be preallocated as p x p matrix)
     * @param[in]   multipliers  (optional) Vector of multipliers for scaling individual Hessian terms [dim_f x 1].
     */
    virtual void computeJacobianAndHessian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                           Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                           const double* multipliers = nullptr) = 0;

    /**
     * @brief Compute Jacobian and Hessian of a desired function (overload which accepts a slightly different callback function)
     *
     * Refer to the documentation of computeJacobianAndHessian().
     * This method takes type \c Eigen::VectorXd& into account for the function evaluation callback
     * rather than \c Eigen::Ref<Eigen::VectorXd>.
     *
     * @note The Jacobian is also scaled by the multipliers if provided.
     *
     * @see computeJacobianAndHessian computeJacobian computeHessian
     *
     * @param[in]   inc_fun      Function callback to the increment operator function
     * @param[in]   eval_fun     Function callback to the function evaluation
     * @param[out]  jacobian     The resulting Jacobian matrix (warning: \c jacobian must be preallocated as n x p matrix)
     * @param[out]  hessian      The resulting Hessian matrix (warning: \c hessian must be preallocated as p x p matrix)
     * @param[in]   multipliers  (optional) Vector of multipliers for scaling individual Hessian terms [dim_f x 1].
     */
    virtual void computeJacobianAndHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                            Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                            const double* multipliers = nullptr) = 0;

#ifdef MESSAGE_SUPPORT
    //! Export selected class and parameters to a given message
    virtual void toMessage(messages::FiniteDifferences& message) const {}
    //! Import selected class and parameters from a given message (optionally pass any issues to the caller).
    virtual void fromMessage(const messages::FiniteDifferences& message, std::stringstream* issues = nullptr) {}
#endif
};

using FiniteDifferencesFactory = Factory<FiniteDifferencesInterface>;
#define FACTORY_REGISTER_FINITE_DIFFERENCES(type) FACTORY_REGISTER_OBJECT(type, FiniteDifferencesInterface)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_FINITE_DIFFERENCES_INTERFACE_H_

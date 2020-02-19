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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_LINEAR_BENCHMARK_SYSTEMS_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_LINEAR_BENCHMARK_SYSTEMS_H_

#include <corbo-systems/system_dynamics_interface.h>

#include <corbo-numerics/matrix_utilities.h>

namespace corbo {

/**
 * @brief System dynamics for a series of integrators (continuous-time)
 *
 * @ingroup systems
 *
 * An continuous-time integrator chain of order \f$ p \f$ is defiend as:
 * \f[
 *      T x^{(p)} = u
 * \f]
 * in which superscript \f$ (p) \f$ defines the p-th derivative of x w.r.t. time.
 * \f$ T \f$ denotes a time constant.
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SerialIntegratorSystem : public SystemDynamicsInterface
{
 public:
    //! Default constructor (do not forget to set the dimension)
    SerialIntegratorSystem() {}
    //! Construct ingerator system with given order/dimension
    explicit SerialIntegratorSystem(int dimension) : _dimension(dimension) {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<SerialIntegratorSystem>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return true; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return _dimension; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.rows() == _dimension);
        assert(x.rows() == f.rows() && "SerialIntegratorSystem::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        const int num_integrator_states = _dimension - 1;

        if (_dimension > 1) f.head(num_integrator_states) = x.segment(1, num_integrator_states);  // x^(1:n-1) = x^(2:n)
        // add last equation containing control u
        f[num_integrator_states] = u[0] / _time_constant;  // x^(n) = u / T
    }

    // access parameters

    //! Get current integrator chain dimension / order of the system
    const int& getDimension() const { return _dimension; }
    //! Set integrator dimension (p >= 1)
    void setDimension(int dim) { _dimension = dim; }
    //! Get time constant T of the integrator
    const double& getTimeConstant() const { return _time_constant; }
    //! Set Time constant T of the integrator
    void setTimeConstant(double time_constant) { _time_constant = time_constant; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_serial_integrators()->set_dimension(_dimension);
        message.mutable_serial_integrators()->set_time_constant(_time_constant);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        _dimension     = message.serial_integrators().dimension();
        _time_constant = message.serial_integrators().time_constant();
    }
#endif

 private:
    int _dimension        = 1;
    double _time_constant = 1.0;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(SerialIntegratorSystem)

class ParallelIntegratorSystem : public SystemDynamicsInterface
{
 public:
    //! Default constructor (do not forget to set the dimension)
    ParallelIntegratorSystem() {}
    //! Construct ingerator system with given order/dimension
    explicit ParallelIntegratorSystem(int dimension) : _dimension(dimension) {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<ParallelIntegratorSystem>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return true; }

    // implements interface method
    int getInputDimension() const override { return _dimension; }
    // implements interface method
    int getStateDimension() const override { return _dimension; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.rows() == _dimension);
        assert(x.rows() == f.rows() && "ParallelIntegratorSystem::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f = _time_constant * u;
    }

    // access parameters

    //! Get current integrator chain dimension / order of the system
    const int& getDimension() const { return _dimension; }
    //! Set integrator dimension (p >= 1)
    void setDimension(int dim) { _dimension = dim; }
    //! Get time constant T of the integrator
    const double& getTimeConstant() const { return _time_constant; }
    //! Set Time constant T of the integrator
    void setTimeConstant(double time_constant) { _time_constant = time_constant; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_parallel_integrators()->set_dimension(_dimension);
        message.mutable_parallel_integrators()->set_time_constant(_time_constant);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        _dimension     = message.parallel_integrators().dimension();
        _time_constant = message.parallel_integrators().time_constant();
    }
#endif

 private:
    int _dimension        = 1;
    double _time_constant = 1.0;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(ParallelIntegratorSystem)

class LinearStateSpaceModel : public SystemDynamicsInterface
{
 public:
    //! Default constructor (do not forget to set the dimension)
    LinearStateSpaceModel() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<LinearStateSpaceModel>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return true; }

    // implements interface method
    int getInputDimension() const override { return _B.cols(); }
    // implements interface method
    int getStateDimension() const override { return _A.rows(); }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "LinearStateSpaceModel::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f = _A * x + _B * u;
    }

    // access parameters

    //! Set Time constant T of the integrator
    void setParameters(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B)
    {
        _A = A;
        _B = B;
    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        messages::LinearStateSpaceModel* msg = message.mutable_linear_state_space();

        // A
        msg->mutable_a()->Resize(_A.size(), 0);
        Eigen::Map<Eigen::VectorXd>(msg->mutable_a()->mutable_data(), _A.size()) = _A;

        // B
        msg->mutable_b()->Resize(_B.size(), 0);
        Eigen::Map<Eigen::VectorXd>(msg->mutable_b()->mutable_data(), _B.size()) = _B;
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        const messages::LinearStateSpaceModel& msg = message.linear_state_space();

        if (msg.a_size() == 0 || msg.b_size() == 0)
        {
            *issues << "LinearStateSpaceModel: Please specify a proper system matrix A and input matrix B.\n";
            return;
        }

        // A
        int numel_a = msg.a_size();
        if (!is_square(numel_a))
        {
            *issues << "LinearStateSpaceModel: system matrix A is not square.\n";
            return;
        }
        int p = std::sqrt(numel_a);
        _A    = Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.a().data(), p, p);

        // B
        int numel_b = msg.b_size();

        double q = numel_b / p;
        if (p * q != numel_b)
        {
            *issues << "LinearStateSpaceModel: input matrix B is not a valid matrix (number of columns must be " << p << ").\n";
            return;
        }
        _B = Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.b().data(), p, (int)q);
    }
#endif

 private:
    Eigen::MatrixXd _A = Eigen::MatrixXd::Zero(1, 1);
    Eigen::MatrixXd _B = Eigen::MatrixXd::Zero(1, 1);
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(LinearStateSpaceModel)

/**
 * @brief System dynamics for a second order integrator (discrete-time)
 *
 * @ingroup systems
 *
 * An 2nd-order discrete-time integrator chain is defiend as:
 * \f[
 *      x_{k+1} = \begin{bmatrix} 1 & dt   \\ 0 & 1 \end{bmatrix} x_k
 *              + \begin{bmatrix} 0.5 dt^2 \\ 0     \end{bmatrix} u_k
 * \f]
 * in which \f$ dt \f$ denotes the sampling time.
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class DoubleIntegratorDiscreteTime : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    DoubleIntegratorDiscreteTime() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<DoubleIntegratorDiscreteTime>(); }

    // implements interface method
    bool isContinuousTime() const override { return false; }
    // implements interface method
    bool isLinear() const override { return true; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.rows() == 2);
        assert(x.rows() == f.rows() && "IntegratorSystemCont::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[0] + _dt * x[1] + 0.5 * _dt * _dt * u[0];
        f[1] = x[1] + _dt * u[0];
    }

    //! Get current sampling time dt
    const double& getDt() const { return _dt; }
    //! Set sampling time dt
    void setDt(double dt) { _dt = dt; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);
        message.mutable_double_integrator_discrete_time()->set_dt(_dt);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issue) override
    {
        SystemDynamicsInterface::fromMessage(message, issue);
        _dt = message.double_integrator_discrete_time().dt();
    }
#endif

 private:
    double _dt = 1.0;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(DoubleIntegratorDiscreteTime)

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_LINEAR_BENCHMARK_SYSTEMS_H_

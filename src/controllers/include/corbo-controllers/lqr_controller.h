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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_LQR_CONTROLLER_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_LQR_CONTROLLER_H_

#include <corbo-controllers/controller_interface.h>
#include <corbo-systems/system_dynamics_interface.h>
#include <memory>

namespace corbo {

/**
 * @brief LQR controller
 *
 * @ingroup controllers
 *
 * Implementation of a LQR controller.
 *
 * A continuous-time LQR requires linear system dynamics
 * \f[
 *      \dot{x} = A x + B u
 * \f]
 * and minimizes the following infinite horizon cost function
 * \f[
 *      J = \int_0^\infty \big( \tilde{x}^T Q \tilde{x} + \tilde{u}^T R \tilde{u} \big) dt
 * \f]
 * with $\f \tilde{x} = x - x_{ref} \f$ and \f$ \tilde{u} = u - u_{ref} \f$
 *
 * The solution leads to a gain matrix \f$ K \f$ and the control law is \f$ u = - K x \f$.
 *
 * For discrete-time system dynamics f\$ x_{k+1} = A x_k + B u_k \f$
 * the cost function is defined as:
 * \f[
 *      J = \sum_k^\infty \big(  \tilde{x}_k^T Q \tilde{x}_k + \tilde{u}_k^T R \tilde{u}_k \big)
 * \f]
 * with $\f \tilde{x}_k = x_k - x_{ref} \f$ and \f$ \tilde{u}_k = u_k - u_{ref} \f$
 *
 * Both continuous-time and discrete-time systems are supported.
 * Depending on the system type either the continuous-time or discrete-time
 * Riccati equation is solved.
 *
 * Note, for nonlinear system dynamics, the dynamics equation provided via
 * a SystemDynamicsInterface are linearized at \f$ x_{ref} \f$ and \f$ u_{ref} \f$.
 *
 * For considering constraints on controls and/or states as well as for general nonlinear system
 * dynamics you might refer to the PredictiveController.
 *
 * @see ControllerInterface PidController PredictiveController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class LqrController : public ControllerInterface
{
 public:
    using Ptr = std::shared_ptr<LqrController>;

    LqrController();
    explicit LqrController(SystemDynamicsInterface::Ptr system_model);

    // implements interface method
    int getControlInputDimension() const override { return _system_model ? _system_model->getInputDimension() : 0; }
    // implements interface method
    int getStateDimension() const override { return _system_model ? _system_model->getStateDimension() : 0; }
    // implements interface method
    bool hasPiecewiseConstantControls() const override { return false; }
    // implements interface method
    bool providesFutureControls() const override { return false; }
    // implements interface method
    bool providesFutureStates() const override { return false; }

    // implements interface method
    ControllerInterface::Ptr getInstance() const override { return std::make_shared<LqrController>(); }

    //! Set system model
    void setSystemModel(SystemDynamicsInterface::Ptr system_model);

    /**
     * @brief Set cost function weights
     * @param[in] R   Positive definite weight matrix for state deviation [getStateDimension() x getStateDimension()]
     * @param[in] Q   Positive definite matrix for control input deviation [getControlInputDimension() x getControlInputDimension()]
     * @return true if the matrices are set successfully, false otherwise (e.g. if not positive definite)
     */
    bool setWeights(const Eigen::Ref<const Eigen::MatrixXd>& R, const Eigen::Ref<const Eigen::MatrixXd>& Q);

    //! Set pos. def. weight matrix for penalizing state deviation [getStateDimension() x getStateDimension()]
    bool setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R);
    //! Set pos. def. weight matrix for penalizing control input deviation [getControlInputDimension() x getControlInputDimension()]
    bool setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q);

    // implements interface method
    bool initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                    const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref = nullptr) override;

    // implements interface method
    bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
              TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
              ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
              ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") override;

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    ControllerStatistics::Ptr getStatistics() const override { return {}; }

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::LqrController& message) const;
    void fromMessage(const corbo::messages::LqrController& message, std::stringstream* issues = nullptr);
    // implements interface method
    void toMessage(corbo::messages::Controller& message) const override { toMessage(*message.mutable_lqr_controller()); }
    // implements interface method
    void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.lqr_controller(), issues);
    }
#endif

    // implements interface method
    void reset() override;

    //! Specify whether the state error should be published via signal
    void setPublishError(bool publish) { _publish_error = publish; }

 private:
    SystemDynamicsInterface::Ptr _system_model;
    bool _discrete_time = false;

    bool _initialized = false;

    bool _publish_error = true;

    double _dt = 0.1;

    Eigen::MatrixXd _A;
    Eigen::MatrixXd _B;

    Eigen::MatrixXd _R;
    Eigen::MatrixXd _Q;

    Eigen::MatrixXd _S;

    Eigen::MatrixXd _K;
};

FACTORY_REGISTER_CONTROLLER(LqrController)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_LQR_CONTROLLER_H_

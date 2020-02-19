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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_SIMPLE_STATE_CONTROLLER_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_SIMPLE_STATE_CONTROLLER_H_

#include <corbo-controllers/controller_interface.h>
#include <memory>

namespace corbo {

/**
 * @brief State feedback controller wigh feedback gain matrix K
 *
 * @ingroup controllers
 *
 * The controller implements the following control law: \f$ u = - K x \f$.
 * Herby the number of rows in K corresponds to state dimension and
 * the number of columns to the control input dimension.
 *
 * @see ControllerInterface LqrController PidController PredictiveController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SimpleStateController : public ControllerInterface
{
 public:
    SimpleStateController() = default;

    // implements interface method
    int getControlInputDimension() const override { return _K.rows(); }
    // implements interface method
    int getStateDimension() const override { return _V.cols() > 0 ? _V.cols() : _K.cols(); }
    // implements interface method
    bool hasPiecewiseConstantControls() const override { return false; }
    // implements interface method
    bool providesFutureControls() const override { return false; }
    // implements interface method
    bool providesFutureStates() const override { return false; }

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<SimpleStateController>(); }

    //! Set feedback gain matrix K [control input dimension x state dimension]
    void setGainMatrixK(const Eigen::Ref<const Eigen::MatrixXd>& K);

    //! Set reference filter matrix V [control input dimension x output dimension]
    void setFilterMatrixV(const Eigen::Ref<const Eigen::MatrixXd>& V);

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

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::Controller& message) const override;
    // implements interface method
    void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) override;
#endif

    // implements interface method
    void reset() override;

    //! Specify whether the state error should be published via signal
    void setPublishError(bool publish) { _publish_error = publish; }

 private:
    bool _publish_error = true;
    Eigen::MatrixXd _K;
    Eigen::MatrixXd _V;
};

FACTORY_REGISTER_CONTROLLER(SimpleStateController)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_SIMPLE_STATE_CONTROLLER_H_

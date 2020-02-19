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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_STEP_RESPONSE_GENERATOR_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_STEP_RESPONSE_GENERATOR_H_

#include <corbo-controllers/controller_interface.h>
#include <memory>

namespace corbo {

/**
 * @brief Step Response Generator
 *
 * @ingroup controllers
 *
 * Generates a constant control input signal.
 *
 * @remarks The state dimension needs to be set using setStateDimension()
 *          due to compatibility purposes even if it is obsolete
 *
 * @see ControllerInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class StepResponseGenerator : public ControllerInterface
{
 public:
    StepResponseGenerator() = default;

    // implements interface method
    int getControlInputDimension() const override { return _u.size(); }
    // implements interface method
    int getStateDimension() const override { return _state_dim; }
    // implements interface method
    bool hasPiecewiseConstantControls() const override { return true; }
    // implements interface method
    bool providesFutureControls() const override { return false; }
    // implements interface method
    bool providesFutureStates() const override { return false; }

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<StepResponseGenerator>(); }
    //! Return a newly created shared instance of the implemented class (static method)
    static Ptr getInstanceStatic() { return std::make_shared<StepResponseGenerator>(); }

    void setControl(const Eigen::VectorXd& u_step) { _u = u_step; }
    const Eigen::VectorXd& getControl() const { return _u; }
    void setStateDimension(int state_dim) { _state_dim = state_dim; }

    // implements interface method
    bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
              TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
              ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
              ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") override
    {
        if (!u_sequence) return false;
        u_sequence->clear();
        u_sequence->add(0.0, _u);
        return true;
    }

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override {}

#ifdef MESSAGE_SUPPORT
    void toMessage(messages::StepResponseGenerator& message) const
    {
        message.set_state_dim(_state_dim);
        message.mutable_u_step()->Resize(_u.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_u_step()->mutable_data(), _u.size()) = _u;
    }
    void fromMessage(const messages::StepResponseGenerator& message, std::stringstream* issues = nullptr)
    {
        _state_dim = message.state_dim();
        if (message.u_step_size() > 0)
            _u = Eigen::Map<const Eigen::VectorXd>(message.u_step().data(), message.u_step_size());
        else
            _u.resize(0);
    }

    // implements interface method
    void toMessage(messages::Controller& message) const override { toMessage(*message.mutable_step_response_generator()); }
    // implements interface method
    void fromMessage(const messages::Controller& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.step_response_generator(), issues);
    }
#endif

    // implements interface method
    void reset() override {}

 private:
    Eigen::VectorXd _u;
    int _state_dim = 0;
};

FACTORY_REGISTER_CONTROLLER(StepResponseGenerator)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_STEP_RESPONSE_GENERATOR_H_

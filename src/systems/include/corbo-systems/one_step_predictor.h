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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_ONE_STEP_PREDICTOR_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_ONE_STEP_PREDICTOR_H_

#include <corbo-numerics/integrator_interface.h>
#include <corbo-systems/system_dynamics_interface.h>
#include <corbo-systems/time_value_buffer.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/systems/one_step_predictor.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief OneStepPredictor
 *
 * @ingroup systems
 *
 * Predict plant output for a single step, e.g. useful for CPU compenation in MPC control.
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class OneStepPredictor
{
 public:
    using StateVector   = Eigen::VectorXd;
    using ControlVector = Eigen::VectorXd;

    //! Default constructor
    OneStepPredictor() = default;

    double getDeadTime();

    //! initialize the predictor
    bool initialize();

    //! Predict x1 using t0, x0, u and dt (alias between x0 and x1 allowed)
    void predict(const Eigen::Ref<const StateVector>& x0, std::vector<std::pair<double, ControlVector>> u_seq, double dt, Eigen::Ref<StateVector> x1);

    //! Set the system dynamics of the simulated plant
    void setSystemDynamics(SystemDynamicsInterface::Ptr dynamics);
    //! Set a numerical integrator for continuous-time dynamics
    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

#ifdef MESSAGE_SUPPORT
    //! Export to message
    virtual void toMessage(messages::OneStepPredictor& message) const;
    //! Import from message
    virtual void fromMessage(const messages::OneStepPredictor& message, std::stringstream* issues = nullptr);
#endif

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    bool _initialized = false;
};

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_ONE_STEP_PREDICTOR_H_

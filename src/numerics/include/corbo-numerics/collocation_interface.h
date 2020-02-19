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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_COLLOCATION_INTERFACE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_COLLOCATION_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>
#include <corbo-numerics/dynamics_eval_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/numerics/collocation.pb.h>
#endif

#include <functional>
#include <memory>

namespace corbo {

/**
 * @brief Interface class for collocation based system dynamics evaluation
 *
 * @ingroup numerics
 *
 * This base class provides an interface for approximating
 * continuous-time dynamics between two consecutive points.
 *
 * Online references:
 * http://www.control.lth.se/media/Education/DoctorateProgram/2011/OptimizationWithCasadi/lecture4b_short_slides.pdf
 * https://mec560sbu.github.io/2016/09/30/direct_collocation/
 * http://epubs.siam.org/doi/pdf/10.1137/16M1062569
 *
 * @remark This interface is provided with factory support of its parent class
 *         (DynamicsEvalInterface).
 *
 * @see ForwardDiffCollocation BackwardDiffCollocation MidpointDiffCollocation
 *      CrankNicolsonDiffCollocation CubicSplineCollocation DynamicsEvalInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class CollocationInterface : public DynamicsEvalInterface
{
 public:
    using Ptr  = std::shared_ptr<CollocationInterface>;
    using UPtr = std::unique_ptr<CollocationInterface>;

    using StateVector = Eigen::VectorXd;
    using InputVector = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~CollocationInterface() {}
    //! Return a newly allocated instances of the inherited class.
    DynamicsEvalInterface::Ptr getInstance() const override = 0;

    /**
     * @brief Compute differentiation error (system dynamics)
     *
     * \f[
     *      e =  f(x,u,t) - \dot{x}
     * \f].
     *
     * @param[in]  x1      Initial state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  u1      Constant control input vector [SystemDynamicsInterface::getInputDimension() x 1]
     * @param[in]  x2      Final state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  dt      Time interval length
     * @param[in]  system  System dynamics object
     * @param[out] error   Resulting error [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
     */
    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt,
                                   const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> error) override = 0;

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }

        states.push_back(x1);
        controls.push_back(u1);

        if (dt < 1e-15) return true;

        int n    = range.getNumInRange();
        double t = range.getStart();
        for (int i = 1; i < n; ++i, t += range.getStep())
        {
            double tau = (t - range.getStart()) / dt;
            states.push_back(x1 + tau * (x2 - x1));  // linear state interpolation by default
            controls.push_back(u1);                  // constant controls per default
        }
        if (range.includeEnd())
        {
            states.push_back(x2);
            controls.push_back(u1);
        }
        return true;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    virtual void toMessage(corbo::messages::Collocation& message) const {}
    // Implements interface method
    virtual void fromMessage(const corbo::messages::Collocation& message, std::stringstream* issues = nullptr) {}

    // Implements interface method
    void toMessage(corbo::messages::DynamicsEval& message) const override { toMessage(*message.mutable_collocation()); }
    // Implements interface method
    void fromMessage(const corbo::messages::DynamicsEval& message, std::stringstream* issues = nullptr) override
    {
        if (message.has_collocation()) fromMessage(message.collocation());
    }
#endif
};

#define FACTORY_REGISTER_COLLOCATION(type) FACTORY_REGISTER_DYNAMICS_EVAL(type)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_COLLOCATION_INTERFACE_H_

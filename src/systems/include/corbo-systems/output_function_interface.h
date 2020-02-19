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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_OUTPUT_FUNCTION_INTERFACE_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_OUTPUT_FUNCTION_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/systems/output_functions.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Interface class for system output functions
 *
 * @ingroup systems
 *
 * This class specifies methods that are required to be implemented by specific
 * subclasses in order to allow their general utilization in a variety of control tasks.
 *
 * A system output is defined by an algebraic equation
 * \f[
 *      y = c(x)
 * \f]
 * and hence maps the state vector of a SystemDynamicsInterface
 * to the actual system output.
 *
 * @remark This interface is provided with factory support (OutputFunctionFactory).
 *
 * @see SystemDynamicsInterface PlantInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Check if we need getStateDimension() in this interface
 */
class SystemOutputInterface
{
 public:
    using Ptr = std::shared_ptr<SystemOutputInterface>;

    using StateVector  = Eigen::VectorXd;
    using OutputVector = Eigen::VectorXd;

    //! Default destructor
    virtual ~SystemOutputInterface() = default;

    //! Get dimension of the system output y
    virtual int getOutputDimension() const = 0;

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    /**
     * @brief Evaluate the system output equation
     *
     * This method defines the mapping \f$ c : \mathbb{R}^p \to \mathbb{R}^q \f$
     * with \f$ p \f$ as the dimension of the state vector x and \f$ q \f$ the dimension
     * of the system output y.
     *
     * @param[in]  x  State vector [SystemDynamicsInterface::getStateDimension x 1]
     * @param[out] y  Output vector [getOutputDimension() x 1]
     */
    virtual void output(const StateVector& x, OutputVector& y) = 0;

    virtual void reset() {}

#ifdef MESSAGE_SUPPORT
    //! Export output function to message
    virtual void toMessage(corbo::messages::OutputFunction& message) const {}
    //! Import output function from message
    virtual void fromMessage(const corbo::messages::OutputFunction& message, std::stringstream* issues = nullptr) {}
#endif
};

using OutputFunctionFactory = Factory<SystemOutputInterface>;
#define FACTORY_REGISTER_OUTPUT_FUNCTION(type) FACTORY_REGISTER_OBJECT(type, SystemOutputInterface)

/**
 * @brief Return full state vector as system output
 *
 * @ingroup systems
 *
 * This output function is defined as \f$ y = x \f$.
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class FullStateSystemOutput : public SystemOutputInterface
{
 public:
    // implements interface method
    Ptr getInstance() const override { return std::make_shared<FullStateSystemOutput>(); }

    // implements interface method
    int getOutputDimension() const override { return property::INHERITED; }
    // implements interface method
    void output(const StateVector& x, OutputVector& y) override { y = x; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::OutputFunction& message) const override { message.mutable_full_state_system_output(); }
#endif
};
FACTORY_REGISTER_OUTPUT_FUNCTION(FullStateSystemOutput)

/**
 * @brief Return first state vector component as system output
 *
 * @ingroup systems
 *
 * This output function is defined as \f$ y = [1 0 \dotsc 0] x \f$.
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class FirstStateSystemOutput : public SystemOutputInterface
{
 public:
    // implements interface method
    Ptr getInstance() const override { return std::make_shared<FirstStateSystemOutput>(); }
    // implements interface method
    int getOutputDimension() const override { return 1; }
    // implements interface method
    void output(const StateVector& x, OutputVector& y) override { y[0] = x[0]; }
#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::OutputFunction& message) const override { message.mutable_first_state_system_output(); }
#endif
};
FACTORY_REGISTER_OUTPUT_FUNCTION(FirstStateSystemOutput)

/**
 * @brief Return first state vector component as system output
 *
 * @ingroup systems
 *
 * This output function is defined as \f$ y = [0 0 \dotsc 0 1] x \f$.
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class LastStateSystemOutput : public SystemOutputInterface
{
 public:
    // implements interface method
    Ptr getInstance() const override { return std::make_shared<LastStateSystemOutput>(); }
    // implements interface method
    int getOutputDimension() const override { return 1; }
    // implements interface method
    void output(const StateVector& x, OutputVector& y) override { y[0] = x.tail(1)[0]; }
#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::OutputFunction& message) const override { message.mutable_last_state_system_output(); }
#endif
};
FACTORY_REGISTER_OUTPUT_FUNCTION(LastStateSystemOutput)

/**
 * @brief Linear system output function
 *
 * @ingroup systems
 *
 * This output function is defined as \f$ y = C x \f$
 * with \f$ C \f$ as linear ouput matrix [SystemDynamicsInterface::getStateDimension() x getOutputDimension()].
 *
 * @see SystemDynamicsInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Message import/output not yet implemented
 */
class LinearStateSystemOutput : public SystemOutputInterface
{
 public:
    // implements interface method
    Ptr getInstance() const override { return std::make_shared<LastStateSystemOutput>(); }
    // implements interface method
    int getOutputDimension() const override { return 1; }
    // implements interface method
    void output(const StateVector& x, OutputVector& y) override { y.noalias() = _mat_c * x; }
#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::OutputFunction& message) const override { message.mutable_last_state_system_output(); }
#endif

    //! Get linear system matrix [SystemDynamicsInterface::getStateDimension() x getOutputDimension()]
    const Eigen::MatrixXd& getLinearMatrixC() const { return _mat_c; }
    //! Set linear system matrix [SystemDynamicsInterface::getStateDimension() x getOutputDimension()]
    void setLinearMatrixC(const Eigen::Ref<const Eigen::MatrixXd>& matric_c) { _mat_c = matric_c; }

 private:
    Eigen::MatrixXd _mat_c;
};
FACTORY_REGISTER_OUTPUT_FUNCTION(LinearStateSystemOutput)

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_OUTPUT_FUNCTION_INTERFACE_H_

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
#ifndef SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_INTERFACE_H_
#define SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_INTERFACE_H_

#include <corbo-core/global.h>
#include <corbo-core/signal_target_interface.h>
#include <corbo-tasks/environment.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/tasks/tasks.pb.h>
#endif

#include <memory>
#include <string>

namespace corbo {

/**
 * @brief Interface class for tasks
 *
 * @ingroup tasks
 *
 * Possible task implementations can be the closed-loop control of a plant,
 * open-loop control, benchmarking, ...
 *
 * Usually, tasks are called with an Environment to facilitate the
 * initialization and verification of commonly used control architectures.
 * An Environment contains a plant, observer and controller.
 * But a particular task does not need to rely on an Environment only but can contain
 * more objects or replace objects (e.g. multiple controllers, ...).
 *
 * @remark This interface is provided with factory support (TaskFactory).
 *
 * @see Environment
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class TaskInterface
{
 public:
    using Ptr = std::shared_ptr<TaskInterface>;

    //! Virtuel destructor
    virtual ~TaskInterface() {}
    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    /**
     * @brief Perform task
     * @param[in]      environment    Standard environment (plant, controller, observer)
     * @param[in,out]  signal_target  Target for occuring signals [optional]
     */
    virtual void performTask(Environment& environment, SignalTargetInterface* signal_target = nullptr, std::string* msg = nullptr,
                             const std::string& ns = "") = 0;

    /**
     * @brief Check if the environment and other settings satisfy all requirements for the given task
     *
     * This function can be called in order to check if all components and models are appropriate,
     * e.g. if all input and output dimensions are chosen adequately.
     *
     * Note, Environment::verify() might be invoked in order to check if controller, plant and observer
     * are specified correctly and if they have
     * matching dimensions.
     *
     * @param[in]      environment    Standard environment (plant, controller, observer)
     * @param[out]     msg            The string contains issue messages and hints if available [optional]
     * @returns true if verification was successfull, false otherwise.
     */
    virtual bool verify(const Environment& environment, std::string* msg = nullptr) const = 0;

    /**
     * @brief Retrieve available signals from the task
     *
     * Register a-priori known signals at the signal target.
     * Registration is optional.
     * Note, during performTask() further signals might occur without
     * registration (in case the they are not known in advance or the implementation lacks a proper registration).
     * @param[in,out] signal_target   Target for occuring signals [optional]
     */
    virtual void getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target, const std::string& ns = "") const {}

    //! Reset task state
    virtual void reset() = 0;

#ifdef MESSAGE_SUPPORT
    //! Export task settings to message
    virtual void toMessage(messages::Task& message) const {}
    //! Import task settings from message
    virtual void fromMessage(const messages::Task& message, std::stringstream* issues = nullptr) {}
#endif
};

using TaskFactory = Factory<TaskInterface>;
#define FACTORY_REGISTER_TASK(type) FACTORY_REGISTER_OBJECT(type, TaskInterface)

}  // namespace corbo

#endif  // SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_INTERFACE_H_

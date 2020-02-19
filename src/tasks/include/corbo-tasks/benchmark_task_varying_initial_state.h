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

#ifndef SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_VARYING_INITIAL_STATE_H_
#define SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_VARYING_INITIAL_STATE_H_

#include <corbo-tasks/task_interface.h>

#include <corbo-tasks/task_open_loop_control.h>

#include <memory>
#include <string>

namespace corbo {

/**
 * @brief BenchmarkTaskVaryingInitialState
 *
 * @ingroup tasks
 *
 * @see TaskInterface Environment
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class BenchmarkTaskVaryingInitialState : public TaskInterface
{
 public:
    using Ptr = std::shared_ptr<BenchmarkTaskVaryingInitialState>;

    //! Default constructor
    BenchmarkTaskVaryingInitialState();

    // implements interface method
    TaskInterface::Ptr getInstance() const override { return std::make_shared<BenchmarkTaskVaryingInitialState>(); }

    // implements interface method
    void performTask(Environment& environment, SignalTargetInterface* signal_target = nullptr, std::string* msg = nullptr,
                     const std::string& ns = "") override;

    // implements interface method
    bool verify(const Environment& environment, std::string* msg = nullptr) const override;

    // implements interface method
    void getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    // implements interface method
    void reset() override;

#ifdef MESSAGE_SUPPORT
    void toMessage(messages::BenchmarkTaskVaryingInitialState& message) const;
    void fromMessage(const messages::BenchmarkTaskVaryingInitialState& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Task& message) const override { toMessage(*message.mutable_benchmark_varying_initial_state()); }
    // implements interface method
    void fromMessage(const messages::Task& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.benchmark_varying_initial_state(), issues);
    }
#endif

 private:
    TaskInterface::Ptr _main_task;

    Eigen::VectorXd _x0_default;

    double _x01_start = 0;
    double _x01_end   = 0;
    int _x01_n        = 2;

    double _x02_start = 0;
    double _x02_end   = 0;
    int _x02_n        = 2;

    double _wait_time = 0;
};

FACTORY_REGISTER_TASK(BenchmarkTaskVaryingInitialState)

}  // namespace corbo

#endif  // SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_VARYING_INITIAL_STATE_H_

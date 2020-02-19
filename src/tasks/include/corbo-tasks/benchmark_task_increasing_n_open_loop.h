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

#ifndef SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_INCREASING_N_OPEN_LOOP_H_
#define SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_INCREASING_N_OPEN_LOOP_H_

#include <corbo-tasks/task_interface.h>

#include <corbo-tasks/task_open_loop_control.h>

#include <memory>
#include <string>

namespace corbo {

/**
 * @brief BenchmarkTaskIncreasingHorizonOpenLoop
 *
 * @ingroup tasks
 *
 * @remarks The underlying controller must support open-loop control!
 *
 * @see TaskInterface Environment
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class BenchmarkTaskIncreasingHorizonOpenLoop : public TaskInterface
{
 public:
    using Ptr = std::shared_ptr<BenchmarkTaskIncreasingHorizonOpenLoop>;

    //! Default constructor
    BenchmarkTaskIncreasingHorizonOpenLoop();

    // implements interface method
    TaskInterface::Ptr getInstance() const override { return std::make_shared<BenchmarkTaskIncreasingHorizonOpenLoop>(); }

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
    void toMessage(messages::BenchmarkTaskIncreasingHorizonOpenLoop& message) const;
    void fromMessage(const messages::BenchmarkTaskIncreasingHorizonOpenLoop& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Task& message) const override { toMessage(*message.mutable_benchmark_increasing_n_open_loop()); }
    // implements interface method
    void fromMessage(const messages::Task& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.benchmark_increasing_n_open_loop(), issues);
    }
#endif

 private:
    OpenLoopControlTask::Ptr _open_loop_task;

    int _n_start                     = 1;
    int _n_end                       = 50;
    int _n_step                      = 1;
    int _repetitions                 = 1;
    double _initial_tf               = 1.0;
    double _wait_time                = 1e-6;
    int _shooting_num_u_per_interval = 1;
    bool _publish_task_signals       = false;
};

FACTORY_REGISTER_TASK(BenchmarkTaskIncreasingHorizonOpenLoop)

}  // namespace corbo

#endif  // SRC_TASKS_INCLUDE_CORBO_TASKS_BENCHMARK_TASK_INCREASING_N_OPEN_LOOP_H_

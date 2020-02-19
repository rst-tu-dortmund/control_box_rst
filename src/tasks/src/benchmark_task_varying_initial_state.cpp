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

#include <corbo-tasks/benchmark_task_varying_initial_state.h>

#include <corbo-tasks/task_closed_loop_control.h>
#include <corbo-tasks/task_open_loop_control.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>
#include <corbo-core/signals.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>

#include <memory>
#include <string>
#include <thread>

namespace corbo {

BenchmarkTaskVaryingInitialState::BenchmarkTaskVaryingInitialState() {}

void BenchmarkTaskVaryingInitialState::getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target,
                                                           const std::string& ns) const
{
}

void BenchmarkTaskVaryingInitialState::performTask(Environment& environment, SignalTargetInterface* signal_target, std::string* msg,
                                                   const std::string& ns)
{
    if (!_main_task && _x01_n < 1) return;

    Eigen::VectorXd x0 = _x0_default;

    if (x0.size() != environment.getController()->getStateDimension())
    {
        if (msg) *msg += "BenchmarkTaskVaryingInitialState: dimension mismatch between x0_default and controller.\n";
    }

    if (x0.size() < 2)
    {
        _x02_n = 1;
    }

    double x01_dir  = _x01_end - _x01_start;
    double x01_step = _x01_n < 2 ? x01_dir : x01_dir / ((double)_x01_n - 1);
    double x02_dir  = _x02_end - _x02_start;
    double x02_step = _x02_n < 2 ? x02_dir : x02_dir / ((double)_x02_n - 1);

    Duration sleeper(_wait_time);

    int idx = 0;
    for (int i = 0; i < _x01_n; ++i)
    {
        x0[0] = _x01_start + (double)i * x01_step;

        for (int j = 0; j < _x02_n; ++j)
        {
            if (x0.size() > 1) x0[1] = _x02_start + (double)j * x02_step;

            PRINT_INFO("=============================");
            PRINT_INFO("run: " << idx + 1 << "/" << _x01_n * _x02_n);
            PRINT_INFO("=============================");

            std::string ns_bench = ns + "bench_" + std::to_string(idx) + "/";

            environment.reset();
            environment.getPlantPtr()->setState(x0);

            _main_task->reset();
            _main_task->performTask(environment, signal_target, msg, ns_bench);

            sleeper.sleep();

            ++idx;
        }
    }
}

bool BenchmarkTaskVaryingInitialState::verify(const Environment& environment, std::string* msg) const
{
    bool ret_val = true;

    if (!_main_task)
    {
        if (msg) *msg += "BenchmarkTaskVaryingInitialState(): no main task specified.\n";
        return false;
    }

    if (!_main_task->verify(environment, msg)) return false;

    int dim_x = environment.getController()->getStateDimension();

    if (_x01_n <= 0)
    {
        *msg += "_x01_n > 0 required.\n";
        ret_val = false;
    }

    if (_x02_n <= 0 && dim_x < 2)
    {
        *msg += "_x01_n > 0 required if state dimension is > 1.\n";
        ret_val = false;
    }

    if (_x0_default.size() != dim_x)
    {
        *msg += "Dimension of x0_default does not match state dimension.\n";
        ret_val = false;
    }

    return ret_val;
}

void BenchmarkTaskVaryingInitialState::reset()
{
    if (_main_task) _main_task->reset();
}

#ifdef MESSAGE_SUPPORT
void BenchmarkTaskVaryingInitialState::toMessage(corbo::messages::BenchmarkTaskVaryingInitialState& message) const
{
    //    if (_main_task) _main_task->toMessage(*message.mutable_main_task());
    //
    message.set_x01_start(_x01_start);
    message.set_x01_end(_x01_end);
    message.set_x01_n(_x01_n);

    message.set_x02_start(_x02_start);
    message.set_x02_end(_x02_end);
    message.set_x02_n(_x02_n);

    ClosedLoopControlTask::Ptr closed_loop_task = std::dynamic_pointer_cast<ClosedLoopControlTask>(_main_task);
    if (closed_loop_task)
        closed_loop_task->toMessage(*message.mutable_closed_loop_control_task());
    else
    {
        OpenLoopControlTask::Ptr open_loop_task = std::dynamic_pointer_cast<OpenLoopControlTask>(_main_task);
        if (open_loop_task) open_loop_task->toMessage(*message.mutable_open_loop_control_task());
    }

    // x0 default
    if (_x0_default.size() > 0)
    {
        message.mutable_x0_default()->Resize(_x0_default.size(), false);
        Eigen::Map<Eigen::VectorXd>(message.mutable_x0_default()->mutable_data(), _x0_default.size()) = _x0_default;
    }

    message.set_wait_time(_wait_time);
}
void BenchmarkTaskVaryingInitialState::fromMessage(const corbo::messages::BenchmarkTaskVaryingInitialState& message, std::stringstream* issues)
{
    switch (message.main_task_case())
    {
        case messages::BenchmarkTaskVaryingInitialState::kClosedLoopControlTask:
        {
            ClosedLoopControlTask::Ptr closed_loop_task = std::make_shared<ClosedLoopControlTask>();
            closed_loop_task->fromMessage(message.closed_loop_control_task(), issues);
            _main_task = closed_loop_task;
            break;
        }
        case messages::BenchmarkTaskVaryingInitialState::kOpenLoopControlTask:
        {
            OpenLoopControlTask::Ptr open_loop_task = std::make_shared<OpenLoopControlTask>();
            open_loop_task->fromMessage(message.open_loop_control_task(), issues);
            _main_task = open_loop_task;
            break;
        }
        default:
        {
            _main_task.reset();
            if (issues) *issues << "BenchmarkTaskVaryingInitialState: selected task not implemented.\n";
        }
    };

    _x01_start = message.x01_start();
    _x01_end   = message.x01_end();
    _x01_n     = message.x01_n();

    _x02_start = message.x02_start();
    _x02_end   = message.x02_end();
    _x02_n     = message.x02_n();

    _x0_default = Eigen::Map<const Eigen::VectorXd>(message.x0_default().data(), message.x0_default_size());

    _wait_time = message.wait_time();
}
#endif

}  // namespace corbo

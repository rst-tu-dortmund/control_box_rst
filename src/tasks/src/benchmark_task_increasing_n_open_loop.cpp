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

#include <corbo-tasks/benchmark_task_increasing_n_open_loop.h>

#include <corbo-communication/utilities.h>
#include <corbo-controllers/predictive_controller.h>
#include <corbo-core/console.h>
#include <corbo-core/signals.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_shooting_grid_base.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/shooting_grid_base.h>

#include <corbo-optimal-control/structured_ocp/structured_optimal_control_problem.h>

#include <memory>
#include <string>
#include <thread>

namespace corbo {

BenchmarkTaskIncreasingHorizonOpenLoop::BenchmarkTaskIncreasingHorizonOpenLoop() {}

void BenchmarkTaskIncreasingHorizonOpenLoop::getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target,
                                                                 const std::string& ns) const
{
}

void BenchmarkTaskIncreasingHorizonOpenLoop::performTask(Environment& environment, SignalTargetInterface* signal_target, std::string* msg,
                                                         const std::string& ns)
{
    if (!_open_loop_task || _n_end < _n_start || _repetitions < 1) return;

    StructuredOptimalControlProblem::Ptr ocp;
    PredictiveController::Ptr controller = std::dynamic_pointer_cast<PredictiveController>(environment.getControllerPtr());
    if (controller)
    {
        ocp = std::dynamic_pointer_cast<StructuredOptimalControlProblem>(controller->getOptimalControlProblem());
    }
    else
    {
        PRINT_WARNING("This benchmark is designed for predictive controllers.");
    }
    if (!ocp)
    {
        PRINT_ERROR(
            "We currently support only StructuredOptimalControlProblems.");  // TODO(roesmann): we just need an interface to setN() on the ocp!
    }

    Duration sleeper(_wait_time);

    std::vector<double> _controller_step_times;
    _controller_step_times.reserve(_repetitions);

    int n_step = std::max(1, _n_step);

    for (int n = _n_start; n <= _n_end && ok(); n = n + n_step)
    {
        PRINT_INFO("=============================");
        PRINT_INFO("n = " << n);
        PRINT_INFO("=============================");

        std::string ns_bench = ns + "bench_" + std::to_string(n) + "/";

        for (int i = 0; i < _repetitions && ok(); ++i)
        {
            PRINT_INFO_COND(i > 0 && i % 10 == 0, "Repetition no " << i + 1 << "/" << _repetitions);
            environment.reset();  // reset should not deallocate memory of anything that is available via the APIs.
            if (ocp)
            {
                DiscretizationGridInterface::Ptr grid = ocp->getDiscretizationGrid();
                grid->setN(n, false);
                if (_initial_tf >= 0)
                    grid->setInitialDt(_initial_tf / double(n - 1));
                else
                    PRINT_WARNING_ONCE("initial_tf: is negative, using default value of the chosen grid.");

                // check if we have a shooting grid
                ShootingGridBase::Ptr shooting_grid = std::dynamic_pointer_cast<ShootingGridBase>(grid);
                if (shooting_grid)
                {
                    shooting_grid->setNumControlsPerShootingInterval(
                        _shooting_num_u_per_interval);  // TODO(roesmann): do we really need this here?! (it's based on an old version)
                }
                else
                {
                    NonUniformShootingGridBase::Ptr nu_shooting_grid = std::dynamic_pointer_cast<NonUniformShootingGridBase>(grid);
                    if (nu_shooting_grid) nu_shooting_grid->setNumControlsPerShootingInterval(_shooting_num_u_per_interval);
                }
            }

            // only send all signals in the first reptition
            _open_loop_task->performTask(environment, (i == 0 && _publish_task_signals) ? signal_target : nullptr, msg, ns_bench);
            // get controller statistics
            ControllerStatistics::Ptr ctrl_statistics = environment.getController()->getStatistics();
            if (ctrl_statistics)
            {
                _controller_step_times.push_back(ctrl_statistics->step_time.toSec());
            }
        }
        signal_target->sendIndexedValues(ns + "ctrl_step_times", n, _controller_step_times);
        _controller_step_times.clear();
    }
}

bool BenchmarkTaskIncreasingHorizonOpenLoop::verify(const Environment& environment, std::string* msg) const
{
    bool ret_val = true;

    if (!_open_loop_task)
    {
        if (msg) *msg += "BenchmarkTaskIncreasingHorizonOpenLoop(): no OpenLoopControlTask specified.\n";
        return false;
    }

    if (!_open_loop_task->verify(environment, msg)) return false;

    if (_repetitions < 1)
    {
        *msg += "Number of repetitions must be positive.\n";
        ret_val = false;
    }

    if (_n_start < 2 && _n_end < 2)
    {
        *msg += "n_start and n_end must be > 1.\n";
        ret_val = false;
    }

    if (_n_end < _n_start)
    {
        *msg += "n_end >= n_start ist not satisfied.\n";
        ret_val = false;
    }

    return ret_val;
}

void BenchmarkTaskIncreasingHorizonOpenLoop::reset()
{
    if (_open_loop_task) _open_loop_task->reset();
}

#ifdef MESSAGE_SUPPORT
void BenchmarkTaskIncreasingHorizonOpenLoop::toMessage(corbo::messages::BenchmarkTaskIncreasingHorizonOpenLoop& message) const
{
    if (_open_loop_task) _open_loop_task->toMessage(*message.mutable_open_loop_control_task());

    message.set_n_start(_n_start);
    message.set_n_end(_n_end);
    message.set_n_step(_n_step);
    message.set_repetitions(_repetitions);
    message.set_shooting_num_u_per_interval(_shooting_num_u_per_interval);
    message.set_initial_tf(_initial_tf);
    message.set_wait_time(_wait_time);
    message.set_publish_task_signals(_publish_task_signals);
}
void BenchmarkTaskIncreasingHorizonOpenLoop::fromMessage(const corbo::messages::BenchmarkTaskIncreasingHorizonOpenLoop& message,
                                                         std::stringstream* issues)
{
    // open-loop task
    _open_loop_task = std::make_shared<OpenLoopControlTask>();
    _open_loop_task->fromMessage(message.open_loop_control_task(), issues);

    _n_start                     = message.n_start();
    _n_end                       = message.n_end();
    _n_step                      = message.n_step();
    _repetitions                 = message.repetitions();
    _shooting_num_u_per_interval = message.shooting_num_u_per_interval();
    _initial_tf                  = message.initial_tf();
    _wait_time                   = message.wait_time();
    _publish_task_signals        = message.publish_task_signals();
}
#endif

}  // namespace corbo

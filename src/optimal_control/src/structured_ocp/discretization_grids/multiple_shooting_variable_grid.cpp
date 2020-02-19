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

#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_variable_grid.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

void MultipleShootingVariableGrid::setDtBounds(double dt_lb, double dt_ub)
{
    _dt_lb = dt_lb;
    _dt_ub = dt_ub;
}

void MultipleShootingVariableGrid::setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio)
{
    _grid_adapt    = GridAdaptStrategy::TimeBasedSingleStep;
    _n_max         = n_max;
    _dt_hyst_ratio = dt_hyst_ratio;
}

void MultipleShootingVariableGrid::setGridAdaptTimeBasedAggressiveEstimate(int n_max, double dt_hyst_ratio)
{
    _grid_adapt    = GridAdaptStrategy::TimeBasedAggressiveEstimate;
    _n_max         = n_max;
    _dt_hyst_ratio = dt_hyst_ratio;
}

void MultipleShootingVariableGrid::setGridAdaptSimpleShrinkingHorizon() { _grid_adapt = GridAdaptStrategy::SimpleShrinkingHorizon; }

bool MultipleShootingVariableGrid::adaptGrid(bool new_run, NlpFunctions& nlp_fun)
{
    // do not adapt grid in a new run
    if (new_run && !_adapt_first_iter) return false;

    bool changed = false;
    switch (_grid_adapt)
    {
        case GridAdaptStrategy::NoGridAdapt:
        {
            break;
        }
        case GridAdaptStrategy::TimeBasedSingleStep:
        {
            changed = adaptGridTimeBasedSingleStep(nlp_fun);
            break;
        }
        case GridAdaptStrategy::TimeBasedAggressiveEstimate:
        {
            changed = adaptGridTimeBasedAggressiveEstimate(nlp_fun);
            break;
        }
        case GridAdaptStrategy::SimpleShrinkingHorizon:
        {
            changed = adaptGridSimpleShrinkingHorizon(nlp_fun);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("selected grid adaptation strategy not implemented.");
        }
    }
    return changed;
}

bool MultipleShootingVariableGrid::adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun)
{
    PRINT_WARNING_COND_NAMED(!isTimeVariableGrid(), "time based adaptation might only be used with a fixed dt.");

    int n = getN();

    double dt = getDt();
    if (dt > _dt_ref * (1.0 + _dt_hyst_ratio) && n < _n_max)
    {
        resampleTrajectory(n + 1, nlp_fun);
        _n_adapt = n + 1;
        return true;
    }
    else if (dt < _dt_ref * (1.0 - _dt_hyst_ratio) && n > _n_min)
    {
        resampleTrajectory(n - 1, nlp_fun);
        _n_adapt = n - 1;
        return true;
    }
    return false;
}

bool MultipleShootingVariableGrid::adaptGridTimeBasedAggressiveEstimate(NlpFunctions& nlp_fun)
{
    PRINT_WARNING_COND_NAMED(!isTimeVariableGrid(), "time based adaptation might only be used with a fixed dt.");

    int n     = getN();
    double dt = getDt();

    // check if hysteresis is satisfied
    if (dt >= _dt_ref * (1.0 - _dt_hyst_ratio) && dt <= _dt_ref * (1.0 + _dt_hyst_ratio)) return false;

    // estimate number of samples based on the fraction dt/dt_ref.
    // dt is the time difference obtained in a previous solution (with a coarser resp. finer trajectory resolution)
    int new_n = n * (int)std::round(_dt.value() / _dt_ref);

    // bound value
    if (new_n > _n_max)
        new_n = _n_max;
    else if (new_n < _n_min)
        new_n = _n_min;

    if (new_n == n) return false;

    // and now resample
    resampleTrajectory(new_n, nlp_fun);
    _n_adapt = new_n;
    return true;
}

bool MultipleShootingVariableGrid::adaptGridSimpleShrinkingHorizon(NlpFunctions& nlp_fun)
{
    int n = getN();
    if (n > _n_min)
    {
        resampleTrajectory(n - 1, nlp_fun);
        _n_adapt = n - 1;
    }
    return false;
}

#ifdef MESSAGE_SUPPORT
void MultipleShootingVariableGrid::fromMessage(const messages::MultipleShootingVariableGrid& message, std::stringstream* issues)
{
    MultipleShootingGrid::fromMessage(message.multiple_shooting_grid(), issues);

    // xf fixed states
    // if (grid_msg.xf_fixed_size() != _p && issues) *issues << "FullDiscretizationGrid: xf_fixed size does not match state dimension " << _p <<
    // ".\n";
    // TODO(roesmann): we need a warning in the gui if xf_fixed has the wrong dimension.
    //                 maybe we could add a "verifyParameters()" method to all interfaces, so predictive control asks the ocp, the ocp the grid and so
    //                 on
    _xf_fixed = Eigen::Map<const Eigen::Matrix<bool, -1, 1>>(message.xf_fixed().data(), message.xf_fixed_size());

    // dt bounds
    setDtBounds(message.dt_min(), message.dt_max());

    // auto resize
    if (message.has_grid_adapt_strategy())
    {
        if (message.grid_adapt_strategy().has_no_grid_adapt())
        {
            disableGridAdaptation();
        }
        else if (message.grid_adapt_strategy().has_time_based_single_step())
        {
            setGridAdaptTimeBasedSingleStep(message.grid_adapt_strategy().time_based_single_step().n_max(),
                                            message.grid_adapt_strategy().time_based_single_step().dt_hyst_ratio());
        }
        else if (message.grid_adapt_strategy().has_time_based_aggr_estim())
        {
            setGridAdaptTimeBasedAggressiveEstimate(message.grid_adapt_strategy().time_based_aggr_estim().n_max(),
                                                    message.grid_adapt_strategy().time_based_aggr_estim().dt_hyst_ratio());
        }
        else if (message.grid_adapt_strategy().has_simple_shrinking_horizon())
        {
            setGridAdaptSimpleShrinkingHorizon();
        }
    }
    _adapt_first_iter = message.grid_adapt_strategy().adapt_first_iter();

    // others
    _n_min = message.n_min();
}

void MultipleShootingVariableGrid::toMessage(messages::MultipleShootingVariableGrid& message) const
{
    MultipleShootingGrid::toMessage(*message.mutable_multiple_shooting_grid());

    // xf fixed states
    if (_xf_fixed.size() > 0)
    {
        message.mutable_xf_fixed()->Resize(_xf_fixed.size(), false);
        Eigen::Map<Eigen::Matrix<bool, -1, 1>>(message.mutable_xf_fixed()->mutable_data(), _xf_fixed.size()) = _xf_fixed;
    }

    // dt bounds
    message.set_dt_min(_dt_lb);
    message.set_dt_max(_dt_ub);

    // auto resize
    switch (_grid_adapt)
    {
        case GridAdaptStrategy::NoGridAdapt:
        {
            message.mutable_grid_adapt_strategy()->mutable_no_grid_adapt();
            break;
        }
        case GridAdaptStrategy::TimeBasedSingleStep:
        {
            message.mutable_grid_adapt_strategy()->mutable_time_based_single_step()->set_n_max(_n_max);
            message.mutable_grid_adapt_strategy()->mutable_time_based_single_step()->set_dt_hyst_ratio(_dt_hyst_ratio);
            break;
        }
        case GridAdaptStrategy::TimeBasedAggressiveEstimate:
        {
            message.mutable_grid_adapt_strategy()->mutable_time_based_aggr_estim()->set_n_max(_n_max);
            message.mutable_grid_adapt_strategy()->mutable_time_based_aggr_estim()->set_dt_hyst_ratio(_dt_hyst_ratio);
            break;
        }
        case GridAdaptStrategy::SimpleShrinkingHorizon:
        {
            message.mutable_grid_adapt_strategy()->mutable_simple_shrinking_horizon();
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("exporting of the selected auto resize strategy not implemented.");
        }
    }
    message.mutable_grid_adapt_strategy()->set_adapt_first_iter(_adapt_first_iter);

    // others
    message.set_n_min(_n_min);
}
#endif

}  // namespace corbo

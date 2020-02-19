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
 *  Authors: Maximilian Krämer, Christoph Rösmann
 *********************************************************************/

#include <corbo-optimal-control/structured_ocp/discretization_grids/full_discretization_grid_move_blocking_base.h>

namespace corbo {

void FullDiscretizationGridMoveBlockingBase::getStateAndControlTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max) const
{
    if (x_sequence) x_sequence->clear();
    if (u_sequence) u_sequence->clear();

    if (isEmpty()) return;
    assert(isValid());

    PRINT_ERROR_COND_NAMED(t_max < 0, "t_max >= 0 required");

    double dt = getDt();

    if (x_sequence)
    {
        double t = 0;
        for (int i = 0; i < _x_seq.size(); ++i)
        {
            x_sequence->add(t, _x_seq[i].values());
            t += dt;
            if (t > t_max) break;
        }
        if (t <= t_max) x_sequence->add(t, _xf.values());
    }

    if (u_sequence)
    {
        double t = 0;
        for (int i = 0; i < _u_seq.size(); ++i)
        {
            // Use blocking vector to reconstruct valid time series
            for (int j = 0; j < _blocking_vector(i); ++j)
            {
                u_sequence->add(t, _u_seq[i].values());
                t += dt;
                if (t > t_max) break;
            }
        }
        // duplicate last u to have the sampe time stamps as x_sequence
        if (t <= t_max) u_sequence->add(t, _u_seq.back().values());
    }
}

bool FullDiscretizationGridMoveBlockingBase::isValid() const
{
    // Check for consistent blocking vector and correct, implicit u sequence
    return (_u_seq.size() == _blocking_vector.size()) && (_blocking_vector.sum() == _x_seq.size());
}

void FullDiscretizationGridMoveBlockingBase::initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf,
                                                                 ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun)
{
    // clear();  // make sure everything is cleared
    _x_seq.clear();
    _u_seq.clear();
    _xf.clear();
    _active_vertices.clear();

    // check or initalize bounds if empty
    nlp_fun.checkAndInitializeBoundDimensions(x0.size(), uref.getDimension());

    // check x_fixed
    checkAndInitializeXfFixedFlags(x0.size());

    int n_init = _n_adapt > 0 ? _n_adapt : _n_ref;

    int num_intervals = n_init - 1;

    // we initialize the state trajectory linearly
    Eigen::VectorXd dir = xf - x0;
    double dist         = dir.norm();
    if (dist != 0) dir /= dist;
    double step = dist / num_intervals;

    for (int k = 0; k < num_intervals; ++k)
    {
        // add new state by linear interpolation
        _x_seq.emplace_back(x0 + (double)k * step * dir, nlp_fun.x_lb, nlp_fun.x_ub);
    }

    int ref_id = 0;
    for (int k = 0; k < _blocking_vector.size(); ++k)
    {
        // add new control according to blocking vector
        _u_seq.emplace_back(uref.getReferenceCached(ref_id), nlp_fun.u_lb, nlp_fun.u_ub);
        ref_id += _blocking_vector(k);
    }

    // add final state
    _xf.set(xf, nlp_fun.x_lb, nlp_fun.x_ub, _xf_fixed);

    // set first state as fixed
    _x_seq.front().setFixed(true);

    // set dt
    _dt.set(_dt_ref, _dt_lb, _dt_ub, isDtFixedIntended());

    assert(_x_seq.size() > 0);

    // notify graph that vertices are changed
    setModified(true);
}

void FullDiscretizationGridMoveBlockingBase::initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf,
                                                                 ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                                                 NlpFunctions& nlp_fun)
{
    // clear();  // make sure everything is cleared
    _x_seq.clear();
    _u_seq.clear();
    _xf.clear();
    _active_vertices.clear();

    // TODO(roesmann): check if this is ever needed in any pracital realization:
    // if ((x0 - xref.getReferenceCached(0)).norm() > _non_static_ref_dist_threshold)
    // {
    //    initializeSequences(x0, xref.getReferenceCached(getN() - 1), uref,  nlp_fun);
    // }

    // check or initalize bounds if empty
    nlp_fun.checkAndInitializeBoundDimensions(x0.size(), uref.getDimension());

    // check x_fixed
    checkAndInitializeXfFixedFlags(x0.size());

    int n_init = _n_adapt > 0 ? _n_adapt : _n_ref;

    int num_intervals = n_init - 1;

    _x_seq.emplace_back(x0, nlp_fun.x_lb, nlp_fun.x_ub);  // always add the (exact) x0
    for (int k = 1; k < num_intervals; ++k)
    {
        // add new state by linear interpolation
        _x_seq.emplace_back(xref.getReferenceCached(k), nlp_fun.x_lb, nlp_fun.x_ub);
    }

    int ref_id = 0;
    for (int k = 0; k < _blocking_vector.size(); ++k)
    {
        // add new control according to blocking vector
        _u_seq.emplace_back(uref.getReferenceCached(ref_id), nlp_fun.u_lb, nlp_fun.u_ub);
        ref_id += _blocking_vector(k);
    }

    // add final state
    _xf.set(xf, nlp_fun.x_lb, nlp_fun.x_ub, _xf_fixed);

    // set first state as fixed
    _x_seq.front().setFixed(true);

    // set dt
    _dt.set(_dt_ref, _dt_lb, _dt_ub, isDtFixedIntended());

    assert(_x_seq.size() > 1);

    // notify graph that vertices are changed
    setModified(true);
}

void FullDiscretizationGridMoveBlockingBase::warmStartShifting(const Eigen::VectorXd& x0)
{
    // find nearest state to x0 (ideally it is the second one in _x_seq).
    int num_shift = findNearestState(x0);
    if (num_shift <= 0) return;

    if (num_shift > getN() - 2)
    {
        PRINT_ERROR_NAMED("Cannot shift if num_shift > N-2");
        return;
    }

    // the simplest strategy would be to just remove remove x_seq[0], append xf to x_seq and replace xf (num_shift=1)
    // however, this requries to change the structure which always triggers an edge recreation and so the solver must reallocate memory.
    // TOOD(roesmann): for time-optimal stuff, if we are sure that we resize the grid, it is then more efficent to implement the strategy above...

    // shift from end to front:
    for (int i = 0; i < getN() - num_shift; ++i)
    {
        int idx = i + num_shift;
        if (idx == getN() - 1)
        {
            // final state reached
            _x_seq[i].values() = _xf.values();
        }
        else
        {
            _x_seq[i].values() = _x_seq[idx].values();
        }
    }

    if (_warm_start_shift_u)
    {
        int u_idx     = 0;
        int block_idx = 0;
        int i         = 0;
        int k         = 0;

        for (; block_idx < _blocking_vector.size(); ++block_idx)
        {
            for (; k < _blocking_vector.size(); ++k)
            {
                if (u_idx + _blocking_vector(k) > num_shift + i)
                {
                    _u_seq[block_idx].values() = _u_seq[k].values();
                    break;
                }
                else
                {
                    u_idx += _blocking_vector(k);
                }
            }
            if (k == _blocking_vector.size())
            {
                break;
            }

            i += _blocking_vector(block_idx);
        }
        for (int j = block_idx; j < _blocking_vector.size(); ++j)
        {
            _u_seq[j].values() = _u_seq[j - 1].values();
        }
    }

    int idx = getN() - num_shift;
    if (idx < 0)
    {
        PRINT_ERROR_NAMED("idx < 0...");
        return;
    }
    for (int i = 0; i < num_shift; ++i, ++idx)
    {
        // now extrapolate
        assert(idx > 1);

        // linearly extrapolate states
        if (i == num_shift - 1)  // consider xf
        {
            _xf.values() = _x_seq[idx - 2].values() + 2.0 * (_x_seq[idx - 1].values() - _x_seq[idx - 2].values());
        }
        else
        {
            _x_seq[idx].values() = _x_seq[idx - 2].values() + 2.0 * (_x_seq[idx - 1].values() - _x_seq[idx - 2].values());
            // TODO(roesmann) multiply by fraction of last dt and _dt
        }
    }
}

void FullDiscretizationGridMoveBlockingBase::resampleTrajectory(int n_new) { PRINT_ERROR_NAMED("Not supported for move blocking"); }

void FullDiscretizationGridMoveBlockingBase::computeActiveVertices()
{
    _active_vertices.clear();
    assert(isValid());
    int n = getN();

    int u_idx     = 0;
    int block_idx = 0;
    for (int i = 0; i < n - 1; ++i)
    {
        if (!_x_seq[i].isFixed()) _active_vertices.push_back(&_x_seq[i]);
        if ((i == u_idx) && (!_u_seq[block_idx].isFixed()))
        {
            _active_vertices.push_back(&_u_seq[block_idx]);

            u_idx += _blocking_vector(block_idx);
            ++block_idx;  // TODO(kraemer) increment on fixed
        }
    }

    if (!_xf.isFixed()) _active_vertices.push_back(&_xf);
    if (!_dt.isFixed()) _active_vertices.push_back(&_dt);
}

}  // namespace corbo

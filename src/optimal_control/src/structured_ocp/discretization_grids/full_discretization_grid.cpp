/*********************************************************************
 *
 * Software License Agreement
 *
 *  Copyright (c) 2018,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This software is currently not released.
 *  Redistribution and use in source and binary forms,
 *  with or without modification, are prohibited.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Authors: Christoph Rösmann
 *********************************************************************/

#include <corbo-optimal-control/structured_ocp/discretization_grids/full_discretization_grid.h>

#include <corbo-core/console.h>

#include <memory>

namespace corbo {

bool FullDiscretizationGrid::initialize(const Eigen::VectorXd& x0, const Eigen::VectorXd& u0, DynamicsEvalInterface::Ptr dynamics_eval)
{
    return initialize(x0, x0, u0, dynamics_eval);
}

bool FullDiscretizationGrid::initialize(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf, const Eigen::VectorXd& u0,
                                        DynamicsEvalInterface::Ptr /*dynamics_eval*/)
{
    clear();

    DiscretizationGrid::initialize(x0.size(), u0.size());

    if (_num_desired_states < 2)
    {
        PRINT_ERROR("FullDiscretizationGrid::initialize(): number of desired states > 1 required. Please call setHorizon() first.");
        return false;
    }

    int num_intervals = _num_desired_states - 1;

    // we initialize the state trajectory linearly
    Eigen::VectorXd dir = xf - x0;
    double dist         = dir.norm();
    if (dist != 0) dir /= dist;
    double step = dist / num_intervals;

    for (int k = 0; k < num_intervals; ++k)
    {
        // add new shooting interval
        appendShootingInterval(x0 + (double)k * step * dir, u0, _dt);
    }

    // fix start state for optimization
    getShootingIntervalsRaw().front().shooting_node->setFixed(true);

    PartiallyFixedVectorVertex::Ptr xf_node = std::make_shared<PartiallyFixedVectorVertex>(xf, _x_lower, _x_upper);
    if (_xf_fixed.size() == 0)
        xf_node->setFixed(false);
    else
        xf_node->setFixed(_xf_fixed);
    getFinalStateRaw() = xf_node;

    // updateBounds();  // TODO(roesmann): more efficienct

    return true;
}

bool FullDiscretizationGrid::initialize(const Eigen::VectorXd& x0, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                        const Time& t, DynamicsEvalInterface::Ptr dynamics_eval)
{
    DiscretizationGrid::initialize(x0.size(), uref.getDimension());

    clear();

    if (xref.isStatic() && uref.isStatic())
    {
        return initialize(x0, xref.getReferenceCached(0), uref.getReferenceCached(0), dynamics_eval);
    }

    else
    {
        if (_num_desired_states < 2)
        {
            PRINT_ERROR("FullDiscretizationGrid::initialize(): number of desired states > 1 required. Please call setHorizon() first.");
            return false;
        }

        if (xref.isStatic())  // xref.isStatic && !uref.isStatic
        {
            uref.precompute(_dt, _num_desired_states, t);

            clear();

            DiscretizationGrid::initialize(xref.getDimension(), uref.getDimension());

            int num_intervals = _num_desired_states - 1;

            // we initialize the state trajectory linearly
            Eigen::VectorXd dir = xref.getReferenceCached(num_intervals) - x0;
            double dist         = dir.norm();
            if (dist != 0) dir /= dist;
            double step = dist / num_intervals;

            for (int k = 0; k < num_intervals; ++k)
            {
                // add new shooting interval
                appendShootingInterval(x0 + (double)k * step * dir, uref.getReferenceCached(k), _dt);
            }

            // fix start state for optimization
            getShootingIntervalsRaw().front().shooting_node->setFixed(true);

            PartiallyFixedVectorVertex::Ptr xf_node =
                std::make_shared<PartiallyFixedVectorVertex>(xref.getReferenceCached(num_intervals), _x_lower, _x_upper);
            if (_xf_fixed.size() == 0)
                xf_node->setFixed(false);
            else
                xf_node->setFixed(_xf_fixed);
            getFinalStateRaw() = xf_node;
        }

        else  // !xref.isStatic && !uref.isStatic || !xref.isStatic && uref.isStatic
        {
            xref.precompute(_dt, _num_desired_states, t);
            uref.precompute(_dt, _num_desired_states, t);

            if ((x0 - xref.getReferenceCached(0)).norm() > _non_static_ref_dist_threshold)
                return initialize(x0, xref.getReferenceCached(_num_desired_states - 1), uref.getReferenceCached(_num_desired_states - 1),
                                  dynamics_eval);
            else
            {
                clear();

                DiscretizationGrid::initialize(xref.getDimension(), uref.getDimension());

                int num_intervals = _num_desired_states - 1;

                appendShootingInterval(x0, uref.getReferenceCached(0), _dt);

                for (int k = 1; k < num_intervals; ++k)
                {
                    // add new shooting interval
                    appendShootingInterval(xref.getReferenceCached(k), uref.getReferenceCached(k), _dt);
                }

                // fix start state for optimization
                getShootingIntervalsRaw().front().shooting_node->setFixed(true);

                PartiallyFixedVectorVertex::Ptr xf_node =
                    std::make_shared<PartiallyFixedVectorVertex>(xref.getReferenceCached(num_intervals), _x_lower, _x_upper);
                if (_xf_fixed.size() == 0)
                    xf_node->setFixed(false);
                else
                    xf_node->setFixed(_xf_fixed);
                getFinalStateRaw() = xf_node;
            }
        }
    }

    return true;
}

bool FullDiscretizationGrid::update(const Eigen::VectorXd& x0, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, bool new_run,
                                    const Time& t)
{
    if (_shooting_intervals.empty() || !_warm_start)
    {
        if (!new_run)
        {
            resizeTrajectory();  // TODO(roesman)
            setN(getN());
        }
        initialize(x0, xref, uref, t);
    }
    else if (new_run)
    {
        // get current goal
        Eigen::VectorXd xf(x0.rows());
        xref.getReference(Time(getFinalTime()), xf);

        // reinit if goal is far away
        if (_current_xf_values.rows() > 0 && (xf - _current_xf_values).norm() > _warm_start_goal_dist_reinit)
        {
            initialize(x0, xref, uref, t);  // initialize should also clear everything first
            _current_xf_values = xf;
        }
        else
        {
            // warm-start!

            // update goal
            // _xf->values() = xf;

            if (!_prune_trajectory)
            {
                // shift values and extrapolate
                if (!shiftGridValues(x0)) PRINT_WARNING_NAMED("shifting of GridValues failed!");
            }
            else
            {
                pruneTrajectory(x0);
            }

            if (_shooting_intervals.empty())
            {
                PRINT_ERROR_NAMED("warm-start failed. No grid intervals found. Reinitializing...");
                initialize(x0, xref, uref, t);
            }

            // update start
            _shooting_intervals.front().shooting_node->values() = x0;
            // int num_erased = pruneTrajectory(x0);

            // recreate erased shooting nodes at the end of the horizon if prune_trajectory is turned off)
            // if (!_prune_trajectory) extrapolateTrajectory(num_erased);

            // store xf values, we might decide later whether to reinitialize (even in case we have a moving horizon controller)
            if (_current_xf_values.rows() == 0) _current_xf_values = xf;
        }

        if (!new_run) resizeTrajectory();
    }

    return isModified();
}

int FullDiscretizationGrid::pruneTrajectory(const Eigen::VectorXd& x0)
{
    int nearest_idx = findNewInitialShootingInterval(x0);

    // prune trajectory at the beginning
    if (nearest_idx > 0)
    {
        // nearest_idx is equal to the number of samples to be removed (since it counts from 0 ;-) )
        _shooting_intervals.erase(_shooting_intervals.begin(),
                                  _shooting_intervals.begin() + nearest_idx);  // delete first states such that the closest state is the new first one

        // fix new start state for optimization
        getShootingIntervalsRaw().front().shooting_node->setFixed(true);

        setModified(true);
    }

    // update start
    _shooting_intervals.front().shooting_node->values() = x0;

    return nearest_idx;
}

void FullDiscretizationGrid::appendShootingInterval(const Eigen::VectorXd& x, const Eigen::VectorXd& u, double dt_if_not_single)
{
    getShootingIntervalsRaw().emplace_back();
    ShootingInterval& interval = getShootingIntervalsRaw().back();

    // initialize x
    ShootingVertex::Ptr xk_node = std::make_shared<ShootingVertex>(x, _x_lower, _x_upper);  // set first state (k == 0) to fixed
    interval.shooting_node      = xk_node;

    // in full discretization we always have a single control and dt per interval
    interval.num_controls = 1;

    // initialize u
    interval.controls.push_back(std::make_shared<ControlVertex>(u, _u_lower, _u_upper));

    // initialize dt
    if (_single_dt)
    {
        if (_shooting_intervals.size() == 1)  // this is the first interval
            interval.dts.push_back(std::make_shared<DtVertex>(_dt, _dt_lower, _dt_upper, _dt_fixed));
        else
            interval.dts.push_back(getFirstDtVertex());
    }
    else
    {
        interval.dts.push_back(std::make_shared<DtVertex>(dt_if_not_single, _dt_lower, _dt_upper, _dt_fixed));
    }

    setModified(true);
}

void FullDiscretizationGrid::insertShootingInterval(int k, const Eigen::VectorXd& x, const Eigen::VectorXd& u, double dt_if_not_single)
{
    assert(k > 0);
    assert(k <= (int)_shooting_intervals.size());

    if (k == _shooting_intervals.size())
    {
        appendShootingInterval(x, u, dt_if_not_single);
        return;
    }

    auto new_interv_it         = _shooting_intervals.insert(_shooting_intervals.begin() + k, ShootingInterval());
    ShootingInterval& interval = *new_interv_it;

    // initialize x
    ShootingVertex::Ptr xk_node = std::make_shared<ShootingVertex>(x, _x_lower, _x_upper);  // set first state (k == 0) to fixed
    interval.shooting_node      = xk_node;

    // in full discretization we always have a single control and dt per interval
    interval.num_controls = 1;

    // initialize u
    interval.controls.push_back(std::make_shared<ControlVertex>(u, _u_lower, _u_upper));

    // initialize dt
    if (_single_dt)
    {
        if (_shooting_intervals.size() == 1)  // this is the first interval
            interval.dts.push_back(std::make_shared<DtVertex>(_dt, _dt_lower, _dt_upper, _dt_fixed));
        else
            interval.dts.push_back(getFirstDtVertex());
    }
    else
    {
        interval.dts.push_back(std::make_shared<DtVertex>(dt_if_not_single, _dt_lower, _dt_upper, _dt_fixed));
    }

    setModified(true);
}

void FullDiscretizationGrid::eraseShootingInterval(int k)
{
    assert(k > 0);
    assert(k < (int)_shooting_intervals.size());
    _shooting_intervals.erase(_shooting_intervals.begin() + k);

    if (k == 0 && !_shooting_intervals.empty())
    {
        // fix new start state for optimization
        getShootingIntervalsRaw().front().shooting_node->setFixed(true);
    }
    setModified(true);
}

void FullDiscretizationGrid::extrapolateTrajectory(int num_new_intervals)
{
    if (num_new_intervals < 1) return;

    if (_shooting_intervals.empty())
    {
        PRINT_ERROR("FullDiscretizationGrid::extrapolateTrajectory(): cannot extrapolate with less than states");  // at least one interval + xf
        return;
    }

    // add current final state as new shooting interval
    appendShootingInterval(_xf->values(), _shooting_intervals.back().controls.front()->values(), _dt);

    // now add further intervals
    for (int i = 0; i < num_new_intervals; ++i)
    {
        // add new shooting interval
        // TODO(roesmann): switch to zero-order hold close to goal
        auto it_back       = std::prev(_shooting_intervals.end());
        Eigen::VectorXd xi = it_back->shooting_node->values() +
                             (it_back->shooting_node->values() -
                              std::prev(it_back)->shooting_node->values());  // TODO(roesmann) multiply by fraction of last dt and _dt

        if (i < num_new_intervals - 1)
        {
            appendShootingInterval(xi, it_back->controls.front()->values(), _dt);
        }
        else
        {
            // update xf
            _xf->values() = xi;
        }
    }
}

bool FullDiscretizationGrid::shiftGridValues(const Eigen::VectorXd& x0)
{
    int nearest_idx = findNewInitialShootingInterval(x0);
    return shiftGridValues(nearest_idx);
}

bool FullDiscretizationGrid::shiftGridValues(int num_shift)
{
    if (num_shift == 0)
        return true;
    else if (num_shift < 0)
        return false;

    if (num_shift > getN() - 2)
    {
        PRINT_ERROR_NAMED("Cannot shift if num_shift > N-2");
        return false;
    }

    // shift from end to front:
    for (int i = 0; i < getN() - num_shift; ++i)
    {
        int idx = i + num_shift;
        if (idx == getN() - 1)
        {
            assert(_xf);
            // final state reached
            _shooting_intervals[i].shooting_node->values() = _xf->values();
        }
        else
        {
            _shooting_intervals[i].shooting_node->values() = _shooting_intervals[idx].shooting_node->values();

            for (int j = 0; j < _shooting_intervals[i].controls.size(); ++j)
            {
                assert(_shooting_intervals[i].controls.size() == _shooting_intervals[idx].controls.size());
                _shooting_intervals[i].controls[j]->values() = _shooting_intervals[idx].controls[j]->values();
            }

            for (int j = 0; j < _shooting_intervals[i].dts.size(); ++j)
            {
                assert(_shooting_intervals[i].dts.size() == _shooting_intervals[idx].dts.size());
                _shooting_intervals[i].dts[j]->values() = _shooting_intervals[idx].dts[j]->values();
            }

            assert(_shooting_intervals[i].additional_states.empty() && "No additional states allowed in full discretization grid.");

            // _shooting_intervals[i].dt_complete_interval = _shooting_intervals[idx].dt_complete_interval; // should be always the same here
            // _shooting_intervals[i].num_controls         = _shooting_intervals[idx].num_controls; // should be always 1
        }
    }

    int idx = getN() - num_shift;
    if (idx < 0)
    {
        PRINT_ERROR_NAMED("idx < 0...");
        return false;
    }
    for (int i = 0; i < num_shift; ++i, ++idx)
    {
        // now extrapolate
        assert(idx > 1);

        // linearly extrapolate states
        if (i == num_shift - 1)  // consider xf
        {
            _xf->values() = _shooting_intervals[idx - 2].shooting_node->values() +
                            2.0 * (_shooting_intervals[idx - 1].shooting_node->values() - _shooting_intervals[idx - 2].shooting_node->values());
        }
        else
        {
            _shooting_intervals[idx].shooting_node->values() =
                _shooting_intervals[idx - 2].shooting_node->values() +
                2.0 * (_shooting_intervals[idx - 1].shooting_node->values() - _shooting_intervals[idx - 2].shooting_node->values());
            // TODO(roesmann) multiply by fraction of last dt and _dt
        }

        for (int j = 0; j < _shooting_intervals[idx - 1].controls.size(); ++j)
        {
            assert(_shooting_intervals[idx - 1].controls.size() == _shooting_intervals[idx - 2].controls.size());
            _shooting_intervals[idx - 1].controls[j]->values() = _shooting_intervals[idx - 2].controls[j]->values();
        }

        for (int j = 0; j < _shooting_intervals[idx - 1].dts.size(); ++j)
        {
            assert(_shooting_intervals[idx - 1].dts.size() == _shooting_intervals[idx - 2].dts.size());
            _shooting_intervals[idx - 1].dts[j]->values() = _shooting_intervals[idx - 2].dts[j]->values();
        }

        assert(_shooting_intervals[i].additional_states.empty() && "No additional states allowed in full discretization grid.");
    }

    return true;
}

int FullDiscretizationGrid::findNewInitialShootingInterval(const Eigen::VectorXd& x0)
{
    assert(!_shooting_intervals.empty());
    assert(_shooting_intervals.front().shooting_node);

    // same start as before
    if (x0 == _shooting_intervals.front().shooting_node->values()) return 0;

    // find nearest state (using l2-norm) in order to prune the trajectory
    // (remove already passed states, should be 1 if the sampling times of the planning and controller call match, but could be higher if rates
    // differ)

    int num_interv = (int)_shooting_intervals.size();

    double dist_cache = (x0 - _shooting_intervals.front().shooting_node->values()).norm();
    double dist;
    int num_keep_interv = std::max(1, _num_states_min - 1);            // we need to keep at least this number of intervals in the trajectory
    int lookahead       = std::min(num_interv - num_keep_interv, 20);  // max 20 samples

    int nearest_idx = 0;

    for (int i = 1; i <= lookahead; ++i)
    {
        dist = (x0 - _shooting_intervals[i].shooting_node->values()).norm();
        if (dist < dist_cache)
        {
            dist_cache  = dist;
            nearest_idx = i;
        }
        else
            break;
    }

    return nearest_idx;
}

bool FullDiscretizationGrid::resizeTrajectory()
{
    bool success = true;
    switch (_auto_resize)
    {
        case AutoResizeStrategy::NoAutoResize:
        {
            break;
        }
        case AutoResizeStrategy::TimeBased:
        {
            success = resizeTrajectoryTimeBased();
            break;
        }
        case AutoResizeStrategy::RedundantControls:
        {
            success = resizeTrajectoryRedundantControls();
            break;
        }
        default:
        {
            PRINT_ERROR("FullDiscretizationGrid::resizeTrajectory(): selected auto resize strategy not implemented.");
            success = false;
        }
    }
    return success;
}

bool FullDiscretizationGrid::resizeTrajectoryTimeBased()
{
    if (isDtFixed())
    {
        PRINT_WARNING("FullDiscretizationGrid::resizeTrajectoryTimeBased(): time based resize might only be used with an unfixed dt.");
    }

    int n = getN();

    if (hasSingleDt())
    {
        assert(getFirstDtVertexRaw());
        double dt = getFirstDtVertexRaw()->value();
        //    if (cfg->ctrl.estimate_dt)
        //    {
        //        // estimate number of samples based on the fraction dt/dt_ref.
        //        // dt is the time difference obtained in a previous solution (with a coarser resp. finer trajectory resolution)
        //        int new_n = n * (int)std::round(_dt.dt() / _dt_ref);
        //        new_n     = bound(cfg->ctrl.n_min, new_n, cfg->ctrl.n_max);
        //        resampleTrajectory(new_n);
        //        PRINT_INFO("new n: " << new_n << " old n: " << n);
        //    }
        // else  // use a simple linear search
        //{
        if (dt > _dt * (1.0 + _dt_hyst_ratio) && n < _num_states_max)
        {
            resampleTrajectory(n + 1);
        }
        else if (dt < _dt * (1.0 - _dt_hyst_ratio) && n > _num_states_min)
        {
            resampleTrajectory(n - 1);
        }
        //}
    }
    else
    {
        // iterate sequence once and check time differences
        bool changed = false;

        for (int i = 0; i < (int)_shooting_intervals.size() - 1; ++i)
        {
            ShootingInterval& interv      = _shooting_intervals[i];
            ShootingInterval& interv_next = _shooting_intervals[i + 1];

            double dt = interv.dts.front()->value();

            if (dt > _dt * (1.0 + _dt_hyst_ratio) && n < _num_states_max)
            {
                double new_dt = 0.5 * dt;

                insertShootingInterval(i + 1, 0.5 * (interv.shooting_node->values() + interv_next.shooting_node->values()),
                                       interv.controls.front()->values(), new_dt);
                changed = true;

                break;  // TODO(roesmann): better strategy, instead of just searching the first dt?

                ++i;  // skip newly inserted node
            }
            else if (dt < _dt * (1.0 - _dt_hyst_ratio) && n > _num_states_min)
            {
                // if (i < (int)_dt_seq.size()-2)
                {
                    interv_next.dts.front()->value() += dt;
                    eraseShootingInterval(i);

                    changed = true;

                    break;  // TODO(roesmann): better strategy, instead of just searching the first dt?
                }
            }
        }
    }
    return true;
}

bool FullDiscretizationGrid::resizeTrajectoryRedundantControls()
{
    // find approx(u_k, u_{k+1})
    std::vector<std::size_t> non_unique_indices;
    for (std::size_t idx = 0; idx < (int)_shooting_intervals.size() - 1; ++idx)  // never delete the last control
    {
        ShootingInterval& interv      = _shooting_intervals[idx];
        ShootingInterval& interv_next = _shooting_intervals[idx + 1];

        // if ( _ctrl_seq[idx-1].controls().isApprox(_ctrl_seq[idx].controls(), 1e-1 ) )
        // also delete if the time diff of the interval is sufficiently small
        if (interv.dts.front()->value() < 1e-6)
        {
            non_unique_indices.emplace_back(idx);
            // non_unique_indices.emplace_back(idx+1); // since we have zero dt, at least two samples are obsolete
            // ++idx;
            continue;
        }
        // PRINT_INFO("control threshold: " << _controls_similar_threshold);
        // if ( (_ctrl_seq[idx+1].controls()-_ctrl_seq[idx].controls()).isZero(_controls_similar_threshold) )
        //     non_unique_indices.emplace_back(idx);
        Eigen::VectorXd& u_k   = interv.controls.front()->values();
        Eigen::VectorXd& u_kp1 = interv_next.controls.front()->values();
        if (((u_kp1 - u_k).cwiseAbs().array() <= _redundant_ctrl_epsilon).all()) non_unique_indices.emplace_back(idx);
    }

    //   PRINT_INFO("number non-unique: " << non_unique_indices.size());

    bool changed = false;

    int backup_diff = (int)non_unique_indices.size() - _redundant_ctrl_backup;

    if (backup_diff < 0)
    {
        changed     = true;
        backup_diff = std::abs(backup_diff);
        for (int i = 0; i < backup_diff && getN() < _num_states_max; ++i)
        {
            // add new sample
            int dt_max_idx = 0;
            if (getN() > 2)
            {
                // insert inbetween largest gap (largest DT)
                auto end_it = _shooting_intervals.end();
                std::advance(end_it, -1);  // we do not want to a new point to the goal (therefore -2)
                auto max_elem = std::max_element(_shooting_intervals.begin(), end_it, [](const ShootingInterval& si1, const ShootingInterval& si2) {
                    return si1.dts.front()->value() < si2.dts.front()->value();
                });
                if (max_elem == end_it)
                {
                    PRINT_INFO("Invalid time max element in resizeTrajectoryRedundantControls(). break...");
                    break;
                }
                dt_max_idx = std::distance(_shooting_intervals.begin(), max_elem);
            }
            assert(dt_max_idx + 1 < (int)_shooting_intervals.size());
            ShootingInterval& interv      = _shooting_intervals[dt_max_idx];
            ShootingInterval& interv_next = _shooting_intervals[dt_max_idx + 1];

            double new_dt               = 0.5 * interv.dts.front()->value();
            interv.dts.front()->value() = new_dt;
            insertShootingInterval(dt_max_idx + 1, 0.5 * (interv.shooting_node->values() + interv_next.shooting_node->values()),
                                   interv.controls.front()->values(), new_dt);
            // PRINT_INFO("inserted");
        }
    }
    else if (backup_diff > 0)
    {
        changed     = true;
        auto idx_it = non_unique_indices.rbegin();  // erase starting from the last one (reverse iterator)
        for (int i = 0; i < backup_diff && getN() > _num_states_min; ++i)
        {
            int k = (int)*idx_it;
            if (k >= getN() - 2) --k;

            assert(k + 1 < (int)_shooting_intervals.size());
            ShootingInterval& interv      = _shooting_intervals[k];
            ShootingInterval& interv_next = _shooting_intervals[k + 1];

            interv.dts.front()->value() += interv_next.dts.front()->value();
            eraseShootingInterval(k + 1);
            ++idx_it;
            // PRINT_INFO("removed");
        }
    }
    return true;
}

void FullDiscretizationGrid::resampleTrajectory(int n_new)
{
    int n = getN();
    if (n == n_new) return;

    // copy vertices (vertices itself are shared pointers)
    //! @todo(roesmann) More efficient strategy without copying containers at all?
    //!
    TimeSeries::Ptr ts_states_old   = std::make_shared<TimeSeries>();
    TimeSeries::Ptr ts_controls_old = std::make_shared<TimeSeries>();
    getShootingNodeTimeSeries(ts_states_old);
    getControlInputTimeSeries(ts_controls_old);

    TimeSeries::ValuesMatMap states_old   = ts_states_old->getValuesMatrixView();
    TimeSeries::ValuesMatMap controls_old = ts_controls_old->getValuesMatrixView();

    int num_interv = n - 1;

    if (hasSingleDt())
    {
        assert(getFirstDtVertexRaw() != nullptr);
        double dt_old = getFirstDtVertexRaw()->value();
        // compute new time diff
        double dt_new = dt_old * double(n - 1) / double(n_new - 1);

        double t_new;
        int idx_old     = 1;
        double t_old_p1 = dt_old;  // time for old sample with index idx_old (acutally its the subsequent time step w.r.t t_new)

        for (int idx_new = 1; idx_new < n_new - 1; ++idx_new)  // n_new-1 since last state is identical and copied later. idx_new=1, since start
                                                               // sample is already valid (we do not touch it)
        // we allow a small mismatch for the control u_1 and let the optimizer correct it later
        {
            t_new = dt_new * double(idx_new);
            while (t_new > double(idx_old) * dt_old && idx_old < n)
            {
                ++idx_old;
            };  // find idx_old that represents the state subsequent to the new one (w.r.t. time)
            t_old_p1 = double(idx_old) * dt_old;

            const Eigen::VectorXd& x_prev = states_old.col(idx_old - 1);
            const Eigen::VectorXd& x_cur  = (idx_old < n - 1) ? states_old.col(idx_old) : _xf->values();

            if (idx_new < num_interv)
            {
                // states / shooting node
                _shooting_intervals[idx_new].shooting_node->values() = x_prev + (t_new - (t_old_p1 - dt_old)) / dt_old * (x_cur - x_prev);

                // controls
                // TODO(roesmann): we do not have a valid control (u_f), so hold last
                _shooting_intervals[idx_new].controls.front()->values() = controls_old.col(idx_old - 1);
            }
            else
            {
                // add new shooting interval
                appendShootingInterval(x_prev + (t_new - (t_old_p1 - dt_old)) / dt_old * (x_cur - x_prev), controls_old.col(idx_old - 1), _dt);
            }
        }

        // clear invalid states
        if (n_new < n)
        {
            _shooting_intervals.resize(n_new - 1);
        }

        // save new dt
        getFirstDtVertexRaw()->value() = dt_new;

        // set first vertex always to fixed
        // _shooting_intervals.front().shooting_node->setFixed(true); // should not be necessary here (we just change values in the first node)
    }
    else
    {
        PRINT_ERROR("FullDiscretizationGrid::resampleTrajectory(): not yet implemented for multiple dts");
    }

    // notify graph that vertices are changed
    setModified(true);
}

double FullDiscretizationGrid::getFinalTime() const
{
    if (_dt_fixed && _single_dt) return (double)_shooting_intervals.size() * _dt;
    if (_single_dt) return (double)_shooting_intervals.size() * getFirstDt();

    return DiscretizationGrid::getFinalTime();
}

void FullDiscretizationGrid::setGridResizeTimeBased(int num_states_max, double dt_hyst_ratio)
{
    _auto_resize    = AutoResizeStrategy::TimeBased;
    _num_states_max = num_states_max;
    _dt_hyst_ratio  = dt_hyst_ratio;
}

void FullDiscretizationGrid::setGridResizeRedundControls(int num_states_max, int num_backup_nodes, double epsilon)
{
    _auto_resize            = AutoResizeStrategy::RedundantControls;
    _num_states_max         = num_states_max;
    _redundant_ctrl_backup  = num_backup_nodes;
    _redundant_ctrl_epsilon = epsilon;
}

void FullDiscretizationGrid::computeActiveVertices()
{
    _active_vertices.clear();
    if (hasSingleDt())
    {
        for (ShootingInterval& interv : _shooting_intervals)
        {
            if (!interv.shooting_node->isFixed()) _active_vertices.push_back(interv.shooting_node.get());
            if (!interv.controls.front()->isFixed()) _active_vertices.push_back(interv.controls.front().get());
        }
        if (!_xf->isFixed()) _active_vertices.push_back(_xf.get());
        if (!getFirstDtVertexRaw()->isFixed()) _active_vertices.push_back(getFirstDtVertexRaw());
    }
    else
    {
        for (ShootingInterval& interv : _shooting_intervals)
        {
            if (!interv.shooting_node->isFixed()) _active_vertices.push_back(interv.shooting_node.get());
            if (!interv.controls.front()->isFixed()) _active_vertices.push_back(interv.controls.front().get());
            if (!interv.dts.front()->isFixed()) _active_vertices.push_back(interv.dts.front().get());
        }
        if (!_xf->isFixed()) _active_vertices.push_back(_xf.get());
    }
}

void FullDiscretizationGrid::updateDtsFixedFlag()
{
    if (_single_dt && !_shooting_intervals.empty() && !_shooting_intervals.front().dts.empty())
    {
        getFirstDtVertexRaw()->setFixed(_dt_fixed);
    }
    else
    {
        for (ShootingInterval& interv : _shooting_intervals)
        {
            for (ScalarVertex::Ptr& dt_vtx : interv.dts)
            {
                dt_vtx->setFixed(_dt_fixed);
            }
        }
    }
}

void FullDiscretizationGrid::updateDts()
{
    // compute average dt
    double dt = 0;

    for (const ShootingInterval& interv : _shooting_intervals) dt += interv.getTotalTimeInterval();

    if (!_shooting_intervals.empty()) dt /= (double)_shooting_intervals.size();

    if (_single_dt)
    {
        ScalarVertex::Ptr single_dt = std::make_shared<ScalarVertex>(dt, _dt_fixed);
        for (ShootingInterval& interv : _shooting_intervals)
        {
            for (ScalarVertex::Ptr& dt_vtx : interv.dts)
            {
                dt_vtx = single_dt;
            }
        }
    }
    else
    {
        for (ShootingInterval& interv : _shooting_intervals)
        {
            if (interv.dts.empty()) continue;
            double local_dt = dt / (double)interv.dts.size();
            for (ScalarVertex::Ptr& dt_vtx : interv.dts)
            {
                dt_vtx = std::make_shared<ScalarVertex>(local_dt, _dt_fixed);
            }
        }
    }
}

#ifdef MESSAGE_SUPPORT
void FullDiscretizationGrid::fromMessage(const messages::FullDiscretizationGrid& message, std::stringstream* issues)
{
    if (message.n() < 2 && issues) *issues << "FullDiscretizationGrid: Number of states must be greater than or equal 2.\n";

    if (message.dt() <= 0 && issues) *issues << "FullDiscretizationGrid: Dt must be greater than 0.0.\n";

    setN(message.n());
    setDt(message.dt());

    // dt
    setSingleDt(message.single_dt());
    setDtFixed(message.fixed_dt());

    int dim_states   = 0;
    int dim_controls = 0;

    // xf fixed states
    // if (grid_msg.xf_fixed_size() != _p && issues) *issues << "FullDiscretizationGrid: xf_fixed size does not match state dimension " << _p <<
    // ".\n";
    _xf_fixed  = Eigen::Map<const Eigen::Matrix<bool, -1, 1>>(message.xf_fixed().data(), message.xf_fixed_size());
    dim_states = _xf_fixed.rows();

    // bounds
    if (message.x_min_size() != dim_states && issues)
    {
        *issues << "FullDiscretizationGrid: x_min size does not match dimension of xf_fixed " << dim_states << ".\n";
    }
    else
    {
        setLowerStateBounds(Eigen::Map<const Eigen::VectorXd>(message.x_min().data(), message.x_min_size()));
    }

    if (message.x_max_size() != dim_states && issues)
    {
        *issues << "FullDiscretizationGrid: x_max size does not match dimension of xf_fixed " << dim_states << ".\n";
    }
    else
    {
        setUpperStateBounds(Eigen::Map<const Eigen::VectorXd>(message.x_max().data(), message.x_max_size()));
    }

    // if (message.u_min_size() != dim_controls && issues)
    //    *issues << "FullDiscretizationGrid: u_min size does not match control input dimension " << dim_controls << ".\n";
    // else
    setLowerControlBounds(Eigen::Map<const Eigen::VectorXd>(message.u_min().data(), message.u_min_size()));
    dim_controls = message.u_min_size();
    if (message.u_max_size() != dim_controls && issues)
        *issues << "FullDiscretizationGrid: u_max size does not match dimension of u_min " << dim_controls << ".\n";
    // else
    setUpperControlBounds(Eigen::Map<const Eigen::VectorXd>(message.u_max().data(), message.u_max_size()));

    // dt bounds
    setDtBounds(message.dt_min(), message.dt_max());

    // auto resize
    if (message.has_resize_strategy())
    {
        if (message.resize_strategy().has_no_auto_resize())
        {
            disableGridResize();
        }
        else if (message.resize_strategy().has_time_based())
        {
            setGridResizeTimeBased(message.resize_strategy().time_based().n_max(), message.resize_strategy().time_based().dt_hyst_ratio());
        }
        else if (message.resize_strategy().has_redundant_controls())
        {
            setGridResizeRedundControls(message.resize_strategy().redundant_controls().n_max(),
                                        message.resize_strategy().redundant_controls().num_backup(),
                                        message.resize_strategy().redundant_controls().epsilon());
        }
    }

    // others
    _warm_start       = message.warm_start();
    _prune_trajectory = message.prune_trajectory();
    _num_states_min   = message.n_min();

    _include_intermediate_constr = message.intermediate_collocation_constr();

    DiscretizationGrid::initialize(dim_states, dim_controls);
}

void FullDiscretizationGrid::toMessage(messages::FullDiscretizationGrid& message) const
{
    message.set_n(_num_desired_states);
    message.set_dt(_dt);

    // dt
    message.set_fixed_dt(_dt_fixed);
    message.set_single_dt(_single_dt);

    // xf fixed states
    if (_xf_fixed.size() > 0)
    {
        message.mutable_xf_fixed()->Resize(_xf_fixed.size(), false);
        Eigen::Map<Eigen::Matrix<bool, -1, 1>>(message.mutable_xf_fixed()->mutable_data(), _xf_fixed.size()) = _xf_fixed;
    }

    // bounds
    message.mutable_x_min()->Resize(_x_lower.size(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_x_min()->mutable_data(), _x_lower.size()) = _x_lower;
    message.mutable_x_max()->Resize(_x_upper.size(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_x_max()->mutable_data(), _x_upper.size()) = _x_upper;
    message.mutable_u_min()->Resize(_u_lower.size(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_u_min()->mutable_data(), _u_lower.size()) = _u_lower;
    message.mutable_u_max()->Resize(_u_upper.size(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_u_max()->mutable_data(), _u_upper.size()) = _u_upper;

    // dt bounds
    message.set_dt_min(_dt_lower);
    message.set_dt_max(_dt_upper);

    // auto resize
    switch (_auto_resize)
    {
        case AutoResizeStrategy::NoAutoResize:
        {
            message.mutable_resize_strategy()->mutable_no_auto_resize();
            break;
        }
        case AutoResizeStrategy::TimeBased:
        {
            message.mutable_resize_strategy()->mutable_time_based()->set_n_max(_num_states_max);
            message.mutable_resize_strategy()->mutable_time_based()->set_dt_hyst_ratio(_dt_hyst_ratio);
            break;
        }
        case AutoResizeStrategy::RedundantControls:
        {
            message.mutable_resize_strategy()->mutable_redundant_controls()->set_n_max(_num_states_max);
            message.mutable_resize_strategy()->mutable_redundant_controls()->set_num_backup(_redundant_ctrl_backup);
            message.mutable_resize_strategy()->mutable_redundant_controls()->set_epsilon(_redundant_ctrl_epsilon);
            break;
        }
        default:
        {
            PRINT_ERROR("FullDiscretizationGrid::toMessage(): exporting of the selected auto resize strategy not implemented.");
        }
    }

    // others
    message.set_n_min(_num_states_min);
    message.set_warm_start(_warm_start);
    message.set_prune_trajectory(_prune_trajectory);
    message.set_intermediate_collocation_constr(_include_intermediate_constr);
}
#endif

}  // namespace corbo

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

#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_full_discretization_grid_base.h>

#include <corbo-optimal-control/structured_ocp/edges/finite_differences_collocation_edges.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

GridUpdateResult NonUniformFullDiscretizationGridBase::update(const Eigen::VectorXd& x0, ReferenceTrajectoryInterface& xref,
                                                              ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun, OptimizationEdgeSet& edges,
                                                              SystemDynamicsInterface::Ptr dynamics, bool new_run, const Time& t,
                                                              ReferenceTrajectoryInterface* sref, const Eigen::VectorXd* prev_u, double prev_u_dt,
                                                              ReferenceTrajectoryInterface* xinit, ReferenceTrajectoryInterface* uinit)
{
    assert(x0.size() == dynamics->getStateDimension());
    assert(xref.getDimension() == x0.size());
    assert(uref.getDimension() == dynamics->getInputDimension());

    GridUpdateResult result;

    if (isGridAdaptActive() && !_first_run && !isEmpty())  // TODO(roesmann): this might not be efficient in case warm start is deactivated
    {
        // adapt grid if implemented by subclass
        adaptGrid(new_run, nlp_fun);  // note, isModified() should return true if something was changed!, we could also use the return value...
    }

    int n = std::max(getNRef(), _n_adapt);

    // check if we need to cache the reference trajectory values // TODO(roesmann): we could restrict this to new_run==true only as long as we do not
    // have a grid resize...
    std::vector<double> dts;
    if (!xref.isStatic() || !uref.isStatic() || (sref && !sref->isStatic()) || (xinit && !xinit->isStatic()) || (uinit && !uinit->isStatic()))
    {
        if (!_dt_seq.empty() && _n_adapt == 0)  // TODO(roesmann): check if we add a more robust check
            getDts(dts);
        else
            dts.resize(_n_ref, _dt_ref);
        if (!xref.isStatic() && !xref.isCached(dts, t))
            xref.precompute(dts, t);  // TODO(roesmann): if we have grid adaptation, we might need to chache up to n_max....
        if (!uref.isStatic() && !uref.isCached(dts, t)) uref.precompute(dts, t);

        if (sref && !sref->isStatic() && !sref->isCached(dts, t)) sref->precompute(dts, t);
        if (xinit && !xinit->isStatic() && !xinit->isCached(dts, t)) xinit->precompute(dts, t);
        if (uinit && !uinit->isStatic() && !uinit->isCached(dts, t)) uinit->precompute(dts, t);
    }

    // set previous control which might be utilized by some edges
    if (prev_u && prev_u->size() > 0)
        setPreviousControl(*prev_u, prev_u_dt);
    else
        setPreviousControl(Eigen::VectorXd::Zero(dynamics->getInputDimension()), prev_u_dt);

    // set last control to uref
    setLastControlRef(uref.getReferenceCached(_n_ref - 1));

    // TODO(roesmann) we do not check if bounds in nlp_fun are updated
    // updateBounds(); // calling this everytime is not efficient

    // initialize trajectories if empty or if the goal differs from the last one
    // reinit if goal is far away
    if (isEmpty() ||
        !_warm_start)  // TODO(roesmann): threshold?: if (_cached_xf.rows() > 0 && (xf - _cached_xf).norm() > _warm_start_goal_dist_reinit)
    {
        if (xref.isStatic() && !xinit)
        {
            initializeSequences(x0, xref.getReferenceCached(_n_ref - 1), uinit ? *uinit : uref, nlp_fun);
        }
        else
        {
            initializeSequences(x0, xref.getReferenceCached(_n_ref - 1), xinit ? *xinit : xref, uinit ? *uinit : uref, nlp_fun);
        }
    }
    else
    {
        if (new_run && isMovingHorizonWarmStartActive())
        {
            // warm_start
            warmStartShifting(x0);
        }
        if (new_run)
        {
            // make sure to always overwrite the start with the actual (measured) values
            _x_seq.front().values() = x0;
            // update fixed goal states
            for (int i = 0; i < _xf_fixed.size(); ++i)
            {
                if (_xf_fixed[i]) _xf.values()[i] = xref.getReferenceCached(getN() - 1)[i];
            }
        }
    }

    result.vertices_updated = isModified();  // isModified() returns if the underlying vertex set is updated

    if (new_run || result.updated())  // new run -> new t
    {
        // update NLP functions w.r.t. the current discretization
        result.edges_updated =
            nlp_fun.update(getN(), t.toSec(), xref, uref, sref, hasSingleDt(), x0, dts, this);  // returns true if edge dimensions changed
        // TODO(roesmann): sref is not yet implemented for testing, add this later...
        // TODO(roesmann): we do not yet check if dt was updated, this might affect nlp update
        // TODO(roesmann): maybe, for non-uniform grids, we need a better nlp update method which gets all correct reference time stamps
    }

    if (result.updated())  // vertices or eges updated
    {
        // create grid specific edges
        // TODO(roesmann): for now, we assume that we are the only maintainer of the edge set (which might be not true generally!)
        //                 so we clear the edge set whenever we want
        createEdges(nlp_fun, edges, dynamics);
        result.edges_updated = true;  // now we definitely updated the edgeset
    }

    _first_run = false;
    return result;
}

void NonUniformFullDiscretizationGridBase::initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf,
                                                               ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun)
{
    // clear();  // make sure everything is cleared
    _x_seq.clear();
    _u_seq.clear();
    _dt_seq.clear();
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
        // add new constant control
        _u_seq.emplace_back(uref.getReferenceCached(k), nlp_fun.u_lb, nlp_fun.u_ub);
        // add new time interval
        _dt_seq.emplace_back(_dt_ref, _dt_lb, _dt_ub, isDtFixedIntended());
    }
    // add final state
    _xf.set(xf, nlp_fun.x_lb, nlp_fun.x_ub, _xf_fixed);

    // set first state as fixed
    _x_seq.front().setFixed(true);

    assert(_x_seq.size() > 1);

    // notify graph that vertices are changed
    setModified(true);
}

void NonUniformFullDiscretizationGridBase::initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf,
                                                               ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                                               NlpFunctions& nlp_fun)
{
    // clear();  // make sure everything is cleared
    _x_seq.clear();
    _u_seq.clear();
    _dt_seq.clear();
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
    _u_seq.emplace_back(uref.getReferenceCached(0), nlp_fun.u_lb, nlp_fun.u_ub);
    _dt_seq.emplace_back(_dt_ref, _dt_lb, _dt_ub, isDtFixedIntended());
    for (int k = 1; k < num_intervals; ++k)
    {
        // add new state by linear interpolation
        _x_seq.emplace_back(xref.getReferenceCached(k), nlp_fun.x_lb, nlp_fun.x_ub);
        // add new constant control
        _u_seq.emplace_back(uref.getReferenceCached(k), nlp_fun.u_lb, nlp_fun.u_ub);
        // add new time interval
        _dt_seq.emplace_back(_dt_ref, _dt_lb, _dt_ub, isDtFixedIntended());
    }
    // add final state
    _xf.set(xf, nlp_fun.x_lb, nlp_fun.x_ub, _xf_fixed);

    // set first state as fixed
    _x_seq.front().setFixed(true);

    assert(_x_seq.size() > 1);

    // notify graph that vertices are changed
    setModified(true);
}

void NonUniformFullDiscretizationGridBase::warmStartShifting(const Eigen::VectorXd& x0)
{
    // TODO(roesmann) check if the following code works and is meaningful for a non-uniform grid!! I did just a quick modification..

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
            _x_seq[i].values()  = _x_seq[idx].values();
            _u_seq[i].values()  = _u_seq[idx].values();
            _dt_seq[i].values() = _dt_seq[idx].values();
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
        _u_seq[idx - 1].values()  = _u_seq[idx - 2].values();
        _dt_seq[idx - 1].values() = _dt_seq[idx - 2].values();
    }
}

int NonUniformFullDiscretizationGridBase::findNearestState(const Eigen::VectorXd& x0)
{
    assert(!isEmpty());
    assert(isValid());

    // same start as before
    double first_dist = (x0 - _x_seq.front().values()).norm();
    if (std::abs(first_dist) < 1e-12) return 0;

    // find nearest state (using l2-norm) in order to prune the trajectory
    // (remove already passed states, should be 1 if the sampling times of the planning and controller call match, but could be higher if rates
    // differ)

    int num_interv = getN() - 1;

    double dist_cache = first_dist;
    double dist;
    int num_keep_interv = 1;  // std::max(1, _num_states_min - 1);            // we need to keep at least this number of intervals in the trajectory
    int lookahead       = std::min(num_interv - num_keep_interv, 20);  // max 20 samples

    int nearest_idx = 0;

    for (int i = 1; i <= lookahead; ++i)
    {
        dist = (x0 - _x_seq[i].values()).norm();
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

bool NonUniformFullDiscretizationGridBase::getFirstControlInput(Eigen::VectorXd& u0)
{
    if (isEmpty() || !isValid()) return false;
    assert(!_u_seq.empty());

    u0 = _u_seq.front().values();
    return true;
}

void NonUniformFullDiscretizationGridBase::getDts(std::vector<double>& dts) const
{
    dts.resize(_dt_seq.size());
    for (int i = 0; i < _dt_seq.size(); ++i) dts[i] = _dt_seq[i].value();
}

double NonUniformFullDiscretizationGridBase::getFinalTime() const
{
    double t = 0;
    for (const ScalarVertex& dtv : _dt_seq) t += dtv.value();
    return t;
}

void NonUniformFullDiscretizationGridBase::setNRef(int n)
{
    // clear grid if we change n
    if (n != getN()) clear();
    if (n < 2)
    {
        PRINT_ERROR_NAMED("Number of states must be n>1.");
        _n_ref = 2;
        return;
    }
    _n_ref   = n;
    _n_adapt = 0;
}

bool NonUniformFullDiscretizationGridBase::checkAndInitializeXfFixedFlags(int dim_x)
{
    if (_xf_fixed.size() == 0)
    {
        _xf_fixed.setConstant(dim_x, false);
        return true;
    }
    else if (_xf_fixed.size() == dim_x)
        return true;

    PRINT_ERROR_NAMED("Dimensions mismatch between xf_fixed and xf. Setting xf_fixed to false.");
    _xf_fixed.setConstant(dim_x, false);
    return false;
}

void NonUniformFullDiscretizationGridBase::updateBounds(const NlpFunctions& nlp_fun)
{
    if (isEmpty()) return;
    for (VectorVertex& vtx : _x_seq)
    {
        if (vtx.getDimension() == nlp_fun.x_lb.size())
            vtx.setLowerBounds(nlp_fun.x_lb);
        else
            PRINT_ERROR_NAMED("Cannot update lower state bounds due to dimensions mismatch");

        if (vtx.getDimension() == nlp_fun.x_ub.size())
            vtx.setUpperBounds(nlp_fun.u_ub);
        else
            PRINT_ERROR_NAMED("Cannot update upper state bounds due to dimensions mismatch");
    }

    for (VectorVertex& vtx : _u_seq)
    {
        if (vtx.getDimension() == nlp_fun.u_lb.size())
            vtx.setLowerBounds(nlp_fun.u_lb);
        else
            PRINT_ERROR_NAMED("Cannot update lower control input bounds due to dimensions mismatch");

        if (vtx.getDimension() == nlp_fun.u_ub.size())
            vtx.setUpperBounds(nlp_fun.u_ub);
        else
            PRINT_ERROR_NAMED("Cannot update upper control input bounds due to dimensions mismatch");
    }

    for (ScalarVertex& vtx : _dt_seq)
    {
        vtx.setLowerBound(_dt_lb);
        vtx.setUpperBound(_dt_ub);
    }
}

void NonUniformFullDiscretizationGridBase::clear()
{
    _x_seq.clear();
    _u_seq.clear();
    _dt_seq.clear();
    _xf.clear();
    _active_vertices.clear();
    _first_run = true;
    _n_adapt   = 0;
    setModified(true);
}

void NonUniformFullDiscretizationGridBase::getVertices(std::vector<VertexInterface*>& vertices)
{
    vertices.clear();
    // order doesn't matter here
    for (VectorVertex& vtx : _x_seq) vertices.push_back(&vtx);
    vertices.push_back(&_xf);
    for (VectorVertex& vtx : _u_seq) vertices.push_back(&vtx);
    for (ScalarVertex& vtx : _dt_seq) vertices.push_back(&vtx);
    vertices.push_back(&_u_prev);  // always fixed...
    vertices.push_back(&_u_ref);
    vertices.push_back(&_u_prev_dt);
}

void NonUniformFullDiscretizationGridBase::computeActiveVertices()
{
    _active_vertices.clear();
    assert(isValid());
    int n = getN();

    for (int i = 0; i < n - 1; ++i)
    {
        if (!_x_seq[i].isFixed()) _active_vertices.push_back(&_x_seq[i]);
        if (!_u_seq[i].isFixed()) _active_vertices.push_back(&_u_seq[i]);
        if (!_dt_seq[i].isFixed()) _active_vertices.push_back(&_dt_seq[i]);
    }
    if (!_xf.isFixed()) _active_vertices.push_back(&_xf);
}

void NonUniformFullDiscretizationGridBase::getStateAndControlTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max) const
{
    if (x_sequence) x_sequence->clear();
    if (u_sequence) u_sequence->clear();

    if (isEmpty()) return;
    assert(isValid());

    PRINT_ERROR_COND_NAMED(t_max < 0, "t_max >= 0 required");

    if (x_sequence)
    {
        double t = 0;
        for (int i = 0; i < _dt_seq.size(); ++i)
        {
            x_sequence->add(t, _x_seq[i].values());
            t += _dt_seq[i].value();
            if (t > t_max) break;
        }
        if (t <= t_max) x_sequence->add(t, _xf.values());
    }

    if (u_sequence)
    {
        double t = 0;
        for (int i = 0; i < _dt_seq.size(); ++i)
        {
            u_sequence->add(t, _u_seq[i].values());
            t += _dt_seq[i].value();
            if (t > t_max) break;
        }
        // duplicate last u to have the sampe time stamps as x_sequence
        if (t <= t_max) u_sequence->add(t, _u_seq.back().values());
    }
}

}  // namespace corbo

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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_REFERENCE_TRAJECTORY_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_REFERENCE_TRAJECTORY_H_

#include <corbo-core/factory.h>
#include <corbo-core/signals.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/core/reference_trajectories.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Interface class for reference trajectories
 *
 * @ingroup core
 *
 * This class represents a generic reference trajectory, e.g. for plant controllers.
 * The interface does not distinguish between state and control input references,
 * but is subject to an Eigen::VectorXd type of arbitrary dimension.
 *
 * @remark This interface is provided with factory support (ReferenceTrajectoryFactory).
 *
 * @see StaticReference ZeroReference
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo The interface is not yet completed/stable for non-static references,
 *       e.g. an appropriate pre-computation and time-shifting is not designed yet.
 */
class ReferenceTrajectoryInterface
{
 public:
    using Ptr          = std::shared_ptr<ReferenceTrajectoryInterface>;
    using OutputVector = Eigen::VectorXd;

    virtual Ptr getInstance() const = 0;

    virtual ~ReferenceTrajectoryInterface() = default;

    //! Get access to the associated factory
    static Factory<ReferenceTrajectoryInterface>& getFactory() { return Factory<ReferenceTrajectoryInterface>::instance(); }

    virtual bool isStatic() const = 0;
    virtual bool isZero() const { return false; }
    virtual int getDimension() const = 0;

    virtual void precompute(double dt, int n, Time t)              = 0;
    virtual void precompute(const std::vector<double>& dt, Time t) = 0;

    virtual void getReference(const Time& t, OutputVector& ref) const = 0;

    virtual const OutputVector& getReferenceCached(int k) const = 0;

    virtual const OutputVector& getNextSteadyState(const Time& t) = 0;

    virtual bool isCached(double dt, int n, Time t) const              = 0;
    virtual bool isCached(const std::vector<double>& dt, Time t) const = 0;

#ifdef MESSAGE_SUPPORT
    //! Export reference trajectory to message
    virtual void toMessage(corbo::messages::ReferenceTrajectory& message) const {}
    //! Import reference trajectory from message
    virtual void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) {}
#endif
};

using ReferenceTrajectoryFactory = Factory<ReferenceTrajectoryInterface>;
#define FACTORY_REGISTER_REFERENCE_TRAJECTORY(type) FACTORY_REGISTER_OBJECT(type, ReferenceTrajectoryInterface)

/**
 * @brief Static reference trajectory
 *
 * @ingroup core
 *
 * Stores a static reference trajectory containg a single set-point
 * with arbitrary components.
 *
 * @see ReferenceTrajectoryInterface ZeroReference
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class StaticReference : public ReferenceTrajectoryInterface
{
 public:
    using Ptr = std::shared_ptr<StaticReference>;

    StaticReference() {}
    explicit StaticReference(const Eigen::Ref<const OutputVector>& ref) : _ref(ref) {}

    ReferenceTrajectoryInterface::Ptr getInstance() const override { return std::make_shared<StaticReference>(); }

    bool isStatic() const override { return true; }
    bool isZero() const override { return _ref.isZero(); }
    int getDimension() const override { return (int)_ref.rows(); }

    void precompute(double /*dt*/, int /*n*/, Time /*t*/) override {}
    void precompute(const std::vector<double>& /*dt*/, Time /*t*/) override {}

    void getReference(const Time& /*t*/, OutputVector& ref) const override { ref = _ref; }

    const OutputVector& getReferenceCached(int /*k*/) const override { return _ref; }

    const OutputVector& getNextSteadyState(const Time& t) override { return _ref; }

    void setReference(const Eigen::Ref<const OutputVector>& ref) { _ref = ref; }

    bool isCached(double /*dt*/, int /*n*/, Time /*t*/) const override { return true; }
    bool isCached(const std::vector<double>& /*dt*/, Time /*t*/) const override { return true; }

// import / export support
#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::ReferenceTrajectory& message) const override;
    void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    OutputVector _ref;
};
FACTORY_REGISTER_REFERENCE_TRAJECTORY(StaticReference)

/**
 * @brief Zero reference trajectory
 *
 * @ingroup core
 *
 * Stores a zero reference vector \f$ [0, 0, \dotsc, 0]^T \f$ according to a desired dimension.
 *
 * @see ReferenceTrajectoryInterface StaticReference
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class ZeroReference : public StaticReference
{
 public:
    using Ptr = std::shared_ptr<ZeroReference>;

    ZeroReference() {}
    explicit ZeroReference(int dimension) : StaticReference(Eigen::VectorXd::Zero(dimension)) {}

    ReferenceTrajectoryInterface::Ptr getInstance() const override { return std::make_shared<ZeroReference>(); }

    bool isZero() const override { return true; }

    void setDimension(int dimension) { setReference(Eigen::VectorXd::Zero(dimension)); }

#ifdef MESSAGE_SUPPORT
    // import / export support
    void toMessage(corbo::messages::ReferenceTrajectory& message) const override;
    void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) override;
#endif
};
FACTORY_REGISTER_REFERENCE_TRAJECTORY(ZeroReference)

/**
 * @brief Sine reference trajectory
 *
 * @ingroup core
 *
 * Stores a non-static reference trajectory containg a sinusoide
 * function of the form f= amplitude * sin(ometa * t + offset)
 *
 * @author Rodrigo Velasco (rodrigo.velasco@tu-dortmund.de)
 */
class SineReferenceTrajectory : public ReferenceTrajectoryInterface
{
 public:
    using Ptr = std::shared_ptr<SineReferenceTrajectory>;

    SineReferenceTrajectory() {}
    explicit SineReferenceTrajectory(const double& amplitude, const double& omega, const double offset)
        : _amplitude(amplitude), _omega(omega), _offset(offset)
    {
        // Initialize the cached trajectory with the FIRST trajectory point
        _cached_trajectory.resize(1);
        Time t(0);
        getReference(t, _cached_trajectory[0]);
    }

    ReferenceTrajectoryInterface::Ptr getInstance() const override { return std::make_shared<SineReferenceTrajectory>(); }

    bool isStatic() const override { return false; }
    bool isZero() const override { return false; }
    int getDimension() const override { return 1; }

    bool isCached(double dt, int n, Time t) const override
    {
        if (_cached_dt.empty() || dt != _cached_dt[0] || n != _cached_trajectory.size() || t != _cached_t) return false;
        return true;
    }

    bool isCached(const std::vector<double>& dt, Time t) const override
    {
        if (_cached_dt.empty() || dt.size() != _cached_dt.size() || t != _cached_t) return false;
        for (int i = 0; i < dt.size(); ++i)
        {
            if (std::abs(dt[i] - _cached_dt[i]) < 1e-15) return false;  // TODO(roesmann): floating point comparison?
        }
        return true;
    }

    void precompute(double dt, int n, Time t) override
    {
        _cached_trajectory.resize(n);

        Duration d;

        for (int i = 0; i < n; ++i)
        {
            d.fromSec(dt * double(i));
            getReference(t + d, _cached_trajectory[i]);
        }
        _cached_dt.resize(1);
        _cached_dt[0] = dt;
        _cached_t     = t;
    }

    void precompute(const std::vector<double>& dt, Time t) override
    {
        _cached_trajectory.resize(dt.size() + 1);

        Duration d;

        getReference(t, _cached_trajectory[0]);
        for (int i = 0; i < dt.size(); ++i)
        {
            d.fromSec(dt[i]);
            getReference(t + d, _cached_trajectory[i + 1]);
        }
        _cached_dt = dt;
        _cached_t  = t;
    }

    void getReference(const Time& t, OutputVector& ref) const override
    {
        ref.resize(1);
        ref(0) = _amplitude * std::sin(_omega * t.toSec() + _offset);
    }

    const OutputVector& getReferenceCached(int k) const override
    {
        if (k >= _cached_trajectory.size())
        {
            PRINT_ERROR("SineReferenceTrajectory::getReferenceCached: k is not a valid index for cached reference. Returning zero value");
            return _zero_vector;
        }

        return _cached_trajectory.at(k);
    }

    const OutputVector& getNextSteadyState(const Time& t) override
    {
        PRINT_ERROR_ONCE("SineReferenceTrajectory: No steady state in periodic reference. Returning zero value.");

        return _zero_vector;
    }

    void setParameters(const double& amplitude, const double& omega, const double offset)
    {
        _amplitude = amplitude;
        _omega     = omega;
        _offset    = offset;
    }

 private:
    double _amplitude = 1;
    double _omega     = 1;
    double _offset    = 0;
    std::vector<OutputVector> _cached_trajectory;
    std::vector<double> _cached_dt;
    Time _cached_t;
    OutputVector _zero_vector = Eigen::VectorXd::Zero(1);

#ifdef MESSAGE_SUPPORT
    // import / export support
    void toMessage(corbo::messages::ReferenceTrajectory& message) const override;
    void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) override;
#endif
};
FACTORY_REGISTER_REFERENCE_TRAJECTORY(SineReferenceTrajectory)

/**
 * @brief discrete time reference trajectory
 *
 * @ingroup core
 *
 * Stores a non-static reference trajectory containg a discrete
 * time trajectory.
 *
 * @author Rodrigo Velasco (rodrigo.velasco@tu-dortmund.de)
 */
class DiscreteTimeReferenceTrajectory : public ReferenceTrajectoryInterface
{
 public:
    using Ptr = std::shared_ptr<DiscreteTimeReferenceTrajectory>;

    DiscreteTimeReferenceTrajectory() : _trajectory(std::make_shared<TimeSeries>()) {}
    explicit DiscreteTimeReferenceTrajectory(TimeSeries::Ptr trajectory, TimeSeries::Interpolation interpolation) : _interpolation(interpolation)
    {
        setTrajectory(trajectory);
    }

    ReferenceTrajectoryInterface::Ptr getInstance() const override { return std::make_shared<DiscreteTimeReferenceTrajectory>(); }

    bool isStatic() const override { return _trajectory && _trajectory->getTimeDimension() == 1; }
    bool isZero() const override
    {
        if (!_trajectory || _trajectory->getValueDimension() == 0) return false;

        return std::all_of(_trajectory->getValues().begin(), _trajectory->getValues().end(), [](double i) { return (i < 1e-9 && i > -1e-9); });
    }
    int getDimension() const override { return _trajectory ? _trajectory->getValueDimension() : 0; }

    bool isCached(double dt, int n, Time t) const override
    {
        if (_cached_dt.empty() || dt != _cached_dt[0] || n != _cached_trajectory.size() || t != _cached_t) return false;
        return true;
    }

    bool isCached(const std::vector<double>& dt, Time t) const override
    {
        if (_cached_dt.empty() || dt.size() != _cached_dt.size() || t != _cached_t) return false;
        for (int i = 0; i < dt.size(); ++i)
        {
            if (std::abs(dt[i] - _cached_dt[i]) < 1e-15) return false;  // TODO(roesmann): floating point comparison?
        }
        return true;
    }

    void precompute(double dt, int n, Time t) override
    {
        _cached_trajectory.resize(n);

        Duration d;

        for (int i = 0; i < n; ++i)
        {
            d.fromSec(dt * double(i));
            getReference(t + d, _cached_trajectory[i]);
        }
        _cached_dt.resize(1);
        _cached_dt[0] = dt;
        _cached_t     = t;
    }

    void precompute(const std::vector<double>& dt, Time t) override
    {
        _cached_trajectory.resize(dt.size() + 1);

        Duration d;

        getReference(t, _cached_trajectory[0]);
        for (int i = 0; i < dt.size(); ++i)
        {
            d.fromSec(dt[i]);
            getReference(t + d, _cached_trajectory[i + 1]);
        }
        _cached_dt = dt;
        _cached_t  = t;
    }

    void getReference(const Time& t, OutputVector& ref) const override
    {
        if (!_trajectory || _trajectory->getValueDimension() == 0)
        {
            PRINT_ERROR("DiscreteTimeReferenceTrajectory: trajectory is empty.");
            return;
        }

        double time = t.toSec() - _trajectory->getTimeFromStart();

        ref.resize(_trajectory->getValueDimension());

        if (time <= 0 || _trajectory->getTimeDimension() ==
                             1)  // TODO(roesmann): should we really check for time <= 0? do we also need to take time_from_start into account?
            ref = _trajectory->getValuesMap(0);
        else if (time >= _trajectory->getFinalTime())
            ref = _trajectory->getValuesMap(_trajectory->getTimeDimension() - 1);
        else
            _trajectory->getValuesInterpolate(time, ref, _interpolation, TimeSeries::Extrapolation::ZeroOrderHold);
    }

    const OutputVector& getReferenceCached(int k) const override
    {
        if (k >= _cached_trajectory.size())
        {
            PRINT_ERROR(
                "DiscreteTimeReferenceTrajectory::getReferenceCached: k is not a valid index for cached reference. Returning next steady state");
            return _next_steady_state;
        }

        return _cached_trajectory[k];
    }

    const OutputVector& getNextSteadyState(const Time& t) override { return _next_steady_state; }

    void setTrajectory(TimeSeries::Ptr trajectory, TimeSeries::Interpolation interpolation)
    {
        setTrajectory(trajectory);
        _interpolation = interpolation;
    }

    void setTrajectory(TimeSeries::Ptr trajectory)
    {
        if (!trajectory || trajectory->getValueDimension() == 0)
        {
            PRINT_ERROR("DiscreteTimeReferenceTrajectory must be initialized with at least one point in trajectory");
            return;
        }

        _trajectory        = trajectory;
        _next_steady_state = trajectory->getValuesMap(trajectory->getTimeDimension() - 1);

        // Initialize the cached trajectory with the FIRST trajectory point
        _cached_trajectory.resize(1);
        _cached_trajectory[0] = trajectory->getValuesMap(0);
    }

    void setTimeFromStart(const Time& time_from_start)
    {
        if (_trajectory) _trajectory->setTimeFromStart(time_from_start.toSec());
    }

    void setInterpolationMethod(TimeSeries::Interpolation interpolation) { _interpolation = interpolation; }

    const TimeSeries::Ptr& getReferenceTrajectory() { return _trajectory; }

#ifdef MESSAGE_SUPPORT
    // import / export support
    void toMessage(corbo::messages::ReferenceTrajectory& message) const override;
    void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    TimeSeries::Ptr _trajectory;
    std::vector<OutputVector> _cached_trajectory;
    std::vector<double> _cached_dt;
    Time _cached_t;
    OutputVector _next_steady_state;
    TimeSeries::Interpolation _interpolation = TimeSeries::Interpolation::ZeroOrderHold;
};
FACTORY_REGISTER_REFERENCE_TRAJECTORY(DiscreteTimeReferenceTrajectory)

/**
 * @brief discrete time reference trajectory in which the precompute step is not aware of the trajectory
 *
 * @ingroup core
 *
 * Stores a non-static reference trajectory containg a discrete
 * time trajectory. The precompute step assumes that the current value is constant in the future..
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class BlindDiscreteTimeReferenceTrajectory : public ReferenceTrajectoryInterface
{
 public:
    using Ptr = std::shared_ptr<BlindDiscreteTimeReferenceTrajectory>;

    BlindDiscreteTimeReferenceTrajectory() : _trajectory(std::make_shared<TimeSeries>()) {}
    explicit BlindDiscreteTimeReferenceTrajectory(TimeSeries::Ptr trajectory, TimeSeries::Interpolation interpolation) : _interpolation(interpolation)
    {
        setTrajectory(trajectory);
    }

    ReferenceTrajectoryInterface::Ptr getInstance() const override { return std::make_shared<BlindDiscreteTimeReferenceTrajectory>(); }

    bool isStatic() const override { return _trajectory && _trajectory->getTimeDimension() == 1; }
    bool isZero() const override
    {
        if (!_trajectory || _trajectory->getValueDimension() == 0) return false;

        return std::all_of(_trajectory->getValues().begin(), _trajectory->getValues().end(), [](double i) { return (i < 1e-9 && i > -1e-9); });
    }
    int getDimension() const override { return _trajectory ? _trajectory->getValueDimension() : 0; }

    bool isCached(double dt, int n, Time t) const override
    {
        if (_output_cache.size() != getDimension() || t != _cached_t) return false;
        return true;
    }

    bool isCached(const std::vector<double>& dt, Time t) const override
    {
        if (_output_cache.size() != getDimension() || t != _cached_t) return false;
        return true;
    }

    void precompute(double dt, int n, Time t) override
    {
        getReference(t, _output_cache);
        _cached_t = t;
    }

    void precompute(const std::vector<double>& dt, Time t) override
    {
        getReference(t, _output_cache);
        _cached_t = t;
    }

    void getReference(const Time& t, OutputVector& ref) const override
    {
        if (!_trajectory || _trajectory->getValueDimension() == 0)
        {
            PRINT_ERROR("BlindDiscreteTimeReferenceTrajectory: trajectory is empty.");
            return;
        }

        double time = t.toSec() - _trajectory->getTimeFromStart();

        ref.resize(_trajectory->getValueDimension());

        if (time <= 0 || _trajectory->getTimeDimension() == 1)
            ref = _trajectory->getValuesMap(0);
        else if (time >= _trajectory->getFinalTime())
            ref = _trajectory->getValuesMap(_trajectory->getTimeDimension() - 1);
        else
            _trajectory->getValuesInterpolate(time, ref, _interpolation, TimeSeries::Extrapolation::ZeroOrderHold);
    }

    const OutputVector& getReferenceCached(int k) const override { return _output_cache; }

    const OutputVector& getNextSteadyState(const Time& t) override { return _next_steady_state; }

    void setTrajectory(TimeSeries::Ptr trajectory, TimeSeries::Interpolation interpolation)
    {
        setTrajectory(trajectory);
        _interpolation = interpolation;
    }

    void setTrajectory(TimeSeries::Ptr trajectory)
    {
        if (!trajectory || trajectory->getValueDimension() == 0)
        {
            PRINT_ERROR("BlindDiscreteTimeReferenceTrajectory must be initialized with at least one point in trajectory");
            return;
        }

        _trajectory        = trajectory;
        _next_steady_state = trajectory->getValuesMap(trajectory->getTimeDimension() - 1);
    }

    void setTimeFromStart(const Time& time_from_start)
    {
        if (_trajectory) _trajectory->setTimeFromStart(time_from_start.toSec());
    }

    void setInterpolationMethod(TimeSeries::Interpolation interpolation) { _interpolation = interpolation; }

    const TimeSeries::Ptr& getReferenceTrajectory() { return _trajectory; }

#ifdef MESSAGE_SUPPORT
    // import / export support
    void toMessage(corbo::messages::ReferenceTrajectory& message) const override;
    void fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    TimeSeries::Ptr _trajectory;
    OutputVector _output_cache;
    Time _cached_t;
    OutputVector _next_steady_state;
    TimeSeries::Interpolation _interpolation = TimeSeries::Interpolation::ZeroOrderHold;
};
FACTORY_REGISTER_REFERENCE_TRAJECTORY(BlindDiscreteTimeReferenceTrajectory)

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_REFERENCE_TRAJECTORY_H_

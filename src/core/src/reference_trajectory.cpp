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

#include <corbo-core/reference_trajectory.h>

namespace corbo {

#ifdef MESSAGE_SUPPORT
void StaticReference::toMessage(corbo::messages::ReferenceTrajectory& message) const
{
    message.mutable_static_reference()->mutable_xf()->Resize(_ref.rows(), 0);
    Eigen::Map<OutputVector>(message.mutable_static_reference()->mutable_xf()->mutable_data(), _ref.rows()) = _ref;
}

void StaticReference::fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues)
{
    int dim = message.static_reference().xf_size();
    if (dim >= 1)
    {
        _ref.resize(dim);
        for (int i = 0; i < dim; ++i) _ref[i] = message.static_reference().xf(i);
    }
    else if (issues)
        *issues << "StaticReference: dimension must be larger than zero.\n";
}
#endif

#ifdef MESSAGE_SUPPORT
void ZeroReference::toMessage(corbo::messages::ReferenceTrajectory& message) const
{
    message.mutable_zero_reference()->set_dimension(getDimension());
}

void ZeroReference::fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues)
{
    int dim = message.zero_reference().dimension();
    if (dim >= 1)
    {
        setDimension(dim);
    }
    else if (issues)
        *issues << "ZeroReference: dimension must be larger than zero.\n";
}
#endif

#ifdef MESSAGE_SUPPORT
void SineReferenceTrajectory::toMessage(corbo::messages::ReferenceTrajectory& message) const
{
    message.mutable_sine_reference_trajectory()->set_amplitude(_amplitude);
    message.mutable_sine_reference_trajectory()->set_omega(_omega);
    message.mutable_sine_reference_trajectory()->set_offset(_offset);
}

void SineReferenceTrajectory::fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues)
{
    setParameters(message.sine_reference_trajectory().amplitude(), message.sine_reference_trajectory().omega(),
                  message.sine_reference_trajectory().offset());
}
#endif

#ifdef MESSAGE_SUPPORT
void DiscreteTimeReferenceTrajectory::toMessage(corbo::messages::ReferenceTrajectory& message) const
{
    if (_trajectory) _trajectory->toMessage(*message.mutable_discrete_time_reference()->mutable_time_series());
    switch (_interpolation)
    {
        case TimeSeries::Interpolation::ZeroOrderHold:
        {
            message.mutable_discrete_time_reference()->set_interpolation(messages::DiscreteTimeReferenceTrajectory_Interpolation_ZERO_ORDER);
            break;
        }
        case TimeSeries::Interpolation::Linear:
        {
            message.mutable_discrete_time_reference()->set_interpolation(messages::DiscreteTimeReferenceTrajectory_Interpolation_LINEAR);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Unknown interpolation method.");
        }
    }
}

void DiscreteTimeReferenceTrajectory::fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues)
{
    TimeSeries::Ptr time_series = std::make_shared<TimeSeries>();
    time_series->fromMessage(message.discrete_time_reference().time_series());
    setTrajectory(time_series);

    switch (message.discrete_time_reference().interpolation())
    {
        case messages::DiscreteTimeReferenceTrajectory_Interpolation_ZERO_ORDER:
        {
            setInterpolationMethod(TimeSeries::Interpolation::ZeroOrderHold);
            break;
        }
        case messages::DiscreteTimeReferenceTrajectory_Interpolation_LINEAR:
        {
            setInterpolationMethod(TimeSeries::Interpolation::Linear);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Unknown interpolation method.");
        }
    }
}
#endif

#ifdef MESSAGE_SUPPORT
void BlindDiscreteTimeReferenceTrajectory::toMessage(corbo::messages::ReferenceTrajectory& message) const
{
    if (_trajectory) _trajectory->toMessage(*message.mutable_blind_discrete_time_reference()->mutable_time_series());

    switch (_interpolation)
    {
        case TimeSeries::Interpolation::ZeroOrderHold:
        {
            message.mutable_blind_discrete_time_reference()->set_interpolation(
                messages::BlindDiscreteTimeReferenceTrajectory_Interpolation_ZERO_ORDER);
            break;
        }
        case TimeSeries::Interpolation::Linear:
        {
            message.mutable_blind_discrete_time_reference()->set_interpolation(messages::BlindDiscreteTimeReferenceTrajectory_Interpolation_LINEAR);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Unknown interpolation method.");
        }
    }
}

void BlindDiscreteTimeReferenceTrajectory::fromMessage(const corbo::messages::ReferenceTrajectory& message, std::stringstream* issues)
{
    TimeSeries::Ptr time_series = std::make_shared<TimeSeries>();
    time_series->fromMessage(message.blind_discrete_time_reference().time_series());
    setTrajectory(time_series);

    switch (message.blind_discrete_time_reference().interpolation())
    {
        case messages::BlindDiscreteTimeReferenceTrajectory_Interpolation_ZERO_ORDER:
        {
            setInterpolationMethod(TimeSeries::Interpolation::ZeroOrderHold);
            break;
        }
        case messages::BlindDiscreteTimeReferenceTrajectory_Interpolation_LINEAR:
        {
            setInterpolationMethod(TimeSeries::Interpolation::Linear);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Unknown interpolation method.");
        }
    }
}
#endif

}  // namespace corbo

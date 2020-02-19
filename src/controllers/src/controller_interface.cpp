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

#include <corbo-controllers/controller_interface.h>

namespace corbo {

bool ControllerInterface::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                               const Duration& dt, const Time& t, ControllerInterface::ControlVector& u, SignalTargetInterface* signal_target,
                               const std::string& ns)
{
    TimeSeries::Ptr u_sequence;
    TimeSeries::Ptr x_sequence;

    if (step(x, xref, uref, dt, t, u_sequence, x_sequence, signal_target, nullptr, nullptr, nullptr, ns))
    {
        u_sequence->getValuesInterpolate(0.0, u);
        return true;
    }
    else
    {
        return false;
    }
}
}  // namespace corbo

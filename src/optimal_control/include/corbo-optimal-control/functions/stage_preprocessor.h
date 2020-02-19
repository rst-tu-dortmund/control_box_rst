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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_STAGE_PREPROCESSOR_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_STAGE_PREPROCESSOR_H_

#include <corbo-core/factory.h>
#include <corbo-core/reference_trajectory.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/optimal_control/stage_preprocessors.pb.h>
#endif

#include <Eigen/Core>

namespace corbo {

class DiscretizationGridInterface;

class StagePreprocessor
{
 public:
    using Ptr      = std::shared_ptr<StagePreprocessor>;
    using ConstPtr = std::shared_ptr<const StagePreprocessor>;

    StagePreprocessor();

    virtual Ptr getInstance() const                                         = 0;
    virtual void precompute(const Eigen::Ref<const Eigen::VectorXd>& input) = 0;
    virtual bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                        bool single_dt, const Eigen::VectorXd& x0, const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/)
    {
        return true;
    }

    std::vector<Eigen::VectorXd> _vector_data;
    std::vector<Eigen::MatrixXd> _matrix_data;

#ifdef MESSAGE_SUPPORT
    virtual bool fromMessage(const messages::StagePreprocessors& message, std::stringstream* issues) { return true; }
    virtual void toMessage(messages::StagePreprocessors& message) const {}
#endif
};

using StagePreprocessorFactory = Factory<StagePreprocessor>;
#define FACTORY_REGISTER_STAGE_PREPROCESSOR(type) FACTORY_REGISTER_OBJECT(type, StagePreprocessor)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_STAGE_PREPROCESSOR_H_

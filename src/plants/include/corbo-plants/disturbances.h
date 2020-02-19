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

#ifndef SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCES_H_
#define SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCES_H_

#include <corbo-plants/disturbance_interface.h>

#include <memory>
#include <random>

namespace corbo {

class DisturbanceGaussianNoise : public DisturbanceInterface
{
 public:
    using Ptr = std::shared_ptr<DisturbanceGaussianNoise>;

    //! Virtual destructor
    virtual ~DisturbanceGaussianNoise() {}

    // implements interface method
    DisturbanceInterface::Ptr getInstance() const override { return std::make_shared<DisturbanceGaussianNoise>(); }

    // implements interface method
    void disturb(const Time& t, const Eigen::Ref<const Eigen::VectorXd>& values, Eigen::Ref<Eigen::VectorXd> disturbed_values) override;

    void setParameters(const Eigen::Ref<const Eigen::VectorXd>& mean_vec, const Eigen::Ref<const Eigen::VectorXd>& std_vec)
    {
        assert(mean_vec.size() == std_vec.size());
        _mean = mean_vec;
        _std  = std_vec;
    }

    void setSeed(int seed);

    bool checkParameters(int values_dim, std::stringstream* issues) const override;

    // implements interface method
    void reset() override;

    void initializeDistributions();

#ifdef MESSAGE_SUPPORT
    //! Export plant settings to message
    virtual void toMessage(messages::DisturbanceGaussianNoise& message) const;
    //! Import plant settings from message
    virtual void fromMessage(const messages::DisturbanceGaussianNoise& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Disturbance& message) const override { toMessage(*message.mutable_gaussian_noise()); }
    // implements interface method
    void fromMessage(const messages::Disturbance& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.gaussian_noise(), issues);
    }
#endif

 private:
    Eigen::VectorXd _mean;
    Eigen::VectorXd _std;

    int _seed                   = 1;
    std::mt19937 _random_engine = std::mt19937(_seed);
    // TODO(roesmann) should we use multiple engines (1 per distribution) to avoid possible correlatons?

    std::vector<std::normal_distribution<double>> _distributions;
};

FACTORY_REGISTER_DISTURBANCE(DisturbanceGaussianNoise)

}  // namespace corbo

#endif  // SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCES_H_

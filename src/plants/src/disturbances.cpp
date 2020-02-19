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

#include <corbo-plants/disturbances.h>

#include <corbo-core/value_comparison.h>

namespace corbo {

void DisturbanceGaussianNoise::disturb(const Time& t, const Eigen::Ref<const Eigen::VectorXd>& values, Eigen::Ref<Eigen::VectorXd> disturbed_values)
{
    // first copy values
    disturbed_values = values;

    if (values.size() != _mean.size() || values.size() != _std.size())
    {
        PRINT_ERROR_NAMED("Cannot disturb due to dimension mismatch");
        return;
    }

    // check if we already have already initialized our distributions
    if (_distributions.size() != values.size()) initializeDistributions();

    // now add noise
    for (int i = 0; i < values.size(); ++i)
    {
        if (approx_zero(_std[i]))
        {
            if (!approx_zero(_mean[i]))
                disturbed_values[i] += _mean[i];
            else
                continue;  // we skip this component
        }
        disturbed_values[i] += _distributions[i](_random_engine);
    }
}

void DisturbanceGaussianNoise::initializeDistributions()
{
    reset();

    if (_mean.size() != _std.size())
    {
        PRINT_ERROR_NAMED("Cannot initialize distributions due to dimension mismatch between mean and std vector.");
        return;
    }

    _distributions.clear();
    for (int i = 0; i < _mean.size(); ++i)
    {
        _distributions.emplace_back(_mean[i], _std[i]);
    }
}

void DisturbanceGaussianNoise::setSeed(int seed)
{
    _seed = seed;
    reset();  // this should also reset the seed
}

void DisturbanceGaussianNoise::reset()
{
    _random_engine.seed(_seed);
    _distributions.clear();
}

bool DisturbanceGaussianNoise::checkParameters(int values_dim, std::stringstream* issues) const
{
    bool success = true;

    if (_mean.size() != values_dim)
    {
        if (issues)
            *issues << "DisturbanceGaussianNoise: Mean vector dimension (" << _mean.size() << ") does not match required value dimension ("
                    << values_dim << ")." << std::endl;
        success = false;
    }

    if (_std.size() != values_dim)
    {
        if (issues)
            *issues << "DisturbanceGaussianNoise: Std vector dimension (" << _std.size() << ") does not match required value dimension ("
                    << values_dim << ")." << std::endl;
        success = false;
    }
    return success;
}

#ifdef MESSAGE_SUPPORT
void DisturbanceGaussianNoise::toMessage(corbo::messages::DisturbanceGaussianNoise& message) const
{
    // mean
    message.mutable_mean_vec()->Resize(_mean.rows(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_mean_vec()->mutable_data(), _mean.rows()) = _mean;

    // std
    message.mutable_std_vec()->Resize(_std.rows(), 0);
    Eigen::Map<Eigen::VectorXd>(message.mutable_std_vec()->mutable_data(), _std.rows()) = _std;
}

void DisturbanceGaussianNoise::fromMessage(const corbo::messages::DisturbanceGaussianNoise& message, std::stringstream* issues)
{
    // mean
    if (message.mean_vec_size() > 0)
        _mean = Eigen::Map<const Eigen::VectorXd>(message.mean_vec().data(), message.mean_vec_size());
    else
        _mean.resize(0);

    // std
    if (message.std_vec_size() > 0)
        _std = Eigen::Map<const Eigen::VectorXd>(message.std_vec().data(), message.mean_vec_size());
    else
        _std.resize(0);

    initializeDistributions();
}
#endif

}  // namespace corbo

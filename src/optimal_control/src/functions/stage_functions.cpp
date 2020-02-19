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

#include <corbo-optimal-control/functions/stage_functions.h>

namespace corbo {

bool StageFunction::hasIntegralTerms(int k) const { return getIntegralStateControlTermDimension(k) > 0; }
bool StageFunction::hasNonIntegralTerms(int k) const
{
    return getNonIntegralControlDeviationTermDimension(k) > 0 || getNonIntegralControlTermDimension(k) > 0 || getNonIntegralDtTermDimension(k) > 0 ||
           getNonIntegralStateControlDtTermDimension(k) > 0 || getNonIntegralStateControlTermDimension(k) > 0 ||
           getNonIntegralStateTermDimension(k) > 0;
}

int StageFunction::getConcatenatedNonIntegralStateTermDimension(int k, bool lsq_mode) const
{
    if (!lsq_mode)
        return getNonIntegralStateTermDimension(k) + getNonIntegralStateControlTermDimension(k) + getNonIntegralStateControlDtTermDimension(k);

    int dim = 0;
    dim += isLsqFormNonIntegralStateTerm(k) ? 1 : getNonIntegralStateTermDimension(k);

    return dim + getNonIntegralStateControlTermDimension(k) + getNonIntegralStateControlDtTermDimension(k);
}

int StageFunction::getConcatenatedNonIntegralStateControlTermDimension(int k, bool lsq_mode) const
{
    // excluding control deviation
    if (!lsq_mode)
        return getNonIntegralStateTermDimension(k) + getNonIntegralControlTermDimension(k) + getNonIntegralStateControlTermDimension(k) +
               getNonIntegralStateControlDtTermDimension(k);

    int dim = 0;
    dim += isLsqFormNonIntegralStateTerm(k) ? 1 : getNonIntegralStateTermDimension(k);
    dim += isLsqFormNonIntegralControlTerm(k) ? 1 : getNonIntegralControlTermDimension(k);

    return dim + getNonIntegralStateControlTermDimension(k) + getNonIntegralStateControlDtTermDimension(k);
}

void StageFunction::computeConcatenatedNonIntegralStateTerms(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                             const Eigen::Ref<const Eigen::VectorXd>& u_k, double dt_k,
                                                             Eigen::Ref<Eigen::VectorXd> cost, bool lsq_mode) const
{
    assert(cost.size() == getConcatenatedNonIntegralStateTermDimension(k));
    int idx = 0;
    if (getNonIntegralStateTermDimension(k) > 0)
    {
        if (lsq_mode)
        {
            Eigen::VectorXd temp(getNonIntegralStateTermDimension(k));
            computeNonIntegralStateTerm(k, x_k, temp);
            cost[0] = temp.squaredNorm();
            ++idx;
        }
        else
        {
            computeNonIntegralStateTerm(k, x_k, cost.head(getNonIntegralStateTermDimension(k)));
            idx += getNonIntegralStateTermDimension(k);
        }
    }
    if (getNonIntegralStateControlTermDimension(k) > 0)
    {
        computeNonIntegralStateControlTerm(k, x_k, u_k, cost.segment(idx, getNonIntegralStateControlTermDimension(k)));
        idx += getNonIntegralStateControlTermDimension(k);
    }
    if (getNonIntegralStateControlDtTermDimension(k) > 0)
    {
        computeNonIntegralStateControlDtTerm(k, x_k, u_k, dt_k, cost.segment(idx, getNonIntegralStateControlDtTermDimension(k)));
        idx += getNonIntegralStateControlDtTermDimension(k);
    }
    assert(idx == getConcatenatedNonIntegralStateTermDimension(k));
}

void StageFunction::computeConcatenatedNonIntegralStateControlTerms(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                                    const Eigen::Ref<const Eigen::VectorXd>& u_k, double dt_k,
                                                                    Eigen::Ref<Eigen::VectorXd> cost, bool lsq_mode) const
{
    assert(cost.size() == getConcatenatedNonIntegralStateControlTermDimension(k));
    int idx = 0;
    if (getNonIntegralStateTermDimension(k) > 0)
    {
        if (lsq_mode)
        {
            Eigen::VectorXd temp(getNonIntegralStateTermDimension(k));
            computeNonIntegralStateTerm(k, x_k, temp);
            cost[idx++] = temp.squaredNorm();
        }
        else
        {
            computeNonIntegralStateTerm(k, x_k, cost.head(getNonIntegralStateTermDimension(k)));
            idx += getNonIntegralStateTermDimension(k);
        }
    }
    if (getNonIntegralControlTermDimension(k) > 0)
    {
        if (lsq_mode)
        {
            Eigen::VectorXd temp(getNonIntegralControlTermDimension(k));
            computeNonIntegralControlTerm(k, u_k, temp);
            cost[idx++] = temp.squaredNorm();
        }
        else
        {
            computeNonIntegralControlTerm(k, u_k, cost.head(getNonIntegralControlTermDimension(k)));
            idx += getNonIntegralControlTermDimension(k);
        }
    }
    if (getNonIntegralStateControlTermDimension(k) > 0)
    {
        computeNonIntegralStateControlTerm(k, x_k, u_k, cost.segment(idx, getNonIntegralStateControlTermDimension(k)));
        idx += getNonIntegralStateControlTermDimension(k);
    }
    if (getNonIntegralStateControlDtTermDimension(k) > 0)
    {
        computeNonIntegralStateControlDtTerm(k, x_k, u_k, dt_k, cost.segment(idx, getNonIntegralStateControlDtTermDimension(k)));
        idx += getNonIntegralStateControlDtTermDimension(k);
    }
    assert(idx == getConcatenatedNonIntegralStateControlTermDimension(k));
}

}  // namespace corbo

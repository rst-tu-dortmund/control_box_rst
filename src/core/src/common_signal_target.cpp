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

#include <corbo-core/common_signal_target.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace corbo {

void CommonSignalTarget::sendMeasurement(Measurement::ConstPtr measurement)
{
    if (measurement->header.name.empty())
    {
        PRINT_ERROR_NAMED("signal name missing in signal header. Skipping signal.");
        return;
    }
    if (measurement->getValues().empty())
    {
        PRINT_WARNING_NAMED("measurement does not contain any values.");
    }

    std::string signal_name;
    SignalInterface::Ptr signal;
    SignalGroup* group = expandGroups(measurement->header.name, &signal_name, &signal);
    if (!signal)
    {
        // new signal, store all measurements with same name into a TimeSeries signal
        signal         = std::make_shared<TimeSeriesSignal>();
        signal->header = measurement->header;  // copy header to parent object
        group->group_signals.push_back(signal);
    }
    TimeSeriesSignal::Ptr ts = std::dynamic_pointer_cast<TimeSeriesSignal>(signal);
    if (!ts)
    {
        PRINT_ERROR_NAMED("previous signal with same name but different type detected. Skipping signal.");
        return;
    }
    ts->add(*measurement);
}

void CommonSignalTarget::sendTimeSeries(TimeSeriesSignal::Ptr time_series)
{
    if (time_series->header.name.empty())
    {
        PRINT_ERROR_NAMED("signal name missing in signal header. Skipping signal.");
        return;
    }

    std::string signal_name;
    SignalInterface::Ptr signal;
    SignalGroup* group = expandGroups(time_series->header.name, &signal_name, &signal);
    if (!signal)
    {
        // new signal, store all measurements with same name into a TimeSeries signal
        signal         = std::make_shared<TimeSeriesSequenceSignal>();
        signal->header = time_series->header;  // copy header to parent object
        group->group_signals.push_back(signal);
    }
    TimeSeriesSequenceSignal::Ptr sequence = std::dynamic_pointer_cast<TimeSeriesSequenceSignal>(signal);
    if (!sequence)
    {
        PRINT_ERROR_NAMED("previous signal with same name but different type detected. Skipping signal.");
        return;
    }
    sequence->add(time_series->getTimeSeriesPtr());
}

void CommonSignalTarget::sendIndexedValues(IndexedValuesSignal::Ptr indexed_values)
{
    if (indexed_values->header.name.empty())
    {
        PRINT_ERROR_NAMED("signal name missing in signal header. Skipping signal.");
        return;
    }
    if (indexed_values->getValueDimension() == 0)
    {
        PRINT_WARNING_NAMED("value vector is empty.");
    }

    std::string signal_name;
    SignalInterface::Ptr signal;
    SignalGroup* group = expandGroups(indexed_values->header.name, &signal_name, &signal);
    if (!signal)
    {
        // new signal, store all measurements with same name into a IndexedValuesSet signal
        signal         = std::make_shared<IndexedValuesSetSignal>();
        signal->header = indexed_values->header;  // copy header to parent object
        group->group_signals.push_back(signal);
    }
    IndexedValuesSetSignal::Ptr set = std::dynamic_pointer_cast<IndexedValuesSetSignal>(signal);
    if (!set)
    {
        PRINT_ERROR_NAMED("previous signal with same name but different type detected. Skipping signal.");
        return;
    }
    set->add(*indexed_values);
}

void CommonSignalTarget::sendIndexedValuesSet(IndexedValuesSetSignal::Ptr indexed_values_set)
{
    if (indexed_values_set->header.name.empty())
    {
        PRINT_ERROR_NAMED("signal name missing in signal header. Skipping signal.");
        return;
    }

    std::string signal_name;
    SignalInterface::Ptr signal;
    SignalGroup* group = expandGroups(indexed_values_set->header.name, &signal_name, &signal);
    if (!signal)
    {
        // new signal
        group->group_signals.push_back(indexed_values_set);
    }
    else
    {
        PRINT_ERROR_NAMED("Already received a previous signal with the same name. IndexedValuesSetSignal does not support signal stacking.");
    }
}

void CommonSignalTarget::sendMatrix(MatrixSignal::Ptr matrix)
{
    if (matrix->header.name.empty())
    {
        PRINT_ERROR_NAMED("signal name missing in signal header. Skipping signal.");
        return;
    }
    if (matrix->isEmpty())
    {
        PRINT_WARNING_NAMED("signal data is empty.");
    }

    std::string signal_name;
    SignalInterface::Ptr signal;
    SignalGroup* group = expandGroups(matrix->header.name, &signal_name, &signal);
    if (!signal)
    {
        // new signal, store all measurements with same name into a IndexedValuesSet signal
        signal         = std::make_shared<MatrixSetSignal>();
        signal->header = matrix->header;  // copy header to parent object
        group->group_signals.push_back(signal);
    }
    MatrixSetSignal::Ptr set = std::dynamic_pointer_cast<MatrixSetSignal>(signal);
    if (!set)
    {
        PRINT_ERROR_NAMED("previous signal with same name but different type detected. Skipping signal.");
        return;
    }
    set->add(matrix);
}

void CommonSignalTarget::removeSignal(const std::string& name)
{
    SignalGroup* signal_group = getGroup(name);

    if (!signal_group) return;
    if (signal_group == &_signals) return;

    SignalGroup* signal_group_parent = signal_group->parent;

    for (int i = 0; i < (int)signal_group_parent->children.size(); ++i)
    {
        if (&signal_group_parent->children[i] == signal_group)
        {
            signal_group_parent->children.erase(signal_group_parent->children.begin() + i);
            return;
        }
    }

    // TODO(Kraemer): Remove single parents

    PRINT_DEBUG_NAMED("Signal " << name << " not found");
}

CommonSignalTarget::SignalGroup* CommonSignalTarget::getGroup(const std::string& full_group_name)
{
    // parse groups / namespaces:
    std::vector<std::string> groups;
    std::istringstream stream(full_group_name);
    std::string current;
    while (std::getline(stream, current, SIGNAL_NAMESPACE_DELIMITER))
    {
        groups.push_back(current);
    }

    SignalGroup* group = &_signals;

    for (int i = 0; i < (int)groups.size(); ++i)
    {
        //        if (i == (int)groups.size() - 1)
        //        {
        //            return group;
        //        }
        auto it = std::find_if(group->children.begin(), group->children.end(), [&](const SignalGroup& item) { return item.group_name == groups[i]; });
        if (it != group->children.end())
        {
            group = &(*it);
        }
        else
        {
            // PRINT_DEBUG_NAMED("Signal " << full_group_name << " not found");
            return nullptr;
        }
    }
    return group;
}

SignalInterface::Ptr CommonSignalTarget::getSignal(const std::string& full_signal_name)
{
    // remove signal name from group name;
    std::string full_group_name;
    // std::string signal_name;
    std::size_t found = full_signal_name.find_last_of(SIGNAL_NAMESPACE_DELIMITER);
    if (found < full_signal_name.size())
    {
        full_group_name = full_signal_name.substr(0, found);
        // signal_name = full_signal_name.substr(found+1);
    }
    // else
    // {
    // signal_name = full_signal_name;
    // }

    // get leaf group
    SignalGroup* group = getGroup(full_group_name);
    if (group)
    {
        // search signals in group
        for (SignalInterface::Ptr& signal : group->group_signals)
        {
            if (signal->header.name.compare(full_signal_name) == 0)
            {
                // found
                return signal;
            }
        }
    }
    return {};
}

CommonSignalTarget::SignalGroup* CommonSignalTarget::expandGroups(const std::string full_signal_name, std::string* signal_name_out,
                                                                  SignalInterface::Ptr* signal)
{
    // parse groups / namespaces:
    std::vector<std::string> groups;
    std::istringstream stream(full_signal_name);
    std::string current;
    while (std::getline(stream, current, SIGNAL_NAMESPACE_DELIMITER))
    {
        groups.push_back(current);
    }

    // the last string in groups contains the actual signal name
    // now traverse tree and construct children if necessary
    int n              = (int)groups.size();
    SignalGroup* group = &_signals;
    for (int i = 0; i < n; ++i)
    {
        if (i == n - 1)  // final group reached, this item is the signal name
        {
            if (signal_name_out) *signal_name_out = groups.back();
            if (signal)
            {
                // also export signal if found
                auto it =
                    std::find_if(group->group_signals.begin(), group->group_signals.end(),
                                 [&full_signal_name](const SignalInterface::Ptr& sig) { return sig->header.name.compare(full_signal_name) == 0; });
                if (it != group->group_signals.end())
                    *signal = *it;
                else
                    signal->reset();
            }
            return group;
        }
        std::string group_name = groups[i];
        auto it                = std::find_if(group->children.begin(), group->children.end(),
                               [&group_name](const SignalGroup& child) { return child.group_name.compare(group_name) == 0; });
        if (it == group->children.end())
        {
            // group not created yet
            SignalGroup new_group;
            new_group.group_name = group_name;
            new_group.parent     = group;
            group->children.push_back(new_group);
            group = &group->children.back();
        }
        else
        {
            // group found
            group = &(*it);
        }
    }
    return group;
}

void CommonSignalTarget::clear()
{
    _signals.group_signals.clear();
    _signals.children.clear();
}

}  // namespace corbo

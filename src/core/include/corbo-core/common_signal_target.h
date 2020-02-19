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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_COMMON_SIGNAL_TARGET_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_COMMON_SIGNAL_TARGET_H_

#include <corbo-core/signal_target_interface.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace corbo {

/**
 * @brief General signal target to store incoming signals for post-processing
 *
 * @ingroup core signals
 *
 * Signals are stored into an internal tree structure according to their
 * namespaces/groups (refer to SignalInterface for a description of signal namespaces).
 * These signals are then available for any post-processing
 * purposes or file export.
 *
 * @see SignalTargetInterface SignalTargetRPC YamlExporter
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Add possibility to reserve space for data (e.g. for TimeSeries object)
 */
class CommonSignalTarget : public SignalTargetInterface
{
 public:
    //! Tree node containing a group of signals and a set of child tree nodes
    struct SignalGroup
    {
        std::string group_name;                           //!< Identifier for the current group / signal namespace
        std::vector<SignalInterface::Ptr> group_signals;  //!< set of signals in the current group
        std::vector<SignalGroup> children;                //!< Container for children (sub-namespace)
        SignalGroup* parent = nullptr;
    };

    using Ptr = std::shared_ptr<CommonSignalTarget>;

    //!< Virtual destructor
    virtual ~CommonSignalTarget() {}

    // implements interface method
    void registerMeasurement(const std::string& /*unique_name*/, int /*value_dimension*/, const std::vector<std::string>& /*value_labels*/ = {},
                             bool zero_order_hold = false) override
    {
    }

    // implements interface method
    void registerTimeSeries(const std::string& /*unique_name*/, int /*value_dimension*/, bool zero_order_hold = false) override {}

    // implements interface method
    void sendMeasurement(Measurement::ConstPtr measurement) override;

    // implements interface method
    void sendTimeSeries(TimeSeriesSignal::Ptr time_series) override;

    // implements interface method
    void sendIndexedValues(IndexedValuesSignal::Ptr indexed_values) override;

    // implements interface method
    void sendIndexedValuesSet(IndexedValuesSetSignal::Ptr indexed_values_set) override;

    // implements interface method
    void sendMatrix(MatrixSignal::Ptr matrix) override;

    void setTopLevelGroupName(const std::string& name) { _signals.group_name = name; }

    // const SignalGroup* getSignalGroup(const std::string full_signal_name) const { return resolveGroupLevel(full_signal_name); }

    void removeSignal(const std::string& name);

    SignalGroup* getGroup(const std::string& full_group_name);
    SignalInterface::Ptr getSignal(const std::string& full_signal_name);

    bool isEmpty() { return _signals.children.empty(); }

    /**
     * @brief Retrieve read-only-access to the underlying signal-tree
     * @returns Signal tree structure starting with the root node of type CommonSignalTarget::SignalGroup
     */
    const SignalGroup& getSignals() const { return _signals; }

    SignalGroup& getSignalsRef() { return _signals; }

    //! Erase stored signals
    void clear();

 protected:
    /**
     * @brief Parse full signal name, create tree-sub-groups if not found and return related signal group.
     *
     * This methods parses the provded signal name and expands the tree according to the signal namespaces/groups.
     * The signal is not added to its corresponding group in the tree yet, but a pointer to its lowest containing group
     * is returned.
     *
     * @remark Make sure to check whether a signal with the same identifier is already present (refer to the third argument).
     *
     * @param[in]  full_signal_name    Signal name containing all namespaces (separated by '/').
     * @param[out] signal_name_out     Signal name without preceding namespaces/groups [optional]
     * @param[out] signal              Retrieve a previous signal with the same identifier if any, otherwise return null [optional].
     *
     * @returns Pointer to the lowest containing group/namespace of the signal (which might be newly created during expansion).
     */
    SignalGroup* expandGroups(const std::string full_signal_name, std::string* signal_name_out = nullptr, SignalInterface::Ptr* signal = nullptr);

 private:
    SignalGroup _signals;  //!< Root element of the signal tree
};

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_COMMON_SIGNAL_TARGET_H_

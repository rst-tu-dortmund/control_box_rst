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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_DATA_EXPORTER_INTERFACE_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_DATA_EXPORTER_INTERFACE_H_

#include <corbo-core/common_signal_target.h>

#include <corbo-core/factory.h>

#include <memory>
#include <string>

namespace corbo {

/**
 * @brief Interface class for exporting signals
 *
 * @ingroup core data_export
 *
 * @author Maximilian Krämer (maximilian.kraemer@tu-dortmund.de)
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class DataExporterInterface
{
 public:
    using Ptr      = std::shared_ptr<DataExporterInterface>;
    using ConstPtr = std::shared_ptr<const DataExporterInterface>;

    virtual ~DataExporterInterface() {}

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    //! Get access to the associated factory
    static Factory<DataExporterInterface>& getFactory() { return Factory<DataExporterInterface>::instance(); }

    virtual std::string getFormatName() const = 0;
    virtual std::string getFileSuffix() const = 0;

    virtual bool isSupportingSignalGroup() const { return false; }
    virtual bool isSupportingTimeSeriesSignal() const { return false; }
    virtual bool isSupportingTimeSeries() const { return false; }
    virtual bool isSupportingTimeSeriesSequenceSignal() const { return false; }
    virtual bool isSupportingIndexedValuesSignal() const { return false; }
    virtual bool isSupportingIndexedValuesSetSignal() const { return false; }
    virtual bool isSupportingMatrixSignal() const { return false; }
    virtual bool isSupportingMatrixSetSignal() const { return false; }

    virtual bool exportSignalGroup(const std::string& filename, const CommonSignalTarget::SignalGroup& signal_group);
    virtual bool exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal);
    virtual bool exportTimeSeries(const std::string& filename, const TimeSeries& time_series);
    virtual bool exportTimeSeriesSequenceSignal(const std::string& filename, const TimeSeriesSequenceSignal& signal);
    virtual bool exportIndexedValuesSignal(const std::string& filename, const IndexedValuesSignal& signal);
    virtual bool exportIndexedValuesSetSignal(const std::string& filename, const IndexedValuesSetSignal& signal);
    virtual bool exportMatrixSignal(const std::string& filename, const MatrixSignal& signal);
    virtual bool exportMatrixSetSignal(const std::string& filename, const MatrixSetSignal& signal);
};

using DataExporterFactory = Factory<DataExporterInterface>;
#define FACTORY_REGISTER_DATA_EXPORTER(type) FACTORY_REGISTER_OBJECT(type, DataExporterInterface)

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_DATA_EXPORTER_INTERFACE_H_

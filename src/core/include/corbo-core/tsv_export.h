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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_TSV_EXPORT_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_TSV_EXPORT_H_

#include <corbo-core/data_exporter_interface.h>

#include <fstream>
#include <memory>

namespace corbo {

/**
 * @brief Export Tab-separated values (tsv) file
 *
 * @ingroup core data_export
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class TsvExporter : public DataExporterInterface
{
 public:
    using Ptr = std::shared_ptr<TsvExporter>;
    
    TsvExporter() = default;

    DataExporterInterface::Ptr getInstance() const override { return std::make_shared<TsvExporter>(); }
    
    std::string getFormatName() const override { return "Tsv"; }
    std::string getFileSuffix() const override { return ".tsv"; }
    
    bool isSupportingSignalGroup() const override { return false; }
    bool isSupportingTimeSeriesSignal() const override { return true; }
    bool isSupportingTimeSeries() const override { return true; }
    bool isSupportingIndexedValuesSignal() const override { return false; }
    bool isSupportingIndexedValuesSetSignal() const override { return false; }

    bool exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal) override;
    bool exportTimeSeries(const std::string& filename, const TimeSeries& time_series) override;
};
FACTORY_REGISTER_DATA_EXPORTER(TsvExporter);

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_TSV_EXPORT_H_

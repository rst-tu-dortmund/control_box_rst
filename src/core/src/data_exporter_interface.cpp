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

#include <corbo-core/data_exporter_interface.h>

#include <corbo-core/console.h>

namespace corbo {

bool DataExporterInterface::exportSignalGroup(const std::string& filename, const CommonSignalTarget::SignalGroup& signal_group)
{
    PRINT_ERROR_NAMED("SignalGroup is not supported for export.");
    return false;
}

bool DataExporterInterface::exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal)
{
    PRINT_ERROR_NAMED("TimeSeriesSignal is not supported for export.");
    return false;
}

bool DataExporterInterface::exportTimeSeries(const std::string& filename, const TimeSeries& time_series)
{
    PRINT_ERROR_NAMED("TimeSeries is not supported for export.");
    return false;
}

bool DataExporterInterface::exportTimeSeriesSequenceSignal(const std::string& filename, const TimeSeriesSequenceSignal& signal)
{
    PRINT_ERROR_NAMED("TimeSeriesSequenceSignal is not supported for export.");
    return false;
}

bool DataExporterInterface::exportIndexedValuesSignal(const std::string& filename, const IndexedValuesSignal& signal)
{
    PRINT_ERROR_NAMED("IndexedValuesSignal is not supported for export.");
    return false;
}

bool DataExporterInterface::exportIndexedValuesSetSignal(const std::string& filename, const IndexedValuesSetSignal& signal)
{
    PRINT_ERROR_NAMED("IndexedValuesSetSignal is not supported for export.");
    return false;
}

bool DataExporterInterface::exportMatrixSignal(const std::string& filename, const MatrixSignal& signal)
{
    PRINT_ERROR_NAMED("MatrixSignal is not supported for export.");
    return false;
}

bool DataExporterInterface::exportMatrixSetSignal(const std::string& filename, const MatrixSetSignal& signal)
{
    PRINT_ERROR_NAMED("MatrixSetSignal is not supported for export.");
    return false;
}

}  // namespace corbo

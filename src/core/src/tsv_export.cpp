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

#include <corbo-core/tsv_export.h>

#include <fstream>

namespace corbo {

bool TsvExporter::exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal)
{
    if (!signal.getTimeSeries()) return false;
    return exportTimeSeries(filename, *signal.getTimeSeries());
}

bool TsvExporter::exportTimeSeries(const std::string& filename, const TimeSeries& time_series)
{
    std::ofstream file;
    file.open(filename);
    if (!file.is_open()) return false;
    
    if (time_series.getTimeDimension() != time_series.getValuesMatrixView().cols())
    {
        PRINT_ERROR_NAMED("time dimension does not match number of columns in values matrix");
        return false;
    }
    
    // add field names (first line)
    file << "t";
    for (int i = 0; i < time_series.getValueDimension(); ++i)
    {
        // use value label if present
        if (i < time_series.getValueLabels().size())
            file << "\t" << time_series.getValueLabels()[i];
        else
            file << "\tx" << i;
    }
    file << "\n";  // add new line. Note, std:endl also flushes the stream
    
    // now add values
    for (int i = 0; i < time_series.getTimeDimension(); ++i)
    {
        // add time
        file << time_series.getTime()[i];
        // add values
        for (int j = 0; j < time_series.getValueDimension(); ++j)
        {
            file << "\t" << time_series.getValuesMap(i)[j];
        }
        // start new line
        file << "\n";
    }
        
    file.close();
    return true;
}


}  // namespace corbo


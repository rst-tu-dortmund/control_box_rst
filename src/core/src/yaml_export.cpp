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

#ifdef YAML_SUPPORT

#include <corbo-core/yaml_export.h>

namespace corbo {

YamlExporter::YamlExporter() {}

bool YamlExporter::exportSignalGroup(const std::string& filename, const CommonSignalTarget::SignalGroup& signal_group)
{
    YAML::Emitter yaml_emitter;

    emitSignalGroup(signal_group, yaml_emitter);

    // write yaml emitter stream to given file
    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitTimeSeries(signal, yaml_emitter);

    // write yaml emitter stream to given file
    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportTimeSeries(const std::string& filename, const TimeSeries& time_series)
{
    YAML::Emitter yaml_emitter;

    emitTimeSeries(time_series, yaml_emitter);

    // write yaml emitter stream to given file
    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportTimeSeriesSequenceSignal(const std::string& filename, const TimeSeriesSequenceSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitTimeSeriesSequence(signal, yaml_emitter);

    // write yaml emitter stream to given file
    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportIndexedValuesSignal(const std::string& filename, const IndexedValuesSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitIndexedValues(signal, yaml_emitter);

    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportIndexedValuesSetSignal(const std::string& filename, const IndexedValuesSetSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitIndexedValuesSet(signal, yaml_emitter);

    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportMatrixSignal(const std::string& filename, const MatrixSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitMatrix(signal, yaml_emitter);

    return write_file(yaml_emitter, filename);
}

bool YamlExporter::exportMatrixSetSignal(const std::string& filename, const MatrixSetSignal& signal)
{
    YAML::Emitter yaml_emitter;

    emitMatrixSet(signal, yaml_emitter);

    return write_file(yaml_emitter, filename);
}

void YamlExporter::emitSignalGroup(const CommonSignalTarget::SignalGroup& signal_group, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal_group.group_name;
    emitter << YAML::Value;

    emitter << YAML::BeginMap;

    if (!signal_group.group_signals.empty())
    {
        // emitter << YAML::Key << "signals";
        // emitter << YAML::Value << YAML::BeginMap;

        for (const SignalInterface::Ptr& signal : signal_group.group_signals)
        {
            // check we have supported types
            const TimeSeriesSignal* ts_signal = dynamic_cast<const TimeSeriesSignal*>(signal.get());
            if (ts_signal)
            {
                emitTimeSeries(*ts_signal, emitter, false);
                continue;
            }
            const TimeSeriesSequenceSignal* ts_seq_signal = dynamic_cast<const TimeSeriesSequenceSignal*>(signal.get());
            if (ts_seq_signal)
            {
                emitTimeSeriesSequence(*ts_seq_signal, emitter, false);
                continue;
            }
            const IndexedValuesSignal* indexed_values_signal = dynamic_cast<const IndexedValuesSignal*>(signal.get());
            if (indexed_values_signal)
            {
                emitIndexedValues(*indexed_values_signal, emitter, false);
                continue;
            }
            const IndexedValuesSetSignal* indexed_values_set_signal = dynamic_cast<const IndexedValuesSetSignal*>(signal.get());
            if (indexed_values_set_signal)
            {
                emitIndexedValuesSet(*indexed_values_set_signal, emitter, false);
                continue;
            }
            const MatrixSignal* matrix_signal = dynamic_cast<const MatrixSignal*>(signal.get());
            if (matrix_signal)
            {
                emitMatrix(*matrix_signal, emitter, false);
                continue;
            }
            const MatrixSetSignal* matrix_set_signal = dynamic_cast<const MatrixSetSignal*>(signal.get());
            if (matrix_set_signal)
            {
                emitMatrixSet(*matrix_set_signal, emitter, false);
                continue;
            }
            PRINT_WARNING_NAMED("Signal '" << signal->header.name << "' cannot be exported (not implemented).");
            // emitter << "unsupported_signal";
        }

        // emitter << YAML::EndMap;
    }

    if (!signal_group.children.empty())
    {
        // iterate children
        for (const CommonSignalTarget::SignalGroup& child : signal_group.children)
        {
            // start recursion
            emitSignalGroup(child, emitter, false);
        }
    }

    emitter << YAML::EndMap;

    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitTimeSeries(const TimeSeries& time_series, YAML::Emitter& emitter, bool wrap_in_map)
{
    // go through time series signal and fill out yaml emitter
    if (wrap_in_map) emitter << YAML::BeginMap;

    std::vector<std::string> labels = time_series.getValueLabels();
    emitter << YAML::Key << "label";
    emitter << YAML::Value << YAML::Flow << labels;

    emitter << YAML::Key << "time";
    emitter << YAML::Value << YAML::Flow << time_series.getTime();

    emitter << YAML::Key << "rows";
    emitter << YAML::Value << time_series.getValuesMatrixView().rows();

    emitter << YAML::Key << "cols";
    emitter << YAML::Value << time_series.getValuesMatrixView().cols();

    emitter << YAML::Key << "is_column_major";
    // emitter << YAML::Value << YAML::TrueFalseBool << !time_series.getValuesMatrixView().IsRowMajor;
    emitter << YAML::Value << (int)time_series.isColumnMajor();

    emitter << YAML::Key << "values";
    emitter << YAML::Value << YAML::Flow << time_series.getValues();

    emitter << YAML::Key << "time_from_start";
    emitter << YAML::Value << time_series.getTimeFromStart();

    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitTimeSeries(const TimeSeriesSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    if (signal.getTimeSeries())
    {
        const TimeSeries& time_series   = *signal.getTimeSeries();
        std::vector<std::string> labels = time_series.getValueLabels();
        emitter << YAML::Key << "label";
        emitter << YAML::Value << YAML::Flow << labels;

        emitter << YAML::Key << "time";
        emitter << YAML::Value << YAML::Flow << time_series.getTime();

        emitter << YAML::Key << "rows";
        emitter << YAML::Value << time_series.getValuesMatrixView().rows();

        emitter << YAML::Key << "cols";
        emitter << YAML::Value << time_series.getValuesMatrixView().cols();

        emitter << YAML::Key << "is_column_major";
        // emitter << YAML::Value << YAML::TrueFalseBool << !time_series.getValuesMatrixView().IsRowMajor;
        emitter << YAML::Value << (int)!time_series.getValuesMatrixView().IsRowMajor;

        emitter << YAML::Key << "values";
        emitter << YAML::Value << YAML::Flow << time_series.getValues();

        emitter << YAML::Key << "time_from_start";
        emitter << YAML::Value << time_series.getTimeFromStart();
    }
    emitter << YAML::EndMap;
    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitTimeSeriesSequence(const TimeSeriesSequenceSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    if (!signal.isEmpty())
    {
        emitter << YAML::Key << "ts_seq";
        emitter << YAML::Value << YAML::BeginSeq;
        for (const TimeSeries::Ptr& ts : signal.getSequence()->getSequence())
        {
            emitTimeSeries(*ts, emitter, true);
        }
        emitter << YAML::EndSeq;
    }

    emitter << YAML::EndMap;

    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitIndexedValues(const IndexedValuesSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    emitter << YAML::Key << "values";
    emitter << YAML::Value << YAML::Flow << YAML::BeginSeq;

    // start tuple
    emitter << YAML::BeginSeq;

    // insert index
    emitter << signal.getIndex();

    // insert data
    emitter << signal.getValues();

    // end tupel
    emitter << YAML::EndSeq;

    emitter << YAML::EndSeq;
    emitter << YAML::EndMap;
    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitIndexedValuesSet(const IndexedValuesSetSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    // emitter << YAML::Comment("You are using the free version of corbo. The free version is limited to 1 export per hour.");
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    emitter << YAML::Key << "values";
    emitter << YAML::Value << YAML::BeginSeq;

    const std::map<int, std::vector<double>>& map = signal.getData();
    auto it = map.begin();

    while (it != map.end())
    {
        // start tupel
        emitter << YAML::Flow << YAML::BeginSeq;

        // insert index
        emitter << it->first;

        // insert data
        emitter << it->second;
        // emitter << it->second.size();

        // end tupel
        emitter << YAML::EndSeq;
        it++;
    }

    emitter << YAML::EndSeq;
    emitter << YAML::EndMap;
    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitMatrix(const MatrixSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    if (!signal.isEmpty())
    {
        emitter << YAML::Key << "label";
        emitter << YAML::Value << signal.getLabel();

        emitter << YAML::Key << "rows";
        emitter << YAML::Value << signal.getRowDimension();

        emitter << YAML::Key << "cols";
        emitter << YAML::Value << signal.getColDimension();

        emitter << YAML::Key << "is_column_major";
        // emitter << YAML::Value << YAML::TrueFalseBool << !signal.getMatrix().IsRowMajor;
        emitter << YAML::Value << 1;  // see for-loops below

        emitter << YAML::Key << "values";
        emitter << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (int col_idx = 0; col_idx < signal.getMatrix().cols(); ++col_idx)
            for (int row_idx = 0; row_idx < signal.getMatrix().rows(); ++row_idx) emitter << signal.getMatrix()(row_idx, col_idx);

        emitter << YAML::EndSeq;
    }

    emitter << YAML::EndMap;
    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitMatrixSet(const MatrixSetSignal& signal, YAML::Emitter& emitter, bool wrap_in_map)
{
    if (wrap_in_map) emitter << YAML::BeginMap;

    emitter << YAML::Key << signal.header.getShortName();
    emitter << YAML::Value << YAML::BeginMap;

    emitHeader(signal.header, emitter);

    if (!signal.isEmpty())
    {
        emitter << YAML::Key << "mat_seq";
        emitter << YAML::Value << YAML::BeginSeq;
        for (const MatrixSignal::Ptr& mat : signal.getData())
        {
            emitMatrix(*mat, emitter, true);
        }
        emitter << YAML::EndSeq;
    }

    emitter << YAML::EndMap;

    if (wrap_in_map) emitter << YAML::EndMap;
}

void YamlExporter::emitHeader(const SignalHeader& header, YAML::Emitter& emitter)
{
    emitter << YAML::Key << "header";
    emitter << YAML::Value << YAML::BeginMap;

    emitter << YAML::Key << "name";
    emitter << YAML::Value << header.name;

    emitter << YAML::Key << "time_stamp";
    emitter << YAML::Value << header.time.toSec();

    emitter << YAML::EndMap;
}

bool YamlExporter::write_file(const YAML::Emitter& emitter, const std::string& filename)
{
    std::ofstream file;
    file.open(filename);
    if (file.is_open())
    {
        file << emitter.c_str();
        file.close();
        return true;
    }
    else
    {
        return false;
    }
}

}  // namespace corbo

#endif  // YAML_SUPPORT

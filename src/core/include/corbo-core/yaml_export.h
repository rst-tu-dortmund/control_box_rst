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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_YAML_EXPORT_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_YAML_EXPORT_H_

#include <corbo-core/data_exporter_interface.h>
#include <yaml-cpp/yaml.h>

#include <fstream>
#include <memory>

namespace corbo {

/**
 * @brief Export yaml files
 *
 * @ingroup core data_export
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 * @author Maximilian Krämer (maximilian.kraemer@tu-dortmund.de)
 */
class YamlExporter : public DataExporterInterface
{
 public:
    using Ptr = std::shared_ptr<YamlExporter>;

    YamlExporter();

    DataExporterInterface::Ptr getInstance() const override { return std::make_shared<YamlExporter>(); }

    std::string getFormatName() const override { return "Yaml"; }
    std::string getFileSuffix() const override { return ".yaml"; }

    bool isSupportingSignalGroup() const override { return true; }
    bool isSupportingTimeSeriesSignal() const override { return true; }
    bool isSupportingTimeSeries() const override { return true; }
    bool isSupportingTimeSeriesSequenceSignal() const override { return true; }
    bool isSupportingIndexedValuesSignal() const override { return true; }
    bool isSupportingIndexedValuesSetSignal() const override { return true; }
    bool isSupportingMatrixSignal() const override { return true; }
    bool isSupportingMatrixSetSignal() const override { return true; }

    bool exportSignalGroup(const std::string& filename, const CommonSignalTarget::SignalGroup& signal_group) override;
    bool exportTimeSeriesSignal(const std::string& filename, const TimeSeriesSignal& signal) override;
    bool exportTimeSeries(const std::string& filename, const TimeSeries& time_series) override;
    bool exportTimeSeriesSequenceSignal(const std::string& filename, const TimeSeriesSequenceSignal& signal) override;
    bool exportIndexedValuesSignal(const std::string& filename, const IndexedValuesSignal& signal) override;
    bool exportIndexedValuesSetSignal(const std::string& filename, const IndexedValuesSetSignal& signal) override;
    bool exportMatrixSignal(const std::string& filename, const MatrixSignal& signal) override;
    bool exportMatrixSetSignal(const std::string& filename, const MatrixSetSignal& signal) override;

 private:
    void emitSignalGroup(const CommonSignalTarget::SignalGroup& signal_group, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitTimeSeries(const TimeSeries& time_series, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitTimeSeries(const TimeSeriesSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitTimeSeriesSequence(const TimeSeriesSequenceSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitIndexedValues(const IndexedValuesSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitIndexedValuesSet(const IndexedValuesSetSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitMatrix(const MatrixSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);
    void emitMatrixSet(const MatrixSetSignal& signal, YAML::Emitter& emitter, bool wrap_in_map = true);

    void emitHeader(const SignalHeader& header, YAML::Emitter& emitter);

    bool write_file(const YAML::Emitter& emitter, const std::string& filename);
};
FACTORY_REGISTER_DATA_EXPORTER(YamlExporter);

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_YAML_EXPORT_H_

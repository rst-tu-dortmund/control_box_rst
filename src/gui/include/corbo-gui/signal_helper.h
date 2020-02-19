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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_HELPER_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_HELPER_H_

#include <corbo-core/common_signal_target.h>
#include <corbo-core/signals.h>
#include <QHash>
#include <QObject>
#include <QSet>
#include <QString>
#include <map>
#include <memory>

namespace corbo {
namespace gui {

class SignalHelper : public QObject
{
    Q_OBJECT

 public:
    using Ptr      = std::shared_ptr<SignalHelper>;
    using ConstPtr = std::shared_ptr<const SignalHelper>;

    static constexpr const int ALL_VALUES = -1;

    struct SignalData
    {
        QString short_name;
        QString name;
        QStringList namespaces;
        int dimension        = 1;
        bool zero_order_hold = false;
        SignalInterface::Ptr signal;
        int task_id = 0;
        QSet<int> active_indices;
        QStringList value_labels;
    };
    using SignalMap = QHash<QString, SignalData>;

    explicit SignalHelper(QObject* parent = 0) : QObject(parent) {}

    void addSignal(const messages::Signal& signal);
    void addMeasurement(const QString& name, Measurement::ConstPtr measurement);
    void addTimeSeries(const QString& name, TimeSeriesSignal::Ptr ts);
    void addIndexedValues(const QString& name, IndexedValuesSignal::Ptr indexed_values);
    void addIndexedValuesSet(const QString& name, IndexedValuesSetSignal::Ptr indexed_values_set);
    void addMatrixSignal(const QString& name, MatrixSignal::Ptr matrix_signal);
    // void addSignal(const QString& name, const SignalInterface::Ptr signal);

    void removeSignal(const QString& key);
    void removeSignal(const QString& key, int value_idx);
    void clearCurrentSeries() { clearSeries(_id); }
    void clearSeries(int task_id);
    void clearNamespace(const QString& namespace_pattern, int task_id);
    void clearAll();

    SignalData* signalData(const QString& key);
    const SignalData* signalData(const QString& key) const;

    const SignalMap& signalMap() const { return _signal_map; }
    SignalMap& signalMap() { return _signal_map; }

    static void extractNamespace(const QString& name, QStringList& namespaces);
    static QString name2Key(const QString& name, int id);
    static bool key2Name(const QString& key, QString& name, int& id);

    int currentSeriesId() const { return _id; }

    void startNewSeries()
    {
        ++_id;
        emit newSeries(_id);
    }

    const std::map<int, CommonSignalTarget>& getSignalTree() const { return _signal_tree_map; }
    std::map<int, CommonSignalTarget>& getSignalTreeRef() { return _signal_tree_map; }

    void printSignals();

 signals:
    void newSignal(const QString& key, SignalData& signal_data);  // TODO(roesmann) do we actually need write access to the signal data?

    void newMeasurement(const QString& key, Measurement::ConstPtr measurement, SignalData& signal_data, bool first);
    void newTimeSeries(const QString& key, TimeSeriesSignal::ConstPtr ts, SignalData& signal_data, bool first);
    void newIndexedValues(const QString& key, IndexedValuesSignal::ConstPtr indexed_values, SignalData& signal_data, bool first);
    // void newIndexedValuesSet(const QString& key, IndexedValuesSetSignal::ConstPtr indexed_values, SignalData& signal_data, bool first);
    void newMatrixSignal(const QString& key, MatrixSignal::ConstPtr matrix_signal, SignalData& signal_data, bool first);
    void signalRemoved(const QString& key, int value_idx);
    void newSeries(int task_id);

 protected:
    QString name2Key(const QString& name) const;

 private:
    SignalMap _signal_map;
    std::map<int, CommonSignalTarget> _signal_tree_map;

    int _id = 0;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_HELPER_H_

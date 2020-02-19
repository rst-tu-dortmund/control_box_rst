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

#include <corbo-core/console.h>
#include <corbo-gui/signal_helper.h>
#include <corbo-gui/utilities.h>

#include <memory>
#include <string>

namespace corbo {
namespace gui {

void SignalHelper::addSignal(const messages::Signal& signal)
{
    std::stringstream issues;

    if (signal.header().name().empty())
    {
        PRINT_WARNING_NAMED("Received signal without id/namespace. Skipping...");
        return;
    }

    if (signal.has_measurement())
    {
        Measurement::Ptr measurement = std::make_shared<Measurement>();
        measurement->fromMessage(signal, &issues);
        if (!issues.str().empty())
        {
            PRINT_WARNING_NAMED(issues.str());
            return;
        }
        addMeasurement(QString::fromStdString(signal.header().name()), measurement);
    }
    else if (signal.has_time_series())
    {
        TimeSeriesSignal::Ptr ts = std::make_shared<TimeSeriesSignal>();
        ts->fromMessage(signal, &issues);

        if (!issues.str().empty())
        {
            PRINT_WARNING_NAMED(issues.str());
            return;
        }
        addTimeSeries(QString::fromStdString(signal.header().name()), ts);
    }
    else if (signal.has_indexed_values())
    {
        IndexedValuesSignal::Ptr indexed_values = std::make_shared<IndexedValuesSignal>();
        indexed_values->fromMessage(signal, &issues);

        if (!issues.str().empty())
        {
            PRINT_WARNING_NAMED(issues.str());
            return;
        }
        addIndexedValues(QString::fromStdString(signal.header().name()), indexed_values);
    }
    else if (signal.has_indexed_values_set())
    {
        IndexedValuesSetSignal::Ptr indexed_values_set = std::make_shared<IndexedValuesSetSignal>();
        indexed_values_set->fromMessage(signal, &issues);

        if (!issues.str().empty())
        {
            PRINT_WARNING_NAMED(issues.str());
            return;
        }
        addIndexedValuesSet(QString::fromStdString(signal.header().name()), indexed_values_set);
    }
    else if (signal.has_matrix())
    {
        MatrixSignal::Ptr matrix_signal = std::make_shared<MatrixSignal>();
        matrix_signal->fromMessage(signal, &issues);

        if (!issues.str().empty())
        {
            PRINT_WARNING_NAMED(issues.str());
            return;
        }
        addMatrixSignal(QString::fromStdString(signal.header().name()), matrix_signal);
    }
    else
    {
        PRINT_ERROR_NAMED("The gui currently does not support the message type of " << signal.header().name());
    }

    //    else if (signal.has_values())
    //    {
    //        ValuesSignal::Ptr signal_obj = std::make_shared<ValuesSignal>();
    //        signal_obj->fromMessage(signal, &issues);
    //        if (!issues.str().empty())
    //        {
    //            PRINT_WARNING("SignalHelper::addSignal(): " << issues.str());
    //            return;
    //        }
    //        addValuesSignal(QString::fromStdString(signal.header().name()), signal_obj);
    //    }
}

void SignalHelper::addMeasurement(const QString& name, Measurement::ConstPtr measurement)
{
    QString key = name2Key(name);
    auto it     = _signal_map.find(key);

    if (it == _signal_map.end())
    {
        // new signal
        it       = _signal_map.insert(key, SignalData());
        it->name = name;
        extractNamespace(name, it->namespaces);
        it->short_name      = it->namespaces.back();
        it->task_id         = _id;
        it->dimension       = measurement->header.value_dimension;
        it->zero_order_hold = measurement->header.zero_order_hold;
        for (const std::string& label : measurement->getValueLabels()) it->value_labels.push_back(QString::fromStdString(label));

        TimeSeriesSignal::Ptr ts = std::make_shared<TimeSeriesSignal>(measurement->header.value_dimension);
        ts->header               = measurement->header;
        if (!measurement->getValues().empty()) ts->add(*measurement);
        it->signal = ts;
        emit newMeasurement(key, measurement, *it, true);
        emit newSignal(key, *it);
    }
    else
    {
        if (!measurement->getValues().empty())
        {
            std::static_pointer_cast<TimeSeriesSignal>(it->signal)->add(*measurement);
        }
        emit newMeasurement(key, measurement, *it, false);
    }

    if (!measurement->getValues().empty())
    {
        _signal_tree_map[_id].sendMeasurement(measurement);
    }
}

void SignalHelper::addTimeSeries(const QString& name, TimeSeriesSignal::Ptr ts)
{
    QString key = name2Key(name);
    auto it     = _signal_map.find(key);
    if (it == _signal_map.end())
    {
        // new signal
        it       = _signal_map.insert(key, SignalData());
        it->name = name;
        extractNamespace(name, it->namespaces);
        it->short_name = it->namespaces.back();
        it->task_id    = _id;
        it->dimension =
            ts->getTimeSeries() ? ts->getTimeSeries()->getValueDimension() : ts->header.value_dimension;  // TODO(roesmann) unqiue variable?
        it->zero_order_hold = ts->header.zero_order_hold;
        std::vector<std::string> sublabels;
        ts->getValueLabels(sublabels);
        for (const std::string& label : sublabels) it->value_labels.push_back(QString::fromStdString(label));

        // Create parent signal
        TimeSeriesSequenceSignal::Ptr sequence = std::make_shared<TimeSeriesSequenceSignal>(it->dimension);
        sequence->header                       = ts->header;
        if (!ts->isEmpty()) sequence->add(ts->getTimeSeriesPtr());
        it->signal = sequence;
        emit newTimeSeries(key, ts, *it, true);
        emit newSignal(key, *it);
    }
    else
    {
        if (!ts->isEmpty())
        {
            std::static_pointer_cast<TimeSeriesSequenceSignal>(it->signal)->add(ts->getTimeSeriesPtr());
        }
        emit newTimeSeries(key, ts, *it, false);
    }

    if (!ts->isEmpty())
    {
        // add time series to signal tree map
        _signal_tree_map[_id].sendTimeSeries(ts);
    }
}

void SignalHelper::addIndexedValues(const QString& name, IndexedValuesSignal::Ptr indexed_values)
{
    QString key = name2Key(name);
    auto it     = _signal_map.find(key);

    if (it == _signal_map.end())
    {
        // new signal
        it       = _signal_map.insert(key, SignalData());
        it->name = name;
        extractNamespace(name, it->namespaces);
        it->short_name      = it->namespaces.back();
        it->task_id         = _id;
        it->dimension       = 1;  // values are not intended to be shown as separate signal-widgets.
        it->zero_order_hold = indexed_values->header.zero_order_hold;

        IndexedValuesSetSignal::Ptr set = std::make_shared<IndexedValuesSetSignal>();
        set->header                     = indexed_values->header;
        if (!indexed_values->isEmpty()) set->add(*indexed_values);
        it->signal = set;
        emit newIndexedValues(key, indexed_values, *it, true);
        emit newSignal(key, *it);
    }
    else
    {
        if (!indexed_values->isEmpty())
        {
            std::static_pointer_cast<IndexedValuesSetSignal>(it->signal)->add(*indexed_values);
        }
        emit newIndexedValues(key, indexed_values, *it, false);
    }

    if (!indexed_values->isEmpty())
    {
        _signal_tree_map[_id].sendIndexedValues(indexed_values);
    }
}

void SignalHelper::addIndexedValuesSet(const QString& name, IndexedValuesSetSignal::Ptr indexed_values_set)
{
    PRINT_ERROR_NAMED("not yet implemented!");
    /*
    QString key = name2Key(name);
    auto it     = _signal_map.find(key);
    if (it == _signal_map.end())
    {
        // new signal
        it       = _signal_map.insert(key, SignalData());
        it->name = name;
        extractNamespace(name, it->namespaces);
        it->short_name = it->namespaces.back();
        it->task_id    = _id;
        it->dimension =
        ts->getTimeSeries() ? ts->getTimeSeries()->getValueDimension() : ts->header.value_dimension;  // TODO(roesmann) unqiue variable?
        it->zero_order_hold = ts->header.zero_order_hold;
        std::vector<std::string> sublabels;
        ts->getValueLabels(sublabels);
        for (const std::string& label : sublabels) it->value_labels.push_back(QString::fromStdString(label));

        // Create parent signal
        TimeSeriesSequenceSignal::Ptr sequence = std::make_shared<TimeSeriesSequenceSignal>(it->dimension);
        sequence->header                       = ts->header;
        if (!ts->isEmpty()) sequence->add(ts->getTimeSeriesPtr());
        it->signal = sequence;
        emit newIndexedValuesSet(key, indexed_values_set, *it, true);
        emit newSignal(key, *it);
    }
    else
    {
        if (!indexed_values_set->isEmpty())
        {
            std::static_pointer_cast<TimeSeriesSequenceSignal>(it->signal)->add(ts->getTimeSeriesPtr());
        }
        emit newIndexedValuesSet(key, indexed_values_set, *it, false);
    }
     */
}

void SignalHelper::addMatrixSignal(const QString& name, MatrixSignal::Ptr matrix_signal)
{
    QString key = name2Key(name);
    auto it     = _signal_map.find(key);

    if (it == _signal_map.end())
    {
        // new signal
        it       = _signal_map.insert(key, SignalData());
        it->name = name;
        extractNamespace(name, it->namespaces);
        it->short_name      = it->namespaces.back();
        it->task_id         = _id;
        it->dimension       = 1;  // values are not intended to be shown as separate signal-widgets.
        it->zero_order_hold = matrix_signal->header.zero_order_hold;

        MatrixSetSignal::Ptr set = std::make_shared<MatrixSetSignal>();
        set->header              = matrix_signal->header;
        if (!matrix_signal->isEmpty()) set->add(matrix_signal);
        it->signal = set;
        emit newMatrixSignal(key, matrix_signal, *it, true);
        emit newSignal(key, *it);
    }
    else
    {
        if (!matrix_signal->isEmpty())
        {
            std::static_pointer_cast<MatrixSetSignal>(it->signal)->add(matrix_signal);
        }
        emit newMatrixSignal(key, matrix_signal, *it, false);
    }

    if (!matrix_signal->isEmpty())
    {
        _signal_tree_map[_id].sendMatrix(matrix_signal);
    }
}

// void SignalHelper::addSignal(const QString& name, const SignalInterface::Ptr signal)
//{

//    QString key = name2Key(name);
//    auto it     = _signal_map.find(key);

//    if (it != _signal_map.end())
//    {
//        return;   // TODO
//    }
//    //        // signal already known, replace
//    //        it->name = name;
//    //        extractNamespace(name, it->namespaces);
//    //        it->short_name = it->namespaces.back();
//    //        it->signal     = signal;
//    //        it->task_id    = _id;
//    //        emit newSignal(key, *it, true);
//    //    }

//    // check parent type
//    switch (signal->header.parent_type)
//    {
//        case SignalType::NoParent:
//        {
//            break;
//        }
//        default:
//        {
//        }
//    }

//    else
//    {

//        // check if we need an appropriate parent signal
//        if (signal->header.has_parent)
//        {
//            it->signal = signal;
//        }
//        else
//        {
//            // add signal directly
//        }
//        emit newSignal(key, *it, false);
//    }
//}

SignalHelper::SignalData* SignalHelper::signalData(const QString& key)
{
    auto signal_it = _signal_map.find(key);
    if (signal_it != _signal_map.end())
    {
        return &signal_it.value();
    }
    return nullptr;
}
const SignalHelper::SignalData* SignalHelper::signalData(const QString& key) const
{
    auto signal_it = _signal_map.find(key);
    if (signal_it != _signal_map.end())
    {
        return &signal_it.value();
    }
    return nullptr;
}

void SignalHelper::extractNamespace(const QString& name, QStringList& namespaces) { namespaces = name.split(corbo::SIGNAL_NAMESPACE_DELIMITER); }

QString SignalHelper::name2Key(const QString& name) const { return name2Key(name, _id); }
QString SignalHelper::name2Key(const QString& name, int id) { return QString::number(id) + util::SIGNAL_NAMESPACE_PREFIX_DELIMITER + name; }

bool SignalHelper::key2Name(const QString& key, QString& name, int& id)
{
    bool ret_val      = true;
    QStringList token = key.split(util::SIGNAL_NAMESPACE_PREFIX_DELIMITER);
    if (token.size() < 2)
    {
        name    = key;
        id      = 0;
        ret_val = false;
    }
    else
    {
        id   = token.front().toInt(&ret_val);
        name = token.back();
    }
    PRINT_WARNING_COND(!ret_val, "SignalHelper::key2Name(): unknown format in " << key.toStdString());
    return ret_val;
}

void SignalHelper::clearSeries(int task_id)
{
    QMutableHashIterator<QString, SignalData> map_it(_signal_map);
    while (map_it.hasNext())
    {
        map_it.next();

        if (map_it.value().task_id == task_id)
        {
            QString removed_name = map_it.value().name;
            int removed_id       = map_it.value().task_id;
            map_it.remove();
            emit signalRemoved(name2Key(removed_name, removed_id), ALL_VALUES);
        }
    }
    _signal_tree_map.erase(task_id);
}

void SignalHelper::clearNamespace(const QString& namespace_pattern, int task_id)
{
    QMutableHashIterator<QString, SignalData> map_it(_signal_map);
    while (map_it.hasNext())
    {
        map_it.next();

        if (task_id < 0 || map_it.value().task_id != task_id) continue;

        if (map_it.value().name.contains(namespace_pattern))
        {
            QString removed_name = map_it.value().name;
            int removed_id       = map_it.value().task_id;
            map_it.remove();
            emit signalRemoved(name2Key(removed_name, removed_id), ALL_VALUES);
        }
    }
    PRINT_WARNING_NAMED("This function is random");
}

void SignalHelper::clearAll()
{
    QMutableHashIterator<QString, SignalData> map_it(_signal_map);
    while (map_it.hasNext())
    {
        map_it.next();

        QString removed_name = map_it.value().name;
        int removed_id       = map_it.value().task_id;
        map_it.remove();
        emit signalRemoved(name2Key(removed_name, removed_id), ALL_VALUES);
    }
    _signal_tree_map.clear();
}

void SignalHelper::removeSignal(const QString& key)
{
    auto it = _signal_map.find(key);
    if (it != _signal_map.end())
    {
        _signal_map.erase(it);
        emit signalRemoved(key, ALL_VALUES);
    }

    QString name;
    int id;
    key2Name(key, name, id);

    _signal_tree_map[id].removeSignal(name.toStdString());

    if (_signal_tree_map[id].isEmpty())
    {
        _signal_tree_map.erase(id);
    }
}

void SignalHelper::removeSignal(const QString& key, int value_idx)
{
    auto it = _signal_map.find(key);
    if (it != _signal_map.end())
    {
        // remove value_idx from the set of used indices
        auto idx_it = it->active_indices.find(value_idx);
        if (idx_it != it->active_indices.end())
        {
            it->active_indices.erase(idx_it);

            // if this was the last value_idx, remove completely
            if (it->active_indices.isEmpty())
            {
                // _signal_map.erase(it);
                removeSignal(key);
                // emit signalRemoved(key, ALL_VALUES);
            }
            else
                emit signalRemoved(key, value_idx);
        }
    }
}

void SignalHelper::printSignals()
{
    PRINT_INFO("---");
    PRINT_INFO("Signals: ");
    for (auto it = _signal_map.begin(); it != _signal_map.end(); ++it)
    {
        PRINT_INFO(it.key().toStdString());
    }
    PRINT_INFO("---");
}

}  // namespace gui
}  // namespace corbo

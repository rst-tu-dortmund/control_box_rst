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
#include <corbo-gui/color_manager.h>
#include <corbo-gui/scope_widget.h>
#include <corbo-gui/utilities.h>
#include <QVBoxLayout>
#include <QVector>

#include <algorithm>
#include <vector>

namespace corbo {
namespace gui {

ScopeWidget::ScopeWidget(SignalHelper::ConstPtr signal_helper, QWidget* parent) : QWidget(parent), _signal_helper(signal_helper)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);

    _layout = new QVBoxLayout(this);
    _layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

    // setup plot widget
    _plot = new QCustomPlot;
    _plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);
    _plot->setMinimumHeight(225);  // minimum height
    // _plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectAxes);
    _plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectAxes | QCP::iMultiSelect);

    // include this section to fully disable antialiasing for higher performance:
    /*
    _plot->setNotAntialiasedElements(QCP::aeAll);
    QFont font;
    font.setStyleStrategy(QFont::NoAntialias);
    _plot->xAxis->setTickLabelFont(font);
    _plot->yAxis->setTickLabelFont(font);
    _plot->legend->setFont(font);
    */

    _layout->addWidget(_plot);

    setupLegend();

    // connect slots that takes care that when an axis is selected, only that direction can be dragged and zoomed:
    connect(_plot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(scopeMousePress()));
    connect(_plot, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(scopeMouseWheel()));
    connect(_plot, SIGNAL(mouseDoubleClick(QMouseEvent*)), this, SLOT(rescaleAxes()));

    // setup policy and connect slot for context menu popup:
    _plot->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(_plot, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(scopeContextMenuRequest(const QPoint&)));

    // enable drag and drop of signals
    setAcceptDrops(true);

    // for resizing the QWidget
    setMouseTracking(true);
}

ScopeWidget::~ScopeWidget() {}

QSize ScopeWidget::sizeHint() const { return QSize(800, 200); }

void ScopeWidget::setupLegend()
{
    if (!_legend_initialized)
    {
        // _plot->setAutoAddPlottableToLegend(false);
        // _plot->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));   // period as decimal separator and comma as thousand
        // separator
        _plot->legend->setVisible(true);
        QFont legend_font = font();   // start out with MainWindow's font..
        legend_font.setPointSize(8);  // and make a bit smaller for legend
        _plot->legend->setFont(legend_font);
        _plot->legend->setBrush(QBrush(QColor(255, 255, 255, 230)));
        _plot->legend->setIconSize(17, 10);
        QPen border_pen = _plot->legend->borderPen();
        border_pen.setColor(QColor("white"));
        _plot->legend->setBorderPen(border_pen);

        // move legend into the axis rect of the main grid layout.
        // create sub layout to generate a small gap between plot cell border and legend border
        QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
        sub_layout->setMargins(QMargins(5, 0, 5, 5));
        sub_layout->addElement(0, 0, _plot->legend);
        sub_layout->addElement(0, 1, new QCPLayoutElement(_plot));  // dummy element in order to force to fillout the complete width
        sub_layout->setColumnStretchFactor(1, 0.0000001);           // set dummy element width to a neglectable width

        _plot->plotLayout()->addElement(1, 0, sub_layout);

        // change fill order to: left -> right
        _plot->legend->setFillOrder(QCPLegend::foColumnsFirst);
        // set legend's row stretch factor very small so it ends up with minimum height:
        _plot->plotLayout()->setRowStretchFactor(1, 0.001);

        _plot->legend->setSelectableParts(QCPLegend::spItems);  // legend box shall not be selectable, only legend items

        _plot->replot();

        _legend_initialized = true;
    }
}

void ScopeWidget::dragEnterEvent(QDragEnterEvent* event)
{
    if (event->mimeData()->hasFormat("text/plain")) event->acceptProposedAction();
}

void ScopeWidget::dropEvent(QDropEvent* event)
{
    QString message = event->mimeData()->text();

    QString key;
    int value_idx = 0;
    fromValueEncodedKey(message, key, value_idx);

    addSignal(key, value_idx);
    event->acceptProposedAction();
}

void ScopeWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton)
    {
        _resize_drag_start_position = event->pos();
        _resize_drag_start_geometry = _plot->geometry();
    }
}

void ScopeWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (!(event->buttons() & Qt::LeftButton))
    {
        // No drag, just change the cursor and return
        if (event->y() >= height() - 5)
        {
            _drag_start_pos = bottom;
            setCursor(Qt::SizeVerCursor);
        }
        else
        {
            // default cursor
            _drag_start_pos = inactive;
            setCursor(Qt::ArrowCursor);
        }
        return;
    }

    // drag to change height

    // consider minimum height of plot field including margin
    if (event->y() <= _plot->minimumHeight() + _layout->margin() + 20) return;

    switch (_drag_start_pos)
    {
        case bottom:
        {
            QPoint window_coords = mapFrom(window(), QPoint(0, window()->height()));
            if (event->y() < window_coords.y() - 3)  // add a small margin which we assume as window bottom
            {
                // setGeometry(geometry().left(), geometry().top(), width(), event->y());
                setFixedHeight(event->y());
            }
            else
            {
                // cursor is at the bottom of the window (parent widget)
                // increase height as long as drag is active
                setFixedHeight(height() + 2);  // otherwise, I was not able to update the parant layout correctly...
            }

            break;
        }
        default:
            break;
    }
}

void ScopeWidget::addSignal(const QString& signal_key, int value_idx)
{
    // new signal
    const SignalHelper::SignalData* signal = _signal_helper->signalData(signal_key);
    if (signal)
    {
        addSignal(signal_key, value_idx, *signal);
    }
    else
        PRINT_WARNING("ScopeWidget::addSignal(): signal " << signal_key.toStdString() << " not found.");
}

void ScopeWidget::addSignal(const QString& signal_key, int value_idx, const SignalHelper::SignalData& signal_data)
{
    QString value_encoded_key = toValueEncodedKey(signal_key, value_idx);
    auto previous_data        = getActiveSignal(value_encoded_key);
    if (previous_data != _active_signals.end())
    {
        // clear potential old graph with same key
        _active_signals.erase(previous_data);
    }

    SignalData& data     = _active_signals[value_encoded_key];
    data.graph_color     = ColorManager::getColor(_plot->plottableCount());
    data.legend_text     = signal_data.name;
    data.zero_order_hold = signal_data.zero_order_hold;
    if (signal_data.dimension > 1)
        data.legend_text += "/" + ((value_idx < signal_data.value_labels.size()) ? signal_data.value_labels[value_idx] : QString::number(value_idx));
    data.value_idx = value_idx;
    data.task_id   = signal_data.task_id;

    // signal type dependent
    data.plottable = nullptr;

    switch (signal_data.signal->getType())
    {
        case SignalType::TimeSeries:
        {
            const TimeSeries* ts = static_cast<const TimeSeriesSignal*>(signal_data.signal.get())->getTimeSeries();
            data.plottable       = addTimeSeriesGraph(*ts, value_idx, data.graph_color, data.legend_text, signal_data.zero_order_hold, true);
            break;
        }
        case SignalType::TimeSeriesSequence:
        {
            data.ts_sequence = std::static_pointer_cast<const TimeSeriesSequenceSignal>(signal_data.signal)->getSequencePtr();
            updateTimeSeriesSequenceGraph(data);  // this method also adds the current TimeSeriesGraph
            break;
        }
        case SignalType::IndexedValuesSet:
        {
            const IndexedValuesSetSignal* set = static_cast<const IndexedValuesSetSignal*>(signal_data.signal.get());
            data.plottable                    = addBoxPlot(*set, data.graph_color, data.legend_text, true);
            break;
        }
        default:
        {
            QMessageBox::warning(this, tr("Cannot plot Signal"), "Signal type not supported yet.");
        }
    }
}

QCPAbstractPlottable* ScopeWidget::addTimeSeriesGraph(const TimeSeries& time_series, int value_idx, const QColor& color, const QString& legend_text,
                                                      bool zero_order_hold, bool replot)
{
    // if (time_series.timeDimension() == 0) return nullptr;

    if (value_idx < 0 || value_idx >= time_series.getValueDimension())
    {
        return nullptr;
    }

    QVector<double> time;
    for (const double& t : time_series.getTime()) time.push_back(t + time_series.getTimeFromStart());

    QVector<double> values;
    TimeSeries::ValuesMatConstMap values_mat = time_series.getValuesMatrixView();
    for (int t = 0; t < values_mat.cols(); ++t)
    {
        values.push_back(values_mat(value_idx, t));
    }

    QCPGraph* graph = _plot->addGraph();

    graph->setPen(QPen(color));
    // graph->setSelectedPen(QPen(QColor("blue"),2));
    graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 3));
    graph->setAntialiasedFill(false);
    graph->setData(time, values, true);
    graph->setName(legend_text);

    if (zero_order_hold) graph->setLineStyle(QCPGraph::lsStepLeft);  // lsStepLeft

    if (time.empty())
        setGraphActive(graph, false);  // not active resp. no values
    else
        setGraphActive(graph, true);

    graph->rescaleAxes(_plot->graphCount() > 1);
    // if (replot) _plot->replot(QCustomPlot::rpQueuedReplot);
    if (replot) _plot->replot();
    return graph;
}

void ScopeWidget::updateTimeSeriesGraph(SignalData& data, double t, double value, bool enlarge_axis, bool replot)
{
    QCPGraph* graph = dynamic_cast<QCPGraph*>(data.plottable);
    if (!graph)
    {
        PRINT_ERROR_NAMED("cannot update plottable since is not of type OCPGraph");
        return;
    }

    graph->addData(t, value);
    if (enlarge_axis) graph->rescaleAxes(true);
    // TODO(roesmann) do not replot for every update -> create a timer in a separate thread to trigger replot with a slower rate!!
    // if (replot) _plot->replot(QCustomPlot::rpQueuedReplot);
    if (replot) _plot->replot();
}

void ScopeWidget::updateTimeSeriesGraph(SignalData& data, const TimeSeries& time_series, bool enlarge_axis, bool replace_data, bool replot)
{
    if (!data.plottable) return;

    QCPGraph* graph = dynamic_cast<QCPGraph*>(data.plottable);
    if (!graph)
    {
        PRINT_ERROR_NAMED("cannot update plottable since is not of type OCPGraph");
        return;
    }

    if (data.value_idx < 0 || data.value_idx >= time_series.getValueDimension()) return;

    QVector<double> time;
    for (const double& t : time_series.getTime()) time.push_back(t + time_series.getTimeFromStart());

    QVector<double> values;
    TimeSeries::ValuesMatConstMap values_mat = time_series.getValuesMatrixView();
    for (int t = 0; t < values_mat.cols(); ++t)
    {
        values.push_back(values_mat(data.value_idx, t));
    }

    // data.graph->addData(t, value);
    if (replace_data)
    {
        graph->setData(time, values);
    }
    else
    {
        graph->addData(time, values);
    }

    if (graph->data()->isEmpty())
        setGraphActive(graph, false);  // not active resp. no values
    else
        setGraphActive(graph, true);

    if (enlarge_axis) graph->rescaleAxes(true);
    // if (replot) _plot->replot(QCustomPlot::rpQueuedReplot);
    if (replot) _plot->replot();
}

void ScopeWidget::addMeasurement(const QString& key, Measurement::ConstPtr measurement, const SignalHelper::SignalData& data)
{
    // get all graphs that correspond to that signal
    QVector<QHash<QString, ScopeWidget::SignalData>::iterator> active_signals = findActiveSignal(key, SearchType::Exact);
    for (int i = 0; i < active_signals.size(); ++i)
    {
        // first check if the current signal is already initialized or just reserved
        if (active_signals[i].value().plottable == nullptr)
        {
            // just reserved: we have to initialize the new signal
            int value_idx;
            QString singal_key;
            if (!fromValueEncodedKey(active_signals[i].key(), singal_key, value_idx)) continue;
            addSignal(singal_key, value_idx, data);
            return;
        }

        QCPGraph* graph = dynamic_cast<QCPGraph*>(active_signals[i].value().plottable);
        if (!graph) continue;

        if (active_signals[i].value().value_idx < measurement->getValues().size())
        {
            updateTimeSeriesGraph(active_signals[i].value(), measurement->getTime(), measurement->getValues()[active_signals[i].value().value_idx],
                                  false, false);

            if (!isGraphActive(graph)) setGraphActive(graph, true);
        }
        else
        {
            PRINT_ERROR(
                "Received values for plotting a time-series to scope, but the dimension is smaller than the values "
                "dimension.");
        }
    }
    // double dur = (measurement->header.time - _last_signal_header_time).toSec();
    // if (dur > 0.05)
    //{
    //    _last_signal_header_time = measurement->header.time;
    _plot->rescaleAxes(true);
    _plot->replot(QCustomPlot::rpQueuedReplot);
    //}
}

void ScopeWidget::setPreviewTime(double preview_time)
{
    _current_preview_time = preview_time;

    for (auto it = _active_signals.begin(); it != _active_signals.end(); ++it)
    {
        if (it->ts_sequence) updateTimeSeriesSequenceGraph(*it, false);
    }
    _plot->replot();
}

void ScopeWidget::updateTimeSeriesSequenceGraph(SignalData& data, bool replot)
{
    if (!data.ts_sequence) return;

    // we assume that data has already been filled appropriately (everything except graph must be set, e.g. color, value_idx, legend_text,
    // ...).

    TimeSeries::ConstPtr current_ts;

    // find TimeSeries to be plotted (w.r.t. _current_preview_time)
    // we assume that the TimeSeries sequence is already ordered according to its time_from_start members.
    auto found_it = std::find_if(data.ts_sequence->getSequence().rbegin(), data.ts_sequence->getSequence().rend(),
                                 [this](const TimeSeries::Ptr& ts) { return ts->getTimeFromStart() <= _current_preview_time; });  // reverse iterators
    if (found_it != data.ts_sequence->getSequence().rend())
    {
        // found
        current_ts = *found_it;
    }
    else
    {
        // plot first available time series:
        if (data.ts_sequence->getSequence().empty())
        {
            PRINT_ERROR("ScopeWidget::updateTimeSeriesSequenceGraph(): time series sequence is empty.");
            return;
        }
        current_ts = data.ts_sequence->getSequence().front();
    }

    // check if we have plotted a TimeSeries of this sequence before
    if (data.plottable)
    {
        // just update the values in order to keep color, labels and legend key.
        updateTimeSeriesGraph(data, *current_ts, true, true, replot);
    }
    else
    {
        // add a complete new TimeSeriesPlot
        data.plottable = addTimeSeriesGraph(*current_ts, data.value_idx, data.graph_color, data.legend_text, data.zero_order_hold, replot);
    }
}

QCPAbstractPlottable* ScopeWidget::addBoxPlot(const IndexedValuesSetSignal& indexed_values_set, const QColor& color, const QString& legend_text,
                                              bool replot)
{
    if (indexed_values_set.isEmpty()) return nullptr;

    QCPStatisticalBox* box = new QCPStatisticalBox(_plot->xAxis, _plot->yAxis);

    for (const auto& item : indexed_values_set.getData())
    {
        // QVector<double> data = QVector<double>::fromStdVector(item.second);
        std::vector<double> data = item.second;

        // get statistic properties of the data
        double minimum = *std::min_element(data.begin(), data.end());
        double maximum = *std::max_element(data.begin(), data.end());

        int num_quartile1 = data.size() / 4;
        int num_quartile2 = data.size() / 2;
        int num_quartile3 = num_quartile1 + num_quartile2;

        std::nth_element(data.begin(), data.begin() + num_quartile1, data.end());
        std::nth_element(data.begin() + num_quartile1 + 1, data.begin() + num_quartile2, data.end());
        std::nth_element(data.begin() + num_quartile2 + 1, data.begin() + num_quartile3, data.end());

        double lower_quartile = data[num_quartile1];
        double median         = data[num_quartile2];
        double upper_quartile = data[num_quartile3];

        box->addData(item.first, minimum, lower_quartile, median, upper_quartile, maximum);  // no outliers yet
    }

    box->setName(legend_text);
    box->setPen(QPen(color));

    // QBrush box_brush(QColor(60, 60, 255, 100));
    // box_brush.setStyle(Qt::Dense6Pattern); // make it look oldschool
    // box->setBrush(box_brush);

    _plot->rescaleAxes();

    if (replot) _plot->replot();

    return box;
}

void ScopeWidget::removeSignal(const QString& key, int value_idx)
{
    int removed = 0;
    if (value_idx != SignalHelper::ALL_VALUES)
    {
        QString value_encoded_key = toValueEncodedKey(key, value_idx);
        removed                   = _active_signals.remove(value_encoded_key);
    }
    else
    {
        auto iterators = findActiveSignal(key, SearchType::Exact);
        for (auto& it : iterators) _active_signals.erase(it);
        removed = (int)iterators.size();
    }

    if (removed > 0)
    {
        _plot->rescaleAxes(true);
        _plot->replot();
    }
}

void ScopeWidget::initializeTask(int task_id, bool inherit_signals)
{
    if (inherit_signals)
    {
        QString signal_key;
        for (auto it = _active_signals.begin(); it != _active_signals.end(); ++it)
        {
            int value_idx;
            QString signal_key;
            if (!fromValueEncodedKey(it.key(), signal_key, value_idx)) continue;
            QString name;
            int prev_task_id;
            if (!SignalHelper::key2Name(signal_key, name, prev_task_id)) continue;

            if (prev_task_id == task_id - 1)  // only inherit if also obtained in the last task run
            {
                // new task found (this might only happen if search_type == SearchType::IgnorePrefix)
                QString new_key = SignalHelper::name2Key(name, task_id);
                // reserve field which should be initialized in one of the add-signals methods
                SignalData& data = _active_signals[toValueEncodedKey(new_key, value_idx)];
                return;
            }
        }
    }
}

void ScopeWidget::rescaleAxes()
{
    _plot->rescaleAxes(false);
    // _plot->replot(QCustomPlot::rpQueuedReplot);
    _plot->replot();
}

bool ScopeWidget::hasSignal(const QString& key, int value_idx) const
{
    auto it = _active_signals.find(toValueEncodedKey(key, value_idx));
    if (it != _active_signals.end()) return true;

    return false;
}

QHash<QString, ScopeWidget::SignalData>::iterator ScopeWidget::getActiveSignal(const QString& value_encoded_key)
{
    auto it = _active_signals.find(value_encoded_key);
    if (it != _active_signals.end()) return it;

    return _active_signals.end();
}

QHash<QString, ScopeWidget::SignalData>::iterator ScopeWidget::getActiveSignal(const QString& key, int value_idx)
{
    return getActiveSignal(toValueEncodedKey(key, value_idx));
}

QVector<QHash<QString, ScopeWidget::SignalData>::iterator> ScopeWidget::findActiveSignal(const QString& key, SearchType search_type)
{
    QVector<QHash<QString, ScopeWidget::SignalData>::iterator> active_list;

    QString signal_key;
    int value_idx = 0;
    for (auto it = _active_signals.begin(); it != _active_signals.end(); ++it)
    {
        if (fromValueEncodedKey(it.key(), signal_key, value_idx))
        {
            bool found = false;
            switch (search_type)
            {
                case SearchType::Exact:
                {
                    found = (signal_key.compare(key) == 0);
                    break;
                }
                case SearchType::IgnorePrefix:
                {
                    QStringList token_source = key.split(util::SIGNAL_NAMESPACE_PREFIX_DELIMITER);
                    QStringList token_active = signal_key.split(util::SIGNAL_NAMESPACE_PREFIX_DELIMITER);
                    found                    = (token_active.back().compare(token_source.back()) == 0);
                    break;
                }
                default:
                {
                    PRINT_WARNING("ScopeWidget::findActiveSignal(): selected search type not implemented.");
                }
            }
            if (found) active_list.push_back(it);
        }
    }
    return active_list;
}

QString ScopeWidget::toValueEncodedKey(const QString& key, int value_idx) const
{
    return key + util::SIGNAL_NAMESPACE_SUFFIX_DELIMITER + QString::number(value_idx);
}

bool ScopeWidget::fromValueEncodedKey(const QString& value_encoded_key, QString& key, int& value_idx)
{
    QStringList tokens = value_encoded_key.split(util::SIGNAL_NAMESPACE_SUFFIX_DELIMITER);
    PRINT_WARNING_COND(tokens.size() > 2, "QString::fromValueEnvodedKey: more than 2 tokens found...");

    bool ret_val = true;

    key = tokens.front();
    if (tokens.size() > 1)
    {
        bool ok;
        value_idx = tokens.back().toInt(&ok);
        if (!ok)
        {
            PRINT_ERROR("QString::fromValueEnvodedKey: cannot detect value index in " << value_encoded_key.toStdString());
            value_idx = 0;
            ret_val   = false;
        }
    }
    else
    {
        value_idx = 0;
        ret_val   = false;
        PRINT_WARNING("QString::fromValueEnvodedKey: no value_idx found in " << value_encoded_key.toStdString() << ". Setting value_idx to 0");
    }
    return ret_val;
}

void ScopeWidget::scopeMousePress()
{
    // if an axis is selected, only allow the direction of that axis to be dragged
    // if no axis is selectedor or if both axes are selected (multi-select), both directions may be dragged

    if (_plot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) && !_plot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        _plot->axisRect()->setRangeDrag(_plot->xAxis->orientation());
    else if (!_plot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) && _plot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        _plot->axisRect()->setRangeDrag(_plot->yAxis->orientation());
    else
        _plot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);
}

void ScopeWidget::scopeMouseWheel()
{
    // if an axis is selected, only allow the direction of that axis to be zoomed
    // if no axis is selected or if both axes are selected (multi-select), both directions may be zoomed

    if (_plot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) && !_plot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        _plot->axisRect()->setRangeZoom(_plot->xAxis->orientation());
    else if (!_plot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) && _plot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        _plot->axisRect()->setRangeZoom(_plot->yAxis->orientation());
    else
        _plot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
}

void ScopeWidget::scopeContextMenuRequest(const QPoint& point)
{
    // QMenu* menu = new QMenu(this);
    // menu->setAttribute(Qt::WA_DeleteOnClose);
    QMenu menu;

    // context menu on legend requested
    if (_plot->legend->selectTest(point, false) >= 0)
    {
        // find legend for active graphs
        auto it = _active_signals.begin();
        while (it != _active_signals.end())
        {
            QCPPlottableLegendItem* graph_item = _plot->legend->itemWithPlottable(it.value().plottable);
            if (graph_item)
            {
                if (graph_item->selectTest(point, false) >= 0)
                {
                    QAction* remove_action = new QAction(tr("Remove signal"), this);
                    connect(remove_action, &QAction::triggered, [this, it]() {
                        _active_signals.erase(it);
                        _plot->replot();
                    });
                    menu.addAction(remove_action);
                    break;
                }
            }
            ++it;
        }
    }
    else  // general context menu
    {
        QAction* rescale_axes = new QAction(tr("Rescale axes"), this);
        connect(rescale_axes, &QAction::triggered, [this]() { rescaleAxes(); });
        menu.addAction(rescale_axes);

        menu.addSeparator();

        bool signal_options = false;
        if (_plot->selectedGraphs().size() > 0)
        {
            menu.addAction("Remove selected signals", this, SLOT(removeSelectedSignals()));
            signal_options = true;
        }
        if (_plot->graphCount() > 0)
        {
            menu.addAction("Remove all signals", this, SLOT(removeAllSignals()));
            signal_options = true;
        }

        if (signal_options) menu.addSeparator();

        QAction* close_scope = new QAction(tr("Close scope"), this);
        connect(close_scope, &QAction::triggered, [this]() { deleteLater(); });
        menu.addAction(close_scope);
    }

    // menu->popup(_plot->mapToGlobal(point));
    menu.exec(_plot->mapToGlobal(point));
}

void ScopeWidget::removeSelectedSignals()
{
    QMutableHashIterator<QString, SignalData> it(_active_signals);
    bool update = false;
    while (it.hasNext())
    {
        it.next();
        if (it.value().plottable && it.value().plottable->selected())
        {
            it.remove();
            update = true;
        }
    }
    if (update) _plot->replot();
}

void ScopeWidget::removeAllSignals()
{
    _active_signals.clear();
    _plot->replot();
}

bool ScopeWidget::isGraphActive(const QCPGraph* graph) const
{
    if (!graph) return false;

    QVariant active = graph->property("active");

    if (active.isValid() && active.toBool()) return true;
    return false;
}

void ScopeWidget::setGraphActive(QCPGraph* graph, bool active)
{
    // set to active (adjust legend font)
    QCPPlottableLegendItem* item = _plot->legend->itemWithPlottable(graph);
    if (item)
    {
        if (active)
            item->setTextColor(QColor(0, 0, 0));
        else
            item->setTextColor(QColor(0, 0, 0, 50));
        graph->setProperty("active", active);
    }
}

}  // namespace gui
}  // namespace corbo

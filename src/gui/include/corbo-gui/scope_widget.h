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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_WIDGET_H_

#include <corbo-communication/messages/core/signals.pb.h>
#include <corbo-gui/signal_collection_widget.h>
#include <corbo-gui/utilities.h>

#include <qcustomplot.h>

#include <QColor>
#include <QHash>
#include <QString>
#include <QVector>
#include <QWidget>

namespace corbo {
namespace gui {

class ScopeWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit ScopeWidget(SignalHelper::ConstPtr signal_helper, QWidget* parent = 0);
    virtual ~ScopeWidget();

    enum ResizeDragStratPositions { bottom, top, inactive };

    QSize sizeHint() const override;

    bool hasSignal(const QString& key, int value_idx) const;

 public slots:
    void addSignal(const QString& signal_key, int value_idx);
    void addSignal(const QString& signal_key, int value_idx, const SignalHelper::SignalData& signal_data);

    void addMeasurement(const QString& key, Measurement::ConstPtr measurement, const SignalHelper::SignalData& data);
    void setPreviewTime(double preview_time);
    void removeSignal(const QString& key, int value_idx);
    void initializeTask(int task_id, bool inherit_signals);
    void rescaleAxes();
    void replot() { _plot->replot(); }

    void removeSelectedSignals();
    void removeAllSignals();

 protected:
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;

    struct SignalData
    {
        ~SignalData()
        {
            // erase graph from plot
            if (plottable && plottable->parentPlot()) plottable->parentPlot()->removePlottable(plottable);
        }
        int value_idx        = 0;
        int task_id          = -1;
        bool zero_order_hold = false;
        QColor graph_color   = QColor("blue");
        QString legend_text;
        QCPAbstractPlottable* plottable = nullptr;
        TimeSeriesSequence::ConstPtr ts_sequence;  // only set if signal contains a TimeSeriesSequence
    };
    QHash<QString, SignalData>::iterator getActiveSignal(const QString& value_encoded_key);
    QHash<QString, SignalData>::iterator getActiveSignal(const QString& key, int value_idx);
    QVector<QHash<QString, SignalData>::iterator> findActiveSignal(const QString& key, SearchType search_type);
    QString toValueEncodedKey(const QString& key, int value_idx) const;
    bool fromValueEncodedKey(const QString& value_encoded_key, QString& key, int& value_idx);

    bool isGraphActive(const QCPGraph* graph) const;
    void setGraphActive(QCPGraph* graph, bool active);

    void setupLegend();

 protected slots:
    QCPAbstractPlottable* addTimeSeriesGraph(const TimeSeries& time_series, int value_idx, const QColor& color, const QString& legend_text,
                                             bool zero_order_hold, bool replot);
    void updateTimeSeriesGraph(SignalData& data, double t, double value, bool enlarge_axis, bool replot = true);
    void updateTimeSeriesGraph(SignalData& data, const TimeSeries& time_series, bool enlarge_axis, bool replace_data, bool replot = true);
    void updateTimeSeriesSequenceGraph(SignalData& data, bool replot = true);

    QCPAbstractPlottable* addBoxPlot(const IndexedValuesSetSignal& indexed_values_set, const QColor& color, const QString& legend_text,
                                     bool replot = true);

    void scopeMousePress();
    void scopeMouseWheel();
    void scopeContextMenuRequest(const QPoint& point);

 private:
    QCustomPlot* _plot;

    QHash<QString, SignalData> _active_signals;

    SignalHelper::ConstPtr _signal_helper;

    bool _legend_initialized = false;

    Time _last_signal_header_time = Time(0);

    ResizeDragStratPositions _drag_start_pos;
    QPoint _resize_drag_start_position;
    QRect _resize_drag_start_geometry;

    double _current_preview_time = 0.0;  //!< determine which TimeSeries of a TimeSeriesSequence should be plotted

    QVBoxLayout* _layout;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_WIDGET_H_

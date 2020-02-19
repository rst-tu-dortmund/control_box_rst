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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_LABEL_SLIDER_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_LABEL_SLIDER_WIDGET_H_

#include <QLabel>
#include <QSlider>

#include <limits>

namespace corbo {
namespace gui {

class LabelSliderWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit LabelSliderWidget(QWidget* parent = 0) : LabelSliderWidget("", 0.0, 1, 2, parent) {}

    LabelSliderWidget(const QString& label, double default_value, int resolution = 1, int display_decimal_places = 2, QWidget* parent = 0);
    ~LabelSliderWidget() {}

    QSize sizeHint() const override;

    void setLabel(const QString& label) { _label->setText(label); }
    void setValue(double value);

    void setResolution(int resolution) { _resolution = resolution; }

    void setMinMax(double min, double max);

    void setDisplayDecimalPlaces(int num_places) { _decimal_places = num_places; }

    QSlider* widgetSlider() { return _slider; }
    const QSlider* widgetSlider() const { return _slider; }

    QLabel* widgetLabel() { return _label; }
    const QLabel* widgetLabel() const { return _label; }

 signals:
    void valueChanged(double value);

 private slots:
    void updateDisplay(int value);

 private:
    QSlider* _slider;
    QLabel* _label;
    QLabel* _display;
    int _resolution       = 1;
    int _decimal_places   = 2;
    double _min_dbl_value = std::numeric_limits<double>::min();
    double _max_dbl_value = std::numeric_limits<double>::max();
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_LABEL_SLIDER_WIDGET_H_

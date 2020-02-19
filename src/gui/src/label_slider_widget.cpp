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

#include <corbo-gui/label_slider_widget.h>
#include <QHBoxLayout>

namespace corbo {
namespace gui {

LabelSliderWidget::LabelSliderWidget(const QString& label, double default_value, int resolution, int display_decimal_places, QWidget* parent)
    : QWidget(parent)
{
    setSizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Preferred);

    QHBoxLayout* box = new QHBoxLayout(this);
    box->setContentsMargins(5, 0, 0, 0);
    box->setAlignment(Qt::AlignLeft);
    _label = new QLabel(label);
    box->addWidget(_label);
    _display = new QLabel;
    _display->setMinimumWidth(15);
    _slider = new QSlider;
    _slider->setOrientation(Qt::Horizontal);
    _resolution     = resolution;
    _decimal_places = display_decimal_places;
    setValue(default_value);  // display is updated as well, must be initialized first
    box->addWidget(_slider);
    connect(_slider, SIGNAL(valueChanged(int)), this, SLOT(updateDisplay(int)));

    _label->setBuddy(_slider);
    box->addWidget(_display, 0, Qt::AlignRight);
}

QSize LabelSliderWidget::sizeHint() const { return QSize(250, 20); }

void LabelSliderWidget::setValue(double value)
{
    int act_value = int(value * _resolution);
    _slider->setValue(act_value);
    updateDisplay(act_value);
}

void LabelSliderWidget::setMinMax(double min, double max)
{
    _min_dbl_value = min;
    _max_dbl_value = max;
    _slider->setMinimum(int(min * _resolution));
    _slider->setMaximum(int(max * _resolution));
}

void LabelSliderWidget::updateDisplay(int value)
{
    double double_val = (double)value / (double)_resolution;
    // bound value, since there could be some rounding errors due to the resolution
    if (double_val < _min_dbl_value)
        double_val = _min_dbl_value;
    else if (double_val > _max_dbl_value)
        double_val = _max_dbl_value;

    _display->setText(QString::number(double_val, 'f', _decimal_places));
    emit valueChanged(double_val);
}

}  // namespace gui
}  // namespace corbo

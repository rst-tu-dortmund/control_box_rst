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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SLIDER_CENTER_STUCK_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SLIDER_CENTER_STUCK_WIDGET_H_

#include <QSlider>

#include <cmath>

namespace corbo {
namespace gui {

class SliderCenterStuck : public QSlider
{
    Q_OBJECT
 public:
    explicit SliderCenterStuck(QWidget* parent = nullptr) : QSlider(parent) { connect(this, SIGNAL(sliderReleased()), this, SLOT(onRelease())); }

    explicit SliderCenterStuck(Qt::Orientation orientation, QWidget* parent = nullptr) : QSlider(orientation, parent)
    {
        connect(this, SIGNAL(sliderReleased()), this, SLOT(onRelease()));
    }

    int getCenter() { return std::round(0.5 * (double)(maximum() + minimum())); }
    void setToCenter() { this->setSliderPosition(getCenter()); }

    void setStuckToCenter(int abs_distance) { _stuck_center_dist = abs_distance; }

 protected slots:
    void onRelease()
    {
        if (abs(value() - getCenter()) < _stuck_center_dist)
        {
            setToCenter();
        }
    }

 private:
    int _stuck_center_dist = 0;
};

/*
class NonWrappingDial : public QDial
{
    Q_OBJECT
 public:
    explicit NonWrappingDial(QWidget* parent = NULL) : QDial(parent)
    {
        connect(this, SIGNAL(actionTriggered(int)), this, SLOT(onAction(int)));
        connect(this, SIGNAL(sliderReleased()), this, SLOT(onRelease()));
    }

    void setStuckToCenter(int abs_distance) { _stuck_center_dist = abs_distance; }

    void setToCenter() { this->setSliderPosition(getCenter()); }

    int getCenter() { return std::round(0.5 * (double)(maximum() + minimum())); }

 protected slots:
    void onAction(int val)
    {
        static const int minDistance = 1;
        if (val == QAbstractSlider::SliderMove)
        {
            if (value() == maximum() && sliderPosition() < maximum() - minDistance)
            {
                this->setSliderPosition(maximum());
            }
            else if (value() == minimum() && sliderPosition() > minimum() + minDistance)
            {
                this->setSliderPosition(minimum());
            }
        }
    }
    void onRelease()
    {
        if (abs(value() - getCenter()) < _stuck_center_dist)
        {
            setToCenter();
        }
    }

 private:
    int _stuck_center_dist = 0;
};
*/

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SLIDER_CENTER_STUCK_WIDGET_H_

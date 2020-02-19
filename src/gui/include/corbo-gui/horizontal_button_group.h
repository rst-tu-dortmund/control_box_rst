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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_HORIZONTAL_BUTTON_GROUP_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_HORIZONTAL_BUTTON_GROUP_H_

#include <QAbstractButton>
#include <QHBoxLayout>
#include <QLabel>

namespace corbo {
namespace gui {

class HoriontalButtonGroup : public QWidget
{
    Q_OBJECT
 public:
    explicit HoriontalButtonGroup(bool exclusive, QWidget* parent = 0) : HoriontalButtonGroup("", exclusive, parent) {}
    explicit HoriontalButtonGroup(const QString& label, bool exclusive = false, QWidget* parent = 0);
    ~HoriontalButtonGroup(){};

    QAbstractButton* addButton(bool checked, const QString& cb_label = "");

    void setLabel(const QString& label) { _label->setText(label); }
    virtual QSize sizeHint() const;
    QLabel* widgetLabel() { return _label; }
    const QLabel* widgetLabel() const { return _label; }
    QAbstractButton* widgetButton(int idx);
    QVector<bool> buttonStates() const;
    void setButtonStates(const QVector<bool>& states);
    int noButtons() const { return _boxes.size(); }
    bool isExclusive() const { return _exclusive; }

 private:
    QList<QAbstractButton*> _boxes;
    QLabel* _label;
    QHBoxLayout* _layout;
    bool _exclusive = false;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_HORIZONTAL_BUTTON_GROUP_H_

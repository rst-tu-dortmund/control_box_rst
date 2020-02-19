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

#include <corbo-gui/horizontal_button_group.h>
#include <QCheckBox>
#include <QRadioButton>

namespace corbo {
namespace gui {

HoriontalButtonGroup::HoriontalButtonGroup(const QString& label, bool exclusive, QWidget* parent) : QWidget(parent)
{
    setSizePolicy(QSizePolicy::Policy::Maximum, QSizePolicy::Policy::Maximum);

    _exclusive = exclusive;
    _layout    = new QHBoxLayout(this);
    _layout->setContentsMargins(0, 0, 0, 0);
    _layout->setAlignment(Qt::AlignLeft);
    _label = new QLabel(label);
    _layout->addWidget(_label);
    _layout->addSpacing(15);
}

QAbstractButton* HoriontalButtonGroup::addButton(bool checked, const QString& cb_label)
{
    if (!_exclusive)
    {
        QCheckBox* cb = new QCheckBox(cb_label);
        cb->setChecked(checked);
        _layout->addWidget(cb);
        _boxes.push_back(cb);
        return cb;
    }
    QRadioButton* cb = new QRadioButton(cb_label);
    cb->setChecked(checked);
    _layout->addWidget(cb);
    _boxes.push_back(cb);
    return cb;
}

QAbstractButton* HoriontalButtonGroup::widgetButton(int idx)
{
    if (idx >= _boxes.size()) return nullptr;
    return _boxes[idx];
}
QVector<bool> HoriontalButtonGroup::buttonStates() const
{
    QVector<bool> states;
    for (QAbstractButton* btn : _boxes)
    {
        states.push_back(btn->isChecked());
    }
    return states;
}
void HoriontalButtonGroup::setButtonStates(const QVector<bool>& states)
{
    if (states.size() != noButtons())
    {
        // PRINT_INFO("Cannot set button states: number of buttons does not correspond to number of input states");
        return;
    }
    int idx = 0;
    for (bool state : states)
    {
        _boxes[idx]->setChecked(state);
        ++idx;
    }
}

QSize HoriontalButtonGroup::sizeHint() const { return QSize(50, 5); }

}  // namespace gui
}  // namespace corbo

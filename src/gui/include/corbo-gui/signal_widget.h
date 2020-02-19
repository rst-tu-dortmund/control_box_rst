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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_WIDGET_H_

#include <QHBoxLayout>
#include <QMenu>
#include <QString>
#include <QWidget>

namespace corbo {
namespace gui {

class SignalWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit SignalWidget(const QString& title, const QString& key, int value_idx, QWidget* parent = 0);
    virtual ~SignalWidget();

    QSize sizeHint() const override;

    QString key() const { return _key; }
    int valueIdx() const { return _value_idx; }

 protected:
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;

 private:
    QHBoxLayout* _layout;

    QString _key;
    int _value_idx = 0;

    QMenu* _context_menu;

    QPoint _drag_start_position;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_WIDGET_H_

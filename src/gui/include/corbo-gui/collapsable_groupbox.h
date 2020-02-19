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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_COLLAPSABLE_GROUPBOX_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_COLLAPSABLE_GROUPBOX_H_

#include <QGroupBox>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidget>

namespace corbo {
namespace gui {

// based on http://stackoverflow.com/questions/32476006/how-to-make-an-expandable-collapsable-section-widget-in-qt
// but we are using a groupbox instead of a scrollwidget and fixed sizes (otherwise we get problems with nested groups)
class CollapsableGroupBox : public QWidget
{
    Q_OBJECT

 public:
    explicit CollapsableGroupBox(const QString& title, QWidget* parent = 0);

    bool isCollapsed() { return _collapsed; }

    QGroupBox* groupBox() { return _content; }

 public slots:
    void setCollapsed(bool collapsed);

 protected:
    void createTitle(const QString& title);
    void createContentArea();

 private:
    QVBoxLayout* _layout;
    QToolButton* _button;
    QGroupBox* _content;

    bool _collapsed;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_COLLAPSABLE_GROUPBOX_H_

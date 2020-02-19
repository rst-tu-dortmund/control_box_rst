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

#include <corbo-gui/collapsable_groupbox.h>
#include <QFrame>
#include <QHBoxLayout>
#include <QPropertyAnimation>
#include <QToolButton>
#include <QVBoxLayout>

#include <corbo-core/console.h>

namespace corbo {
namespace gui {

CollapsableGroupBox::CollapsableGroupBox(const QString& title, QWidget* parent) : QWidget(parent)
{
    _layout = new QVBoxLayout(this);
    _layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    _layout->setContentsMargins(0, 0, 0, 0);
    _layout->setSpacing(0);

    createTitle(title);
    createContentArea();

    setCollapsed(false);

    connect(_button, &QToolButton::clicked, [this](bool checked) { setCollapsed(!checked); });
}

void CollapsableGroupBox::createTitle(const QString& title)
{
    QHBoxLayout* title_layout = new QHBoxLayout;
    title_layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    title_layout->setContentsMargins(0, 0, 0, 0);

    // create button plus label
    _button = new QToolButton;
    _button->setStyleSheet("QToolButton { border: none; }");
    _button->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    _button->setArrowType(Qt::ArrowType::RightArrow);
    _button->setText(title);
    _button->setCheckable(true);
    _button->setChecked(false);
    title_layout->addWidget(_button);

    // create horizontal line
    QFrame* hline = new QFrame;
    hline->setFrameShape(QFrame::HLine);
    hline->setFrameShadow(QFrame::Sunken);
    hline->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);
    title_layout->addWidget(hline);

    _layout->addLayout(title_layout);
}

void CollapsableGroupBox::createContentArea()
{
    _content = new QGroupBox;
    _content->setFlat(true);
    _content->setStyleSheet("QGroupBox { border: none; }");

    // set collapsed
    _content->hide();

    _layout->addWidget(_content);
}

void CollapsableGroupBox::setCollapsed(bool collapsed)
{
    if (_button->isChecked() == collapsed)
    {
        _button->setChecked(!collapsed);
        return;  // retrigger signal to synchronize button status if this method has been called programmatically
    }

    _collapsed = collapsed;

    _button->setArrowType(_collapsed ? Qt::ArrowType::RightArrow : Qt::ArrowType::DownArrow);

    // int collapsed_height = sizeHint().height() - _content->maximumHeight();
    if (_collapsed)
    {
        _content->hide();
    }
    else
    {
        _content->show();
    }
}

}  // namespace gui
}  // namespace corbo

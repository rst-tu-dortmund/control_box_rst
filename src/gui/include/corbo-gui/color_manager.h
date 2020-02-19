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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_COLOR_MANAGER_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_COLOR_MANAGER_H_

#include <QColor>
#include <vector>

namespace corbo {
namespace gui {

class ColorManager
{
 public:
    static void addColor(const QColor& color) { _colors.push_back(color); }

    static const QColor& getColor(int index)
    {
        index = index % (int)_colors.size();
        return _colors.at(index);
    }

    static QColor getColor(int index, int alpha)
    {
        index        = index % (int)_colors.size();
        QColor color = _colors.at(index);
        color.setAlpha(alpha);
        return color;
    }

 private:
    static std::vector<QColor> _colors;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_COLOR_MANAGER_H_

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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_TEXT_STYLE_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_TEXT_STYLE_H_

#include <iostream>
#include <sstream>
#include <string>

namespace corbo {

/**
 * @brief Color codes for colored text in consoles/terminals (if supported)
 *
 * see http://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
 * and http://misc.flogisoft.com/bash/tip_colors_and_formatting
 */
enum class TextColorCode {
    FG_DEFAULT       = 39,
    FG_RED           = 31,
    FG_GREEN         = 32,
    FG_BLUE          = 34,
    FG_BLACK         = 30,
    FG_YELLOW        = 33,
    FG_MAGENTA       = 35,
    FG_CYAN          = 36,
    FG_LIGHT_GRAY    = 37,
    FG_DARK_GRAY     = 90,
    FG_LIGHT_RED     = 91,
    FG_LIGHT_GREEN   = 92,
    FG_LIGHT_YELLOW  = 93,
    FG_LIGHT_BLUE    = 94,
    FG_LIGHT_MAGENTA = 95,
    FG_LIGHT_CYAN    = 96,
    FG_WHITE         = 97,
    BG_RED           = 41,
    BG_GREEN         = 42,
    BG_BLUE          = 44,
    BG_DEFAULT       = 49
};

//! Output stream operator to change the color of subsequent stream output
inline std::ostream& operator<<(std::ostream& os, TextColorCode code)
{
#ifdef __linux__
    return os << "\033[" << static_cast<int>(code) << "m";
#else  // not supported
    return os;
#endif
}

//! Return input text as colored text for terminal output (note: the terminal must support colored text)
inline std::string coloredText(const std::string& text, TextColorCode color)
{
#ifdef __linux__
    std::stringstream ss;
    ss << color << text << TextColorCode::FG_DEFAULT;
    return ss.str();
#else  // not supported
    return text;
#endif
}
}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_TEXT_STYLE_H_

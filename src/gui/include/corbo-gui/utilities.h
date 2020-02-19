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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_UTILITIES_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_UTILITIES_H_

#include <corbo-core/types.h>
#include <corbo-core/utilities.h>

#include <QFrame>
#include <QRegExp>
#include <QString>

#include <corbo-core/console.h>

#include <cmath>
#include <string>
#include <type_traits>
#include <vector>

namespace corbo {
namespace gui {

enum class SearchType { Exact, IgnorePrefix };

namespace util {

constexpr const char SIGNAL_NAMESPACE_PREFIX_DELIMITER[] = "__";
constexpr const char SIGNAL_NAMESPACE_SUFFIX_DELIMITER[] = "::";

inline double qstring_to_double(const QString& string, bool* ok = nullptr)
{
    if (ok) *ok = true;
    // TODO(roesmann): create dictionary, e.g. with QHash
    if (string.compare("inf", Qt::CaseInsensitive) == 0) return CORBO_INF_DBL;
    if (string.compare("-inf", Qt::CaseInsensitive) == 0) return -CORBO_INF_DBL;
    if (string.compare("pi", Qt::CaseInsensitive) == 0) return M_PI;
    if (string.compare("-pi", Qt::CaseInsensitive) == 0) return -M_PI;
    if (string.compare("2pi", Qt::CaseInsensitive) == 0) return M_2_PI;
    if (string.compare("-2pi", Qt::CaseInsensitive) == 0) return -M_2_PI;
    if (string.compare("pi/2", Qt::CaseInsensitive) == 0) return M_PI_2;
    if (string.compare("-pi/2", Qt::CaseInsensitive) == 0) return -M_PI_2;
    if (string.compare("pi/4", Qt::CaseInsensitive) == 0) return M_PI_4;
    if (string.compare("-pi/4", Qt::CaseInsensitive) == 0) return -M_PI_4;
    return corbo::util::clamp(string.toDouble(ok), -CORBO_INF_DBL, CORBO_INF_DBL);
}

inline int qstring_to_int(const QString& string, bool* ok = nullptr)
{
    if (ok) *ok = true;
    if (string.compare("inf", Qt::CaseInsensitive) == 0) return CORBO_INF_INT;
    if (string.compare("-inf", Qt::CaseInsensitive) == 0) return -CORBO_INF_INT;
    return corbo::util::clamp(string.toInt(ok), -CORBO_INF_INT, CORBO_INF_INT);
}

inline bool qstring_to_bool(const QString& string, bool* ok = nullptr)
{
    if (ok) *ok = true;
    if (string.compare("true", Qt::CaseInsensitive) == 0) return true;
    if (string.compare("false", Qt::CaseInsensitive) == 0) return false;
    if (string.compare("yes", Qt::CaseInsensitive) == 0) return true;
    if (string.compare("no", Qt::CaseInsensitive) == 0) return false;
    if (string.compare("y", Qt::CaseInsensitive) == 0) return true;
    if (string.compare("n", Qt::CaseInsensitive) == 0) return false;
    return (bool)string.toInt(ok);
}

inline void qstring_to_container(const QString& string, std::vector<double>& values, bool* ok = nullptr)
{
    QRegExp rx("[, ]");  // match comma or space
    QStringList list = string.split(rx, QString::SkipEmptyParts);
    values.clear();
    bool ok_i   = true;
    if (ok) *ok = true;
    for (int i = 0; i < list.size(); ++i)
    {
        values.push_back(qstring_to_double(list.at(i), &ok_i));
        if (ok) *ok = *ok && ok_i;
    }
}

inline void qstring_to_container(const QString& string, std::vector<int>& values, bool* ok = nullptr)
{
    QRegExp rx("[, ]");  // match comma or space
    QStringList list = string.split(rx, QString::SkipEmptyParts);
    values.clear();
    bool ok_i   = true;
    if (ok) *ok = true;
    for (int i = 0; i < list.size(); ++i)
    {
        values.push_back(qstring_to_int(list.at(i), &ok_i));
        if (ok) *ok = *ok && ok_i;
    }
}

inline void qstring_to_container(const QString& string, std::vector<bool>& values, bool* ok = nullptr)
{
    QRegExp rx("[, ]");  // match comma or space
    QStringList list = string.split(rx, QString::SkipEmptyParts);
    values.clear();
    bool ok_i   = true;
    if (ok) *ok = true;
    for (int i = 0; i < list.size(); ++i)
    {
        values.push_back(qstring_to_bool(list.at(i), &ok_i));
        if (ok) *ok = *ok && ok_i;
    }
}

inline QString double_to_qstring(double value)
{
    if (value >= CORBO_INF_DBL) return "inf";
    if (value <= -CORBO_INF_DBL) return "-inf";
    if (value == M_PI) return "pi";
    if (value == -M_PI) return "-pi";
    if (value == M_2_PI) return "2pi";
    if (value == -M_2_PI) return "-2pi";
    if (value == M_PI_2) return "pi/2";
    if (value == -M_PI_2) return "-pi/2";
    if (value == M_PI_4) return "pi/4";
    if (value == -M_PI_4) return "-pi/4";
    return QString::number(value);
}

inline QString int_to_qstring(int value)
{
    if (value >= CORBO_INF_INT) return "inf";
    if (value <= -CORBO_INF_INT) return "-inf";
    return QString::number(value);
}

template <typename Iterator>
QString double_container_to_qstring(Iterator first, Iterator last)
{
    static_assert(std::is_same<typename Iterator::value_type, double>::value, "container of value type 'double' required");

    if (first == last) return QString("");

    QString values_string;
    for (Iterator it = first; it != last;)
    {
        values_string += double_to_qstring(*it);
        ++it;
        if (it != last) values_string += ", ";
    }
    return values_string;
}

template <typename Iterator>
QString int_container_to_qstring(Iterator first, Iterator last)
{
    static_assert(std::is_same<typename Iterator::value_type, int>::value, "container of value type 'int' required");

    if (first == last) return QString();

    QString values_string;
    for (Iterator it = first; it != last;)
    {
        values_string += int_to_qstring(*it);
        ++it;
        if (it != last) values_string += ", ";
    }
    return values_string;
}

template <typename Iterator>
QString bool_container_to_qstring(Iterator first, Iterator last)
{
    static_assert(std::is_same<typename Iterator::value_type, bool>::value, "container of value type 'bool' required");

    if (first == last) return QString();

    QString values_string;
    for (Iterator it = first; it != last;)
    {
        values_string += QString::number(*it);  // number implicitly casts to 0 or 1
        ++it;
        if (it != last) values_string += ", ";
    }
    return values_string;
}

inline QFrame* create_horizontal_line()
{
    QFrame* hline = new QFrame;
    hline->setFrameShape(QFrame::HLine);
    hline->setFrameShadow(QFrame::Sunken);
    hline->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);
    return hline;
}

}  // namespace util

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_UTILITIES_H_

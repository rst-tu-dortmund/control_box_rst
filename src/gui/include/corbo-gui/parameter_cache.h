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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_CACHE_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_CACHE_H_

#include <google/protobuf/message.h>

#include <list>
#include <memory>
#include <string>
#include <unordered_map>

namespace corbo {
namespace gui {

class ParameterCache
{
 public:
    using CacheMap = std::unordered_map<std::string, std::unique_ptr<google::protobuf::Message>>;

    void toCache(const std::string& key, const google::protobuf::Message& message);

    std::unique_ptr<google::protobuf::Message> fromCache(const std::string& key) const;

    bool hasCache(const std::string& key) const;

    void erase(const std::string& key) { _params.erase(key); }

    static std::string nestedNameListToKey(const std::list<std::string>& nested_name);

 private:
    CacheMap _params;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_CACHE_H_

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

#include <corbo-gui/parameter_cache.h>

#include <list>
#include <memory>
#include <string>

namespace corbo {
namespace gui {

void ParameterCache::toCache(const std::string& key, const google::protobuf::Message& message)
{
    std::unique_ptr<google::protobuf::Message>& param = _params[key];
    param                                             = std::unique_ptr<google::protobuf::Message>(message.New());
    param->CopyFrom(message);
}

std::unique_ptr<google::protobuf::Message> ParameterCache::fromCache(const std::string& key) const
{
    auto param_it = _params.find(key);
    if (param_it != _params.end())
    {
        std::unique_ptr<google::protobuf::Message> param(param_it->second->New());
        param->CopyFrom(*param_it->second);
        return param;
    }
    return std::unique_ptr<google::protobuf::Message>();
}

bool ParameterCache::hasCache(const std::string& key) const { return _params.find(key) != _params.end(); }

std::string ParameterCache::nestedNameListToKey(const std::list<std::string>& nested_name)
{
    std::string key;
    auto it       = nested_name.begin();
    auto prev_end = std::prev(nested_name.end());
    for (; it != nested_name.end(); ++it)
    {
        key += *it;
        if (it != prev_end) key += "/";
    }
    return key;
}

}  // namespace gui
}  // namespace corbo

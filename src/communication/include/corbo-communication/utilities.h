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

#ifndef SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_UTILITIES_H_
#define SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_UTILITIES_H_

#include <corbo-core/factory.h>

#ifdef MESSAGE_SUPPORT
#include <google/protobuf/message.h>
#endif

#include <memory>
#include <string>

namespace corbo {

namespace util {

#ifdef MESSAGE_SUPPORT
/**
 * @brief Helper function for other classes to get the typename of an active oneof field
 * @param[in]  message           Protobuf message that contains the related oneof field
 * @param[in]  oneof_name        Name of the oneof
 * @param[out] item_name         Name of the selected item in the oneof
 * @param[in]  include_namespace Should the namespace be included in \c item_name?
 * @returns true if \c item_name could be retrieved, false otherwise.
 */
bool get_oneof_field_type(const google::protobuf::Message& message, const std::string& oneof_name, std::string& item_name,
                          bool include_namespace = false);

/**
 * @brief Helper function for other classes to get the typename of an active oneof field
 *
 * In comparison to get_oneof_field_type() this version check if the given top-level
 * one-of has also a single nested (isolated) one-of.
 * E.g. consider the following scenario:
 * You want to define a couple of Plant implementations for different robots.
 * Instead of adding all robots to the PlantMessage, you might want to create a nested
 * message: RobotPlantMessage.
 * Now should RobotPlantMessage contain a (single) one-of type containing all robot plants.
 *
 * @param[in]  top_message       Top-level protobuf message which is recursively checked
 * @param[in]  oneof_name        Top-level one-of descriptor which is checked for a single child recursively
 * @param[out] item_name         Name of the selected item in the oneof
 * @param[in]  include_namespace Should the namespace be included in \c item_name?
 * @param[in]  depth             Specify how often an isolated one of should be expanded
 * @returns true if \c item_name could be retrieved, false otherwise.
 */
bool get_oneof_field_type_expand_isolated(const google::protobuf::Message& top_message, const google::protobuf::OneofDescriptor* top_oneof,
                                          std::string& item_name, bool include_namespace, int max_depth = 100);

/**
 * @brief Helper function for other classes to get the typename of an active oneof field
 *
 * Note, this is just an overload (refer to the other class for a detailed description).
 * In this case, the top-level oneof is specified by it's field name.
 *
 * @param[in]  top_message       Top-level protobuf message which is recursively checked
 * @param[in]  top_oneof_name    Field name of the top-level one-of descriptor which is checked for a single child recursively
 * @param[out] item_name         Name of the selected item in the oneof
 * @param[in]  include_namespace Should the namespace be included in \c item_name?
 * @param[in]  depth             Specify how often an isolated one of should be expanded
 * @returns true if \c item_name could be retrieved, false otherwise.
 */
bool get_oneof_field_type_expand_isolated(const google::protobuf::Message& message, const std::string& top_oneof_name, std::string& item_name,
                                          bool include_namespace, int max_depth = 100);

/**
 * @brief Helper function that creates classes from one-of fields using their factory
 *
 * This method tries to find the type of the selected oneof class (or a nested type if \c max_search_depth > 0).
 * The class which should be created must provide the following members:
 *
 *     static Factory<Base>& getFactory();
 *     void fromMessage(const Message&  message, std::stringstream* issues);
 *
 * This class also initalizes the parameters of the created class by calling fromMessage().
 *
 * @remarks  It might happen, that the algorithm expands to deep and hence does not get the correct field type;
 *                   in that case you can set the depth with \c max_search_depth
 * @todo        future work could try to trigger the factory recursively
 *
 * @param[in]  message                      Top-level protobuf message which is checked recursively
 * @param[in]  issues                           Issues during parameter initialization might be forwareded to this stringstream
 * @param[in]  top_oneof_name           Field name of the top-level one-of descriptor which is checked for a single child recursively
 * @param[in]  max_search_depth       Specify how often an isolated one of should be expanded
 * @param[in]  include_namespace      Should the namespace be included in \c item_name?

 * @returns the create item if everything is ok, a nullptr-cosntructed item otherwise.
 *
 * @tparam Base         Type of the base class of the obejct to be created by the factory
 * @tparam Message   Protobuf message type (this template type can be deduced from the argument list)
 */
template <typename Base, typename Message>
std::shared_ptr<Base> create_oneof_type_from_factory(const Message& message, std::stringstream* issues = nullptr,
                                                     const std::string& top_oneof_name = "", int max_search_depth = 100,
                                                     bool include_namespace = false)
{
    std::string item_name;
    if (!get_oneof_field_type_expand_isolated(message, top_oneof_name, item_name, include_namespace, max_search_depth)) return {};
    if (item_name.empty()) return {};
    std::shared_ptr<Base> base_ptr = create_from_factory<Base>(item_name);
    if (!base_ptr) return {};
    base_ptr->fromMessage(message, issues);
    return base_ptr;
}

#endif

}  // namespace util

}  // namespace corbo

#endif  // SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_UTILITIES_H_

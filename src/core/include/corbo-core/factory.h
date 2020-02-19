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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_FACTORY_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_FACTORY_H_

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace corbo {

/**
 * @brief Generic factory object
 *
 * @ingroup core
 *
 * A factory is used to create objects based on a known identifier (which is a string in this case).
 * This facilitates the dynamic creation of specific interface implementations without including their
 * actual declaration. Registration of objects is deferred from the library to the linker stage of the
 * executable.
 * This factory is implemented using the Singleton design pattern in order to facilitate the
 * dynamic creation of objects wherever required without the need of sharing the factory
 * among multiple classes and modules.
 *
 * A macro for a convenient object registration is provided which can be invoked in the header file of the particular
 * interface implementation. Hence, the source of the execuable only needs to include the particular header
 * (containing the macro).
 *
 * Note, interfaces must define a method of prototype:
 * @code
 *  virtual std::shared_ptr<Base> getInstance() = 0;
 * @endcode
 * This method must be implemented in all subclasses in order to return a valid instance of the related subclass.
 * Furthermore, the default constructure most not be deleted.
 *
 * @remark It is recommended to keep the class' default constructors as simple and efficient as possible
 *         since a default-constructed instance for each registered object is kept in the factory.
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
template <typename Base>
class Factory
{
 public:
    //!< Retrieve static instance of the factory
    static Factory& instance()
    {
        static Factory factory;
        return factory;
    }

    /**
     * @brief Check if a specified object (name) has been registered
     * @param[in] name    Unique object identifier
     * @returns true if object is registered, false otherwise
     */
    bool hasObject(const std::string& name) { return _object_map.count(name) > 0; }

    /**
     * @brief Create a shared instance of the desired object
     * @param[in] name           Unique object identifier
     * @param[in] print_error    If true, an error is printed if the object is unknown
     * @returns shared pointer to the derived class if \c name is known, an empty shared_ptr otherwise
     * @tparam Derived Specify proper subclass
     */
    template <typename Derived = Base>
    std::shared_ptr<Derived> create(const std::string& name, bool print_error = true) const
    {
        static_assert(std::is_base_of<Base, Derived>::value, "Factory::create(): specified dervied class type is an actual subclass of base.");

        auto obj = _object_map.find(name);
        if (obj == _object_map.end())
        {
            PRINT_ERROR_COND(print_error, "Factory<" << typeid(Base).name() << ">::create():: unkown object '" << name << "'.");
            return std::shared_ptr<Derived>();
        }

        std::shared_ptr<Derived> derived_ptr = std::dynamic_pointer_cast<Derived>(obj->second->getInstance());
        PRINT_ERROR_COND(!derived_ptr && print_error, "Factory::create():: cannot cast base object " << typeid(Base).name() << " to derived object "
                                                                                                     << typeid(Derived).name() << ".");
        return derived_ptr;
    }

    /**
     * @brief Register a new object
     * @param[in] name           Unique object identifier
     * @param[in] object_type    Default constructed object of the given type
     * @returns true if successful
     */
    bool registerObject(const std::string& name, std::shared_ptr<Base> object_type)
    {
        _object_map[name] = object_type;
        return true;
    }

    //! Print registered object names to console
    void printObjects() const
    {
        for (auto it = _object_map.begin(); it != _object_map.end(); ++it) PRINT_INFO("* " << it->first);
    }
    
    const std::unordered_map<std::string, std::shared_ptr<Base>>& getObjectMap() const { return _object_map; }

 private:
    //! map of identifiers and their corresponding objects (default-constructed)
    std::unordered_map<std::string, std::shared_ptr<Base>> _object_map;
};

#define FACTORY_REGISTER_OBJECT_ID(type, base, id)               \
    static const bool corbo_CAT(corbo_CAT(type, __regged), id) = \
        Factory<base>::instance().registerObject(corbo_STRINGIZE(type), std::make_shared<type>());

#define FACTORY_REGISTER_OBJECT(type, base) FACTORY_REGISTER_OBJECT_ID(type, base, 0)

/**
 * @brief Helper function to create new (derived) objects from factory
 *
 * This functino requries the Objects class to provide a getFactory() method.
 * @param[in] name    Unique object identifier
 * @returns shared pointer to the derived class if \c name is known, an empty shared_ptr otherwise
 * @tparam Object     The object type itself or base type to be created
 */
template <typename Object>
std::shared_ptr<Object> create_from_factory(const std::string& name)
{
    return Object::getFactory().template create<Object>(name);
}

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_FACTORY_H_

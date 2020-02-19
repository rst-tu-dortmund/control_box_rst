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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_CONSOLE_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_CONSOLE_H_

#include <corbo-core/text_style.h>
#include <string>

// === Colors ===

#define WARNING_COLOR corbo::TextColorCode::FG_LIGHT_YELLOW
#define ERROR_COLOR corbo::TextColorCode::FG_LIGHT_RED
#define DEFAULT_COLOR corbo::TextColorCode::FG_DEFAULT

// ==== Function name ===

#if defined(__GNUC__)
#define corbo_FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#define corbo_FUNCTION_NAME __FUNCSIG__
#else
#define corbo_FUNCTION_NAME __FUNCTION__
#endif

#define corbo_FUNCTION_NAME_FORMATTED "[" << corbo_FUNCTION_NAME << "] "

// === Debug messages ===

#if defined(NDEBUG) || defined(DISABLE_IO)

#define PRINT_DEBUG(msg)
#define PRINT_DEBUG_ONCE(msg)
#define PRINT_DEBUG_COND(cond, msg)
#define PRINT_DEBUG_COND_ONCE(cond, msg)

#else
//! Print \c msg-stream only if project is compiled in Debug-mode
#define PRINT_DEBUG(msg) std::cout << "Debug: " << msg << std::endl;

//! Print \c msg-stream only once and only if project is compiled in Debug-mode
#define PRINT_DEBUG_ONCE(msg)                           \
    {                                                   \
        static const auto debugOnce = [&] {             \
            std::cout << "Debug: " << msg << std::endl; \
            return true;                                \
        }();                                            \
        (void)debugOnce;                                \
    }  // void cast: avoid compiler warnings since it is unused later

//! Print \c msg-stream only if \c cond == \c true and only if project is compiled in Debug-mode
#define PRINT_DEBUG_COND(cond, msg) \
    if (cond) std::cout << "Debug: " << msg << std::endl;

//! Print \c msg-stream only if \c cond == \c true, only once and only if project is compiled in Debug-mode
#define PRINT_DEBUG_COND_ONCE(cond, msg)                          \
    {                                                             \
        static const auto debugOnce = [&] {                       \
            if (cond) std::cout << "Debug: " << msg << std::endl; \
            return true;                                          \
        }();                                                      \
        (void)debugOnce;                                          \
    }

#endif

#ifdef DISABLE_IO

#define PRINT_INFO(msg)
#define PRINT_INFO_ONCE(msg)
#define PRINT_INFO_COND(cond, msg)
#define PRINT_INFO_COND_ONCE(cond, msg)

#define PRINT_WARNING(msg)
#define PRINT_WARNING_ONCE(msg)
#define PRINT_WARNING_COND(cond, msg)
#define PRINT_WARNING_COND_ONCE(cond, msg)

#define PRINT_ERROR(msg)
#define PRINT_ERROR_ONCE(msg)
#define PRINT_ERROR_COND(cond, msg)
#define PRINT_ERROR_COND_ONCE(cond, msg)

#define PRINT_FATAL(msg)

#define INPUT_STREAM(variable, default_val) variable = default_val;
#define INPUT_STREAM_MSG(msg, variable, default_val) INPUT_STREAM(variable, default_val);
#define INPUT_STREAM_DEFAULT(variable, default_val) variable = default_val;
#define INPUT_STREAM_DEFAULT_MSG(msg, variable, default_val) INPUT_STREAM_DEFAULT(variable, default_val);

#else

#include <iostream>

//! Print \c msg-stream
#define PRINT_INFO(msg) std::cout << "Info: " << msg << std::endl
// #define PRINT_INFO_ARGS(msg, ...) std::printf(msg, ##__VA_ARGS__);

//! Print \c msg-stream only once
#define PRINT_INFO_ONCE(msg)                           \
    {                                                  \
        static const auto infoOnce = [&] {             \
            std::cout << "Info: " << msg << std::endl; \
            return true;                               \
        }();                                           \
        (void)infoOnce;                                \
    }  // void cast: avoid compiler warnings since it is unused later

//! Print \c msg-stream only if \c cond == \c true
#define PRINT_INFO_COND(cond, msg) \
    if (cond) std::cout << "Info: " << msg << std::endl

//! Print \c msg-stream only if \c cond == \c true, only once
#define PRINT_INFO_COND_ONCE(cond, msg)                          \
    {                                                            \
        static const auto infoOnce = [&] {                       \
            if (cond) std::cout << "Info: " << msg << std::endl; \
            return true;                                         \
        }();                                                     \
        (void)infoOnce;                                          \
    }

//! Print \c msg-stream
#define PRINT_WARNING(msg) std::cout << WARNING_COLOR << "Warning: " << msg << DEFAULT_COLOR << std::endl
// #define PRINT_INFO_ARGS(msg, ...) std::printf(msg, ##__VA_ARGS__);

//! Print \c msg-stream only once
#define PRINT_WARNING_ONCE(msg)                                                             \
    {                                                                                       \
        static const auto infoOnce = [&] {                                                  \
            std::cout << WARNING_COLOR << "Warning: " << msg << DEFAULT_COLOR << std::endl; \
            return true;                                                                    \
        }();                                                                                \
        (void)infoOnce;                                                                     \
    }  // void cast: avoid compiler warnings since it is unused later

//! Print \c msg-stream only if \c cond == \c true
#define PRINT_WARNING_COND(cond, msg) \
    if (cond) std::cout << WARNING_COLOR << "Warning: " << msg << DEFAULT_COLOR << std::endl

//! Print \c msg-stream only if \c cond == \c true, only once
#define PRINT_WARNING_COND_ONCE(cond, msg)                                                            \
    {                                                                                                 \
        static const auto warningOnce = [&] {                                                         \
            if (cond) std::cout << WARNING_COLOR << "Warning: " << msg << DEFAULT_COLOR << std::endl; \
            return true;                                                                              \
        }();                                                                                          \
        (void)warningOnce;                                                                            \
    }

//! Print \c msg-stream as error msg
#define PRINT_ERROR(msg) std::cerr << ERROR_COLOR << "Error: " << msg << DEFAULT_COLOR << std::endl

//! Print \c msg-stream only once
#define PRINT_ERROR_ONCE(msg)                                                           \
    {                                                                                   \
        static const auto infoOnce = [&] {                                              \
            std::cerr << ERROR_COLOR << "Error: " << msg << DEFAULT_COLOR << std::endl; \
            return true;                                                                \
        }();                                                                            \
        (void)infoOnce;                                                                 \
    }  // void cast: avoid compiler warnings since it is unused later

//! Print \c msg-stream only if \c cond == \c true
#define PRINT_ERROR_COND(cond, msg) \
    if (cond) std::cerr << ERROR_COLOR << "Error: " << msg << DEFAULT_COLOR << std::endl

//! Print \c msg-stream only if \c cond == \c true, only once
#define PRINT_ERROR_COND_ONCE(cond, msg)                                                          \
    {                                                                                             \
        static const auto errorOnce = [&] {                                                       \
            if (cond) std::cerr << ERROR_COLOR << "Error: " << msg << DEFAULT_COLOR << std::endl; \
            return true;                                                                          \
        }();                                                                                      \
        (void)errorOnce;                                                                          \
    }

//! Print \c msg-stream as fatal error msg and exit program execution
#define PRINT_FATAL(msg)                                                                                                    \
    {                                                                                                                       \
        std::cerr << ERROR_COLOR << "Fatal error: " << msg << " Stopping program execution." << DEFAULT_COLOR << std::endl; \
        exit(1);                                                                                                            \
    }

//! Print message in debug mode with warning-color-code
#define PRINT_DEBUG_WARN(msg) PRINT_DEBUG(WARNING_COLOR << msg << DEFAULT_COLOR)

#define INPUT_STREAM(variable, default_val) std::cin >> variable

#define INPUT_STREAM_MSG(msg, variable, default_val) \
    std::cout << "Interaction: " << msg << ". ";     \
    std::cin >> variable

// Never use unsigned variables
#define INPUT_STREAM_DEFAULT(variable, default_val)                                               \
    {                                                                                             \
        std::cout << "[" << default_val << "] ";                                                  \
        std::string name;                                                                         \
        std::getline(std::cin, name);                                                             \
        if (name.empty())                                                                         \
            variable = default_val;                                                               \
        else                                                                                      \
        {                                                                                         \
            try                                                                                   \
            {                                                                                     \
                variable = static_cast<decltype(variable)>(std::stod(name));                      \
            }                                                                                     \
            catch (const std::invalid_argument& ia)                                               \
            {                                                                                     \
                std::cerr << "Invalid argument (no number found). Using default..." << std::endl; \
                variable = default_val;                                                           \
            }                                                                                     \
        }                                                                                         \
    }

// Never use unsigned variables
#define INPUT_STREAM_DEFAULT_MSG(msg, variable, default_val) \
    std::cout << "Interaction: " << msg << ". ";             \
    INPUT_STREAM_DEFAULT(variable, default_val)

#endif

// Extend print macros by named version (which also print the function signature)
#define PRINT_DEBUG_NAMED(msg) PRINT_DEBUG(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_DEBUG_ONCE_NAMED(msg) PRINT_DEBUG_ONCE(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_DEBUG_COND_NAMED(cond, msg) PRINT_DEBUG_COND(cond, corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_DEBUG_COND_ONCE_NAMED(cond, msg) PRINT_DEBUG_COND_ONCE(cond, corbo_FUNCTION_NAME_FORMATTED << msg)

#define PRINT_INFO_NAMED(msg) PRINT_INFO(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_INFO_ONCE_NAMED(msg) PRINT_INFO_ONCE(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_INFO_COND_NAMED(cond, msg) PRINT_INFO_COND(cond, corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_INFO_COND_ONCE_NAMED(cond, msg) PRINT_INFO_COND_ONCE(cond, corbo_FUNCTION_NAME_FORMATTED << msg)

#define PRINT_WARNING_NAMED(msg) PRINT_WARNING(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_WARNING_ONCE_NAMED(msg) PRINT_WARNING_ONCE(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_WARNING_COND_NAMED(cond, msg) PRINT_WARNING_COND(cond, corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_WARNING_COND_ONCE_NAMED(cond, msg) PRINT_WARNING_COND_ONCE(cond, corbo_FUNCTION_NAME_FORMATTED << msg)

#define PRINT_ERROR_NAMED(msg) PRINT_ERROR(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_ERROR_ONCE_NAMED(msg) PRINT_ERROR_ONCE(corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_ERROR_COND_NAMED(cond, msg) PRINT_ERROR_COND(cond, corbo_FUNCTION_NAME_FORMATTED << msg)
#define PRINT_ERROR_COND_ONCE_NAMED(cond, msg) PRINT_ERROR_COND_ONCE(cond, corbo_FUNCTION_NAME_FORMATTED << msg)

#define PRINT_FATAL_NAMED(msg) PRINT_FATAL(corbo_FUNCTION_NAME_FORMATTED << msg)

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_CONSOLE_H_

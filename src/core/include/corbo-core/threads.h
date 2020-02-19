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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_THREADS_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_THREADS_H_

#include <thread>

#if !defined(DISABLE_IO)
#include <cstring>
#include <iostream>
#endif

#ifdef __unix__
#include <pthread.h>  // required for set_thread_scheduling
#endif

namespace corbo {

namespace threads {

/**
 * @brief Change the scheduling priority of a thread
 * @param priority  Priority value between 1 and 99
 * @returns true if the scheduling options has been changed, false otherwise.
 */
inline bool set_thread_scheduling(std::thread& th, int priority)  // priority 1-99
{
#ifdef __unix__
    sched_param sch_params;
    sch_params.sched_priority = priority;
    if (pthread_setschedparam(th.native_handle(), SCHED_RR, &sch_params))
    {
#if !defined(DISABLE_IO)
        std::cerr << "Failed to set Thread scheduling : " << std::strerror(errno) << std::endl;
#endif
        return false;
    }
    return true;
#else
    return false;
#endif
}

}  // namespace threads

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_THREADS_H_

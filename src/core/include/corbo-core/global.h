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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_GLOBAL_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_GLOBAL_H_

namespace corbo {

/**
 * @brief global method to check whether to proceed or cancel the current action
 *
 * Time consuming parts of the program should check whether the current state is ok.
 * If not, they should return immediately. This allows other processes to interrupt
 * resepctively preempt long running processes.
 *
 * E.g. the execution of closed-loop control tasks might be configured for a long time duration,
 * but the user wants to stop by setting setOk(bool ok) to false. The tasks can check for
 * ok() frequently and abort its execution in this case.
 *
 * @return true if everything is still ok, false if an abort of the current action is requested.
 */
bool ok();

/**
 * @brief Change the current program state to ok or request the termination of the current action
 * @details Refer to ok() for more details.
 * @param ok  set to false, if current actions should be terminated (requires actions to check for ok()).
 */
void setOk(bool ok);

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_GLOBAL_H_

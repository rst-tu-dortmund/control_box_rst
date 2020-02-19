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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_TIME_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_TIME_H_

#include <chrono>
#include <thread>

#ifndef DISABLE_IO
#include <iostream>
#include <string>
#endif

#if _WIN32
#include <Windows.h>  // for wait / sleep
#undef min            // this macro conflicts with std::min
#undef max            // this macro conflicts with std::max
#else
#include <unistd.h>  // for wait /sleep
#endif

namespace corbo {

namespace internals {

template <typename T>
struct is_chrono_duration
{
    static constexpr bool value = false;
};

template <typename Rep, typename Period>
struct is_chrono_duration<std::chrono::duration<Rep, Period>>
{
    static constexpr bool value = true;
};

template <typename T>
struct is_chrono_timepoint
{
    static constexpr bool value = false;
};

template <typename Clock>
struct is_chrono_timepoint<std::chrono::time_point<Clock>>
{
    static constexpr bool value = true;
};
}  // end namespace internals

class Time;

/**
 * @brief Representation of time durations
 *
 * @ingroup core
 *
 * This object stores time durations and provides convenient methods
 * and operator overloads for common operations with time and durations.
 *
 * A duration object can be constructed from a value given in seconds
 * or in a generic way from a std::chrono<> type.
 *
 * In comparison to durations, actual time stamps are managed with the
 * Time class.
 *
 * The following arithmetic operators (in combination with Time objects)
 * are provided:
 *
 * Duration + Duration = Duration
 * Duration - Duration = Duration
 * Duration * double = Duration
 * Duration == Duration = bool
 * Duration > Duration == bool
 * Duration < Duration == bool
 *
 * This class is inspired and based on the ROS class ros::Duration
 * http://wiki.ros.org/roscpp/Overview/Time
 *
 * @see Time Rate
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Duration
{
    friend class Time;

 public:
    //! Default constructor
    Duration() {}
    //! Construct duration from integer (seconds)
    explicit Duration(int t) { fromSec(static_cast<double>(t)); }
    //! Construct duration from double (seconds)
    explicit Duration(double t) { fromSec(t); }

    /**
     * @brief Construct duration from std::chrono::duration object.
     *
     * E.g.:
     * @code
     *  Duration my_duration1 = Duration(std::chrono::milliseconds(10));
     *  Duration my_duration2 = Duration(std::chrono::microseconds(50));
     *  Duration my_duration3 = Duration(std::chrono::hours(1));
     * @endcode
     *
     * @param[in]  duration         object of template type ChronoDuration
     * @tparam     ChronoDuration   Template for a specific std::chrono::duration<> type
     */
    template <typename ChronoDuration>
    explicit Duration(const ChronoDuration& duration)
    {
        fromChrono(duration);
    }

    //! Return duration in seconds
    double toSec() const { return std::chrono::duration_cast<std::chrono::duration<double>>(_duration).count(); }
    //! Set duration from seconds (see Duration(double t))
    void fromSec(double t) { _duration = std::chrono::duration_cast<DurationType>(std::chrono::duration<double>(t)); }

    //! Set duration from std::chrono::duration object (see Duration(const ChronoDuration& duration))
    template <typename ChronoDuration>
    void fromChrono(const ChronoDuration& duration)
    {
        static_assert(internals::is_chrono_duration<ChronoDuration>::value,
                      "Template parameter ChronoDuration must be of type std::chrono::duration");
        _duration = std::chrono::duration_cast<DurationType>(duration);
    }

    /**
     * @brief Convert duration to std::chrono::duration object
     *
     * The exact chrono type might be specified as template argument,
     * e.g.:
     * @code
     *   auto dur = Duration(5).toChrono<std::chrono::duration<double>>(); // value in seconds, floating point
     *   auto dur = Duration(5).toChrono<std::chrono::microseconds>();
     *   auto dur = Duration(5).toChrono<std::chrono::milliseconds>();
     * @endcode
     * @returns desired std::chrono::duration<> instance
     */
    template <typename DurType = std::chrono::duration<double>>
    DurType toChrono()
    {
        return std::chrono::duration_cast<DurType>(_duration);
    }

    Time toTime(double basis_time = 0);

    /**
     * @brief Sleep (current thread) for the specified duration
     *
     * @code
     *   Duration(2).sleep();  // pause current thread for 2 seconds
     * @endcode
     */
    void sleep() { std::this_thread::sleep_for(_duration); }

    // operators
    Duration operator+(const Duration& rhs) const { return Duration(_duration + rhs._duration); }
    Duration operator-(const Duration& rhs) const { return Duration(_duration - rhs._duration); }
    Duration operator-() const { return Duration(-_duration); }
    Duration operator*(double scale) const { return Duration(toSec() * scale); }
    Duration& operator+=(const Duration& rhs)
    {
        _duration += rhs._duration;
        return *static_cast<Duration*>(this);
    }
    Duration& operator-=(const Duration& rhs)
    {
        _duration -= rhs._duration;
        return *static_cast<Duration*>(this);
    }
    Duration& operator*=(double scale)
    {
        fromSec(toSec() * scale);
        return *static_cast<Duration*>(this);
    }
    bool operator==(const Duration& rhs) const { return _duration == rhs._duration; }
    inline bool operator!=(const Duration& rhs) const { return !(*static_cast<const Duration*>(this) == rhs); }
    bool operator>(const Duration& rhs) const { return _duration > rhs._duration; }
    bool operator<(const Duration& rhs) const { return _duration < rhs._duration; }
    bool operator>=(const Duration& rhs) const { return _duration >= rhs._duration; }
    bool operator<=(const Duration& rhs) const { return _duration <= rhs._duration; }
#ifndef DISABLE_IO
    friend std::ostream& operator<<(std::ostream& os, const Duration& rhs)
    {
        os << rhs.toSec();
        return os;
    }
#endif

 private:
    using DurationType = std::chrono::nanoseconds;
    DurationType _duration;
};

/**
 * @brief Representation of time stamps
 *
 * @ingroup core
 *
 * This object stores time stamps and provides convenient methods
 * and operator overloads for common operations with time and durations.
 *
 * A time object can be constructed from a value given in seconds
 * or in a generic way from a std::chrono<>::time_point type.
 *
 * In comparison to durations, actual time stamps are managed with the
 * Time class. The internal high-resolution (system) clock is used
 * for precise time measurement.
 *
 * The following arithmetic operators (in combination with Time objects)
 * are provided:
 *
 * Time - Time = Duration
 * Time + Duration = Time
 * Time - Duration = Time
 * Time == Time = bool
 * Time > Time == bool
 * Time < Time == bool
 *
 * This class is inspired and based on the ROS class ros::Time
 * http://wiki.ros.org/roscpp/Overview/Time
 *
 * @see Duration Rate
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Time
{
 public:
    //! Default constructor
    Time() {}

    //! Construct time object from integer (seconds)
    explicit Time(int t) { fromSec(static_cast<double>(t)); }
    //! Construct time object from double (seconds)
    explicit Time(double t) { fromSec(t); }

    /**
     * @brief Construct time from std::chrono::time_point object.
     *
     * @param[in]  timepoint         object of template type ChronoTimePoint
     * @tparam     ChronoTimePoint   Template for a specific std::chrono::time_point<> type
     */
    template <typename ChronoTimePoint>
    explicit Time(const ChronoTimePoint& timepoint)
    {
        fromChrono(timepoint);
    }

    //! Retrieve current system time
    static Time now() { return Time(Clock::now()); }

    //! Sleep (current thread) until the specified timestamp
    static void sleepUntil(const Time& time) { std::this_thread::sleep_until(time._timepoint); }

    //! Cast time stamp to seconds
    double toSec() const { return std::chrono::duration_cast<std::chrono::duration<double>>(_timepoint.time_since_epoch()).count(); }
    //! Set time stamp from seconds
    void fromSec(double t) { _timepoint = TimePoint(std::chrono::duration_cast<TimePoint::duration>(std::chrono::duration<double>(t))); }

    //! Set time stamp from std::chrono::time_stamp type
    template <typename ChronoTimePoint>
    void fromChrono(const ChronoTimePoint& timepoint)
    {
        static_assert(internals::is_chrono_timepoint<ChronoTimePoint>::value,
                      "Template parameter ChronoTimePoint must be of type std::chrono::time_point");
        _timepoint = std::chrono::time_point_cast<Clock::duration>(timepoint);
    }

    // Operators
    Duration operator-(const Time& rhs) const { return Duration(_timepoint - rhs._timepoint); }
    Time operator+(const Duration& rhs) const { return Time(_timepoint + rhs._duration); }
    Time operator-(const Duration& rhs) const { return Time(_timepoint - rhs._duration); }
    Time& operator+=(const Duration& rhs)
    {
        _timepoint += rhs._duration;
        return *static_cast<Time*>(this);
    }
    Time& operator-=(const Duration& rhs)
    {
        _timepoint -= rhs._duration;
        return *static_cast<Time*>(this);
    }
    bool operator==(const Time& rhs) const { return _timepoint == rhs._timepoint; }
    inline bool operator!=(const Time& rhs) const { return !(*static_cast<const Time*>(this) == rhs); }
    bool operator>(const Time& rhs) const { return _timepoint > rhs._timepoint; }
    bool operator<(const Time& rhs) const { return _timepoint < rhs._timepoint; }
    bool operator>=(const Time& rhs) const { return _timepoint >= rhs._timepoint; }
    bool operator<=(const Time& rhs) const { return _timepoint <= rhs._timepoint; }
#ifndef DISABLE_IO
    friend std::ostream& operator<<(std::ostream& os, const Time& rhs)
    {
        os << rhs.toSec();
        return os;
    }
#endif

 private:
    //! Datatype for the clock
    using Clock = std::chrono::high_resolution_clock;
    //! Datatype for the timepoint
    using TimePoint = Clock::time_point;

    TimePoint _timepoint;
};

/**
 * @brief Rate objects can be used to run loops at a desired rate
 *
 * @ingroup core
 *
 * @code
 *   Rate rate(10) // The loop should run at 10 Hz
 *   while (1)
 *   {
 *      // perform some work which is ideally finished before 1/10 s
 *      // in order to remain synchronized
 *      rate.sleep(); // pause to synchronize with 10 Hz
 *   }
 * @endcode
 *
 * This class is inspired and based on the ROS class ros::Time
 * http://wiki.ros.org/roscpp/Overview/Time
 *
 * @see Time Duration
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Rate
{
 public:
    //! Construct rate object from frequency [1/s]
    explicit Rate(double frequency) : _cycle_time(1.0 / frequency) {}
    //! Cosntruct rate object from duration object (interval length)
    explicit Rate(const Duration& dt) : _cycle_time(dt) {}
    //! Sleep for the remaining duration to met the desired frequency (w.r.t previous sleep call)
    bool sleep();
    void reset() { _start = Time::now(); }

    //! Update cycle time without resetting start time
    void updateCycleTime(const Duration& dt) { _cycle_time = dt; }

    //! Return actual execution time of the last cycle
    Duration lastCycleTime() const { return _last_cycle_time; }

 private:
    Time _start = Time::now();
    Duration _cycle_time;
    Duration _last_cycle_time = Duration(0);
};

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_TIME_H_

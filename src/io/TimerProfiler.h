/*
 * TimerProfiler.h
 *
 *  Created on: Apr 9, 2017
 *      Author: Andrei Costinescu
 */

#ifndef SRC_IO_TIMERPROFILER_H_
#define SRC_IO_TIMERPROFILER_H_

#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "utils/Timer.h"

class XMLfileUnits;

/**
@class TimerProfiler
@brief Class for managing timers across the simulation

This class is a manager for all timers in the simulation.
There should only be one instance of this class across the simulation.
This class supports a hierarchical structure of the timers and can output the time of the timers in a more user-friendly manner
*/
class TimerProfiler {
public:

	enum class Displaymode {
		ALL,
		ACTIVE,
		NON_ZERO,
		NONE
	};

	/**
	@brief Constructor of TimerProfiler

	The constructor creates a "virtual" base timer (named "_baseTimer"), which is the top of the timer hierarchy.
	All other inserted timers are sub-timers of this _baseTimer.
	Also reads and initializes all timers in the timer config-file.
	*/
	TimerProfiler();

	/** @brief Read in XML configuration for TimerProfiler.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <programtimers>
	     <displaymode>all|active|non-zero|none</displaymode>
	   </programtimers>
	   \endcode
	 * Display mode explanation:
	 * - all: display all registered timer
	 * - active: display all active timers
	 * - non-zero: display all timers which have non zero time
	 * - none: do not display timers
	 */
	void readXML(XMLfileUnits& xmlconfig);

	void setDisplayMode(Displaymode mode) { _displayMode = mode; }
	Displaymode getDisplayMode() { return _displayMode; }

	/**
	@fn void registerTimer(std::string timerName, std::vector<std::string> parentTimerNames, Timer *timer=nullptr, bool activate=true)
	@brief Adds a timer in the container.
	@param timerName The name of the timer to be added
	@param parentTimerNames The names of the timers of which this timer is a sub-timer of
	@param timer A pointer to the timer to be added in the container or \a nullptr if this is just a "virtual" timer ("manager" of timers)
	@param activate Boolean variable indicating whether the timer should be active or not. All calls to an inactive timer will have no effect on it
	*/
	void registerTimer(std::string timerName, std::vector<std::string> parentTimerNames, Timer *timer=nullptr, bool activate=true);

	/**
	@fn Timer* getTimer(std::string timerName)
	@brief Gets a pointer to the timer with name "timerName".
	@param timerName The name of the requested timer
	@return A pointer to the requested timer or nullptr if the timer is not an actual timer
	*/
	Timer* getTimer(std::string timerName);

	/**
	@fn void activateTimer(std::string timerName)
	@brief Activates the timer "timerName".
	@param timerName The name of the timer to be activated
	*/
	void activateTimer(std::string timerName);

	/**
	@fn void deactivateTimer(std::string timerName)
	@brief Deactivates the timer "timerName".
	@param timerName The name of the timer to be deactivated
	*/
	void deactivateTimer(std::string timerName);

	/**
	@fn void setSyncTimer(std::string timerName, bool sync)
	@brief Sets sync mode of the timer "timerName".
	@param timerName The name of the timer whose sync mode is to be set
	@param sync Boolean variable indicating whether the timer should be synced or not
	*/
	void setSyncTimer(std::string timerName, bool sync);

	/**
	@fn void printTimer(std::string timerName, std::string outputPrefix="")
	@brief Prints time collected by timer "timerName".
	@param timerName The name of the timer to printed
	@param outputPrefix String prefix for the hierarchical print
	*/
	void print(std::string timerName, std::string outputPrefix="");

	/**
	@fn void printTimers(std::string startingTimerName=_baseTimerName, std::string outputPrefix="")
	@brief Prints time collected by timer "timerName" and the time from all of its descendants.
	@param timerName The name of the timer from where to start the print; defaults to the base timer
	@param outputPrefix String prefix for the hierarchical print
	*/
	void printTimers(std::string timerName=_baseTimerName, std::string outputPrefix="");

	/**
	@fn void start(std::string timerName)
	@brief Starts the timer "timerName".
	@param timerName The name of the timer to be started
	*/
	void start(std::string timerName);

	/**
	@fn void stopTimer(std::string timerName)
	@brief Stops the timer "timerName".
	@param timerName The name of the timer to be stopped
	*/
	void stop(std::string timerName);

	/**
	@fn void resetTimer(std::string timerName)
	@brief Resets the timer "timerName".
	@param timerName The name of the timer to be reset
	*/
	void reset(std::string timerName);

	/**
	@fn void resetTimers(std::string timerName=_baseTimerName)
	@brief Resets the timer "timerName" and all of its descendants.
	@param timerName The name of the timer to be deactivated; defaults to the base timer
	*/
	void resetTimers(std::string timerName=_baseTimerName);

	/**
	@fn void readInitialTimersFromFile(std::string fileName)
	@brief Sets up initial configuration of the timers during the simulation.

	Reads timers provided in the timer config-file and inserts them in the container, (de-)activates them and creates the hierarchical structure of the timers.

	@param fileName The name of the timer config-file
	*/
	void readInitialTimersFromFile(std::string fileName);

	/**
	@fn void setOutputString(std::string timerName, std::string outputString)
	@brief Sets the output string for the timer output.
	@param timerName The name of the timer whose output string is to be modified
	@param outputString The string to be displayed whenever the timer is printed
	*/
	void setOutputString(std::string timerName, std::string outputString);

	/**
	@fn void deactivateTimer(std::string timerName)
	@brief Gets the output string for the timer output.
	@param timerName The name of the timer whose output string is requested
	@return The string to be displayed whenever the timer is printed
	*/
	std::string getOutputString(std::string timerName);

	/**
	@fn double getTime(std::string timerName)
	@brief Gets the time collected by timer "timerName".
	@param timerName The name of the timer whose collected time is requested
	@return Time collected by timer "timerName"
	*/
	double getTime(std::string timerName);

	/**
	 @fn void incrementTimerTimestepCounter()
	 @brief Increment the timestep-counter.
	 */
	void incrementTimerTimestepCounter() {
		++_numElapsedIterations;
	}

	/**
	 @fn unsigned long getNumElapsedIterations() const
	 @brief get the number of iterations since the last resetting of all timers.
	 @return Number of iterations since last resetting of timers.
	 */
	unsigned long getNumElapsedIterations() const {
		return _numElapsedIterations;
	}

private:
	/**
	@var static const std::string _baseTimerName
	@brief The name of the base timer
	*/
	static const std::string _baseTimerName;

	/**
	@fn bool _checkTimer(std::string timerName, bool checkActive=true)
	@brief Checks validity of timer.
	Checks if the timer is actually a timer and not a "virtual" one. If the timer is a timer and if checkActive is true the function also checks if the timer is active.
	@param timerName The name of the timer to be checked
	@param checkActive Boolean indicating whether to also check if the timer is active or not
	@return Result of the checks
	*/
	bool _checkTimer(std::string timerName, bool checkActive=true);

	/**
	@fn void _debugMessage(std::string timerName)
	@brief Prints in Log::global_log->debug() information message about timer.
	@param timerName The name of the timer
	*/
	void _debugMessage(std::string timerName);

	/**
	@class _Timer
	@brief Private helper class containing a timer and its hierarchical information
	*/
	class _Timer{
		public:
			/**
			@brief Default constructor of _Timer

			A default constructor is necessary because the container creates for every new key an object initialised with the default constructor.
			*/
			_Timer(std::string timerName="", Timer* timer=nullptr,
					std::vector<std::string>childTimerNames={}, std::vector<std::string>parentTimerNames={},
					std::string outputString=""):
				_timer(timer), _childTimerNames(childTimerNames), _parentTimerNames(parentTimerNames),
				_outputString(outputString), _timerName(timerName) {}

			/**
			@var std::unique_ptr<Timer> _timer
			@brief Pointer to the timer or nullptr if it is a "virtual" timer
			*/
			std::unique_ptr<Timer> _timer;

			/**
			@var std::vector<std::string> _childTimerNames
			@brief Vector of names of all descendants of the current timer
			*/
			std::vector<std::string> _childTimerNames;

			/**
			@var std::vector<std::string> _parentTimerNames
			@brief Vector of names of all parents of the current timer
			*/
			std::vector<std::string> _parentTimerNames;

			/**
			@var std::string _outputString
			@brief String to be displayed whenever the timer is printed
			This should be a message providing context for the usage of the timer.
			*/
			std::string _outputString;

			/**
			@var std::string _timerName
			@brief The name of the current timer
			*/
			std::string _timerName;
	};

	/**
	@var std::map<std::string, _Timer> _timers;
	@brief The container of all timers and their hierarchical information
	A mapping between timer names and the actual timer information.
	If no key matches the ones already in the map a new element will be (default-)constructed and inserted in the container.
	*/
	std::map<std::string, _Timer> _timers;

	unsigned long _numElapsedIterations;

	Displaymode _displayMode;
};

#endif /* SRC_IO_TIMERPROFILER_H_ */

/*
 * estimateRemainingTimeCommand.h
 *
 *  Created on: Jul 1, 2010
 *      Author: hpcdjenz
 */
#ifdef STEEREO
#ifndef __ESTIMATEREMAININGTIMECOMMAND_H__
#define __ESTIMATEREMAININGTIMECOMMAND_H__

#include <steereo/steereoCommand.h>

class EstimateRemainingTimeCommand : public SteereoCommand {
public:
	EstimateRemainingTimeCommand();
	virtual ~EstimateRemainingTimeCommand();

	virtual ReturnType execute ();
	void setParameters (std::list<std::string> params);
	static SteereoCommand* generateNewInstance ();
	bool condition ();
};

#endif /* __ESTIMATEREMAININGTIMECOMMAND_H__ */
#endif /* STEEREO */

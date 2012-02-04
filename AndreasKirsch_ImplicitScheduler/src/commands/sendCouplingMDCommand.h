/*
 * SendCouplingInfoCommand.h
 *
 *  Created on: Apr 23, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#ifndef SENDCOUPLINGINFOCOMMAND_H_
#define SENDCOUPLINGINFOCOMMAND_H_

#include <steereoCommand.h>

class SendCouplingMDCommand: public SteereoCommand {
public:
	SendCouplingMDCommand();
	virtual ~SendCouplingMDCommand();

	virtual ReturnType execute ();
	void setParameters (std::list<std::string> params);

	static SteereoCommand* generateNewInstance ();

	/*static void addData (void* dataPtr) {data.push_back (dataPtr);};
	static void* getData (int index = 0) {return data[index];};
	static std::vector<void*>* getDataVectorPtr () {return &data;};

	static void addDataSize (int size) {dataSize.push_back (size);};
	static int getDataSize (int index = 0) {return dataSize[index];};
	static std::vector<int>* getDataSizeVectorPtr () {return &dataSize;};
*/

	bool condition ();
	void setStepInterval (int interval) {stepInterval = interval;};

private:
	static int startStep;
	static int borderToLook;
	static int outmin, outmax;

	//static std::vector<void*> data;
	//static std::vector<int> dataSize;

	int stepInterval;
};

#endif /* SENDCOUPLINGINFOCOMMAND_H_ */
#endif /* STEEREO */

/*
 * ObserverBase.h
 *
 *  Created on: 31.10.2016
 *      Author: mheinen
 */

#ifndef OBSERVERBASE_H_
#define OBSERVERBASE_H_

class ObserverBase
{
public:
	virtual void set(double dInterfaceMidLeft, double dInterfaceMidRight) = 0;
};

class SubjectBase
{
public:
	virtual void registerObserver(ObserverBase* observer) = 0;
	virtual void deregisterObserver(ObserverBase* observer) = 0;
	virtual void informObserver() = 0;
};

#endif /* OBSERVERBASE_H_ */

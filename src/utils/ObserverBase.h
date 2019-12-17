/*
 * ObserverBase.h
 *
 *  Created on: 31.10.2016
 *      Author: mheinen
 */

#ifndef OBSERVERBASE_H_
#define OBSERVERBASE_H_

class SubjectBase;
class ObserverBase
{
public:
	virtual ~ObserverBase() {}
	virtual void update(SubjectBase* subject) = 0;
};

class SubjectBase
{
public:
	virtual ~SubjectBase() {}
	virtual void registerObserver(ObserverBase* observer) = 0;
	virtual void unregisterObserver(ObserverBase* observer) = 0;
	virtual void informObserver() = 0;
};

#endif /* OBSERVERBASE_H_ */

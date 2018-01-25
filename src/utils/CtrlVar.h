/*
 * CtrlVar.h
 *
 *  Created on: 25.01.2018
 *      Author: mheinen
 */

#ifndef CTRLVAR_H_
#define CTRLVAR_H_

/** Template class for defining a variable triplet: 'target', 'actual' and 'spread'.
 *  This is intended for storing information of a observable like e.g. the temperature in a control region:
 *  its target value, actual value and the spread between those.
 *
 *  You can also combine it with the template class CommVar, to store local and global information:
 * \code{.cpp}
	CtrlVar<CommVar<double> > temperature;
   \endcode
 */
template<typename T>
class CtrlVar
{
public:
	T target;
	T actual;
	T spread;
};

#endif /* CTRLVAR_H_ */

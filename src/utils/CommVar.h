/*
 * CommVar.h
 *
 *  Created on: 25.01.2018
 *      Author: mheinen
 */

#ifndef COMMVAR_H_
#define COMMVAR_H_

/** Template class for defining a variable pair 'local' and 'global', that carry information of the local process
 *  owning only a subdomain of the system, and the global information that is reduced from all processes and
 *  their subdomains.
 */
template<typename T>
class CommVar
{
public:
	T local;
	T global;
};

#endif /* COMMVAR_H_ */

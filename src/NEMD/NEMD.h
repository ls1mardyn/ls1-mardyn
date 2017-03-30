/*
 * NEMD.h
 *
 *  Created on: 30.03.2017
 *      Author: mheinen
 */

#ifndef NEMD_H_
#define NEMD_H_

#include <cstdint>

enum NEMDflags : uint32_t
{
	CHANGE_COMPONENT_AC_TO_N2 = (1u << 0)
};

#endif

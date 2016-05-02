/*
 * CollectiveCommunicationSingleNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */

#pragma once

#include "CollectiveCommunication.h"

/**
 * CollectiveCommunicationSingleNonBlocking extends the CollectiveCommunication by an implementation of nonblocking collectives.
 * @author Steffen Seckler
 */
class CollectiveCommunicationSingleNonBlocking: public CollectiveCommunication {

public:
	/**
	 * Constructor
	 * @param key the key of the collective operation
	 */
	CollectiveCommunicationSingleNonBlocking(int key) {
		_key = key;
	}

	/**
	 * Destructor
	 */
	virtual ~CollectiveCommunicationSingleNonBlocking() {
		assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

private:
	int _key;
};

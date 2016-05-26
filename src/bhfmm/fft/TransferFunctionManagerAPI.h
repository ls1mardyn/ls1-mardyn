/*
 * TransferFunctionManagerAPI.h
 *
 *  Created on: Apr 18, 2016
 *      Author: gallardjm
 */

#ifndef TRANSFERFUNCTIONMANAGER_API_H_
#define TRANSFERFUNCTIONMANAGER_API_H_

class TransferFunctionManagerAPI {
  public:
    virtual ~TransferFunctionManagerAPI() {}
    /**
     * get a transfer function corresponding the the input vector
     * 
     * @param utils::Vector<double,3> r, vector between the source and the target
     * @return FFTDataContainer*, the transfer function
     */
    virtual FFTDataContainer* getTransferFunction(int x, int y, int z, double cell_size_x, double cell_size_y, double cell_size_z) = 0;  
};

#endif

/*
 * TransferFunctionManager.cpp
 *
 *  Created on: Feb 05, 2016
 *      Author: gallardjm
 */

#include "TransferFunctionManager.h"

TransferFunctionManager::TransferFunctionManager(int ord, FFTAccelerationAPI* FFTA, bool verbose)
      : _ord(ord), _verbose(verbose), _asked(0), _builded(0), _FFTAcceleration(FFTA) 
{}

TransferFunctionManager::~TransferFunctionManager() 
{
  if(_verbose)
    std::cout << "TransferFunctionManager stats: TF asked=" << _asked << " | TF builded=" << _builded << std::endl;
}


FFTDataContainer* TransferFunctionManager::getTransferFunction(int x, int y, int z, double cell_size_x, double cell_size_y, double cell_size_z)
{
  _asked++;
  _builded++;
  DummyExpansion temp_M(_ord); //generate an empty expansion
  temp_M.evaluate_M_at_r(((double)x)*cell_size_x,((double)y)*cell_size_y,((double)z)*cell_size_z); //set it at r (= make it a transferfunction expansion of vector r)
  _FFTAcceleration->FFT_initialize_TransferFunction(temp_M); //use the FFTAcceleration to set its FFTDatacontainer
  
  //deep copy and return the data container (the one in the expansion will be destroy by the expansion destruct)
  return temp_M._FFTData->copyContainer(); 
}

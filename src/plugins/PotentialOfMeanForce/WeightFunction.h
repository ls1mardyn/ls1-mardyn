/*
 * Created on Thu Jun 13 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Region.h"

class WeightFunction{
    public:
    
    /*
    Position receibed must be COM of particle*/
    double WeightValue(std::array<double,3>& pos, FPRegion& region);
};
/*
 * Created on Tue Apr 30 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include<string>
#include<vector>

class Subset{
    private:

    std::string name;
    std::vector<int> local_nodes;

    public:

    Subset(std::string Name, std::vector<int> nodes);

    void SetNodes(std::vector<int> list);
    std::vector<int>& GetNodes();
    void SetName(std::string name);
    std::string& GetName();
    
};
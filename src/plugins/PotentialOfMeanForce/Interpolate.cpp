#include "Interpolate.h"

void Interpolate::ReadInRDF(){
    std::string filename;

    filename = "rdf.txt";
    std::ifstream file{filename};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the rdf data file"<<std::endl;
    }
    double n1, n2;

    while(file >> n1 >> n2){
        r_nodes.push_back(n1);
        g_nodes.push_back(n2);
    }
}

std::vector<double>& Interpolate::GetGValues(){
    return this->g_nodes;
}

std::vector<double>& Interpolate::GetRValues(){
    return this->r_nodes;
}


double Interpolate::GetRDFAt(double r){
    
    //need to artificially return a one
    if(r > r_nodes[r_nodes.size()-1]){
        return 1.0;
    }

    double gr =0.0;

    //find between which 2 values
    int low, up;
    up = GetUpperLimit(r);
    if(up ==0){
        return 0.0;
    }
    low = up -1;
    //get indeces gr values
    //interpolate
    gr = LinearInterpolation(low,up,r);
    return gr;
}

double Interpolate::CentralFiniteDifference(double r){
    //need to artificially return a zero since profile is flat
    if(r > r_nodes[r_nodes.size()-1]){
        return 0.0;
    }

    int low, up;
    up = GetUpperLimit(r);
    if(up == 0){
        return 0.0;
    }

    low = up -1;

    double gb,ga,ra,rb;
    ga = g_nodes[low]; ra = r_nodes[low];
    gb = g_nodes[up]; rb = r_nodes[up];

    return (gb-ga)/(rb-ra);
}

double Interpolate::LinearInterpolation(int fa, int fb, double fx){
    double ga, gb, ra, rb;
    ga = g_nodes[fa]; ra = r_nodes[fa];
    gb = g_nodes[fb]; rb = r_nodes[fb];

    double gc =0.0;
    gc = ga + (gb-ga)/(rb-ra)*(fx-ra);

    return gc;

}

int Interpolate::GetUpperLimit(double r){
    
    int i=0;
    int i_max = r_nodes.size()-1;
    while(r>r_nodes[i]){
        i++;
        if(i==i_max) break;
    }
    return i;
}
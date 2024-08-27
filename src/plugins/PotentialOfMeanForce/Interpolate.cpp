#include "Interpolate.h"

Interpolate::Interpolate(double def_val):default_value{def_val}{

}

void Interpolate::SetXValues(std::vector<double>& v){
    this->x_values = v;
}

void Interpolate::SetYValues(std::vector<double>& v){
    this->y_values = v;
}

void Interpolate::ReadInRDF(){
    std::string filename;

    filename = "rdf.txt";
    std::ifstream file{filename};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the rdf data file"<<std::endl;
    }
    double n1, n2;

    while(file >> n1 >> n2){
        x_values.push_back(n1);
        y_values.push_back(n2);
    }
}

std::vector<double>& Interpolate::GetGValues(){
    return this->y_values;
}

std::vector<double>& Interpolate::GetRValues(){
    return this->x_values;
}


double Interpolate::GetRDFAt(double r){
    
    //need to artificially return a one
    if(r > x_values[x_values.size()-1]){
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
    if(r > x_values[x_values.size()-1]){
        return 0.0;
    }

    int low, up;
    up = GetUpperLimit(r);
    if(up == 0){
        return 0.0;
    }

    low = up -1;

    double gb,ga,ra,rb;
    ga = y_values[low]; ra = x_values[low];
    gb = y_values[up]; rb = x_values[up];

    return (gb-ga)/(rb-ra);
}

double Interpolate::LinearInterpolation(int fa, int fb, double fx){
    double ga, gb, ra, rb;
    ga = y_values[fa]; ra = x_values[fa];
    gb = y_values[fb]; rb = x_values[fb];

    double gc =0.0;
    gc = ga + (gb-ga)/(rb-ra)*(fx-ra);

    return gc;

}

int Interpolate::GetUpperLimit(double r){
    
    int i=0;
    int i_max = x_values.size()-1;
    while(r>x_values[i]){
        i++;
        if(i==i_max) break;
    }
    return i;
}
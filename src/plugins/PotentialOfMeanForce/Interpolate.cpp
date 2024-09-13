#include "Interpolate.h"

Interpolate::Interpolate(double def_val):default_value{def_val}{

}

void Interpolate::SetXValues(std::vector<double>& v){
    this->x_values = v;
}

void Interpolate::SetYValues(std::vector<double>& v){
    this->y_values = v;
}

std::vector<double>& Interpolate::GetXValues(){
    return this->y_values;
}

std::vector<double>& Interpolate::GetYValues(){
    return this->x_values;
}


double Interpolate::InterpolateAt(double r){
    
    //need to artificially return a one
    if(r > x_values[x_values.size()-1]){
        return default_value;
    }

    double y_at_r =0.0;

    //find between which 2 values
    int low, up;
    up = GetUpperLimit(r);
    if(up ==0){
        return 0.0;
    }
    low = up -1;
    //get indeces gr values
    //interpolate
    y_at_r = LinearInterpolation(low,up,r);

    if(!std::isfinite(y_at_r)){
        int i=0;
        while(std::isinf(y_values[i])){
            i++;
        }
        return y_values[i];
    }

    return y_at_r;
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
    double ratio =(gb-ga)/(rb-ra);

    if(!std::isfinite(ratio)){
        return 5.0;
    }

    return ratio;
}

double Interpolate::LinearInterpolation(int fa, int fb, double fx){
    double ya, yb, xa, xb;
    ya = y_values[fa]; xa = x_values[fa];
    yb = y_values[fb]; xb = x_values[fb];

    double gc =0.0;
    gc = ya + (yb-ya)/(xb-xa)*(fx-xa);

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
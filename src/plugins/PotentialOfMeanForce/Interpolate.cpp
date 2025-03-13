#include "Interpolate.h"

Interpolate::Interpolate(double def_val):default_value{def_val}{

}

Interpolate::Interpolate(double d, int sz):default_value{d}{
    x_values.resize(sz);
    y_values.resize(sz);
}

void Interpolate::ResizeVectors(int size){
    x_values.resize(size);
    y_values.resize(size);
}

void Interpolate::SetyXYValues(std::vector<double>& x, std::vector<double>& y){
    x_values = x;
    y_values = y;
}

void Interpolate::SetXValues(std::vector<double>& v){
    this->x_values = v;
}

void Interpolate::SetYValues(std::vector<double>& v){
    this->y_values = v;
}

std::vector<double>& Interpolate::GetXValues(){
    return this->x_values;
}

std::vector<double>& Interpolate::GetYValues(){
    return this->y_values;
}


double Interpolate::InterpolateAt(double r){
    
    // Return default value if larger than largest x value
    if(r > x_values[x_values.size()-1]){
        return default_value;
    }

    double y_at_r =0.0;

    //find between which 2 values
    int low, up;
    up = GetUpperLimit(r);

    // if(up <= first_valid && check_for_isfinite){
    //     return y_values[first_valid];
    // }

    if(up ==0){
        return y_values[0];
    }

    low = up -1;
    //get indeces gr values
    //interpolate
    y_at_r = LinearInterpolation(low,up,r);

    return y_at_r;
}

void Interpolate::AddVector(std::vector<double>& v){
    if(v.size()!=y_values.size()){
        Log::global_log->error()<<"[Interpolate]Not same size add"<<std::endl;
    }

    for(int i=0;i<v.size();++i){
        if(!std::isfinite(v[i])) continue;
            y_values[i] += v[i];
    }

}

double Interpolate::CentralFiniteDifference(double r){
    // Return zero since we are at the limit and function is flat
    if(r > _simulation.getcutoffRadius()){
        Log::global_log->warning()<<"Why are we checking beyond r_c="<<x_values[x_values.size()-1]<<"\t r="<<r<<std::endl;
        return 0.0;
    }

    int low, up;
    up = GetUpperLimit(r);

    low = up -1;

    double gb,ga,ra,rb;
    ga = y_values[low]; ra = x_values[low];
    gb = y_values[up]; rb = x_values[up];
    double ratio =(gb-ga)/(2.0*(rb-ra));


    return ratio;
}

double Interpolate::LinearInterpolation(int low, int up, double fx){
    double ya, yb, xa, xb;
    ya = y_values[low]; xa = x_values[low];
    yb = y_values[up]; xb = x_values[up];

    return (ya + ((fx-xa)*(yb-ya)/(xb-xa)));

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

void Interpolate::LinearExtrapolation(){
    // find first valid node
    int idx_valid = -1;
    for (int i = 0; i < y_values.size(); i++) {
        if (!std::isfinite(y_values[i])) continue;
        idx_valid = i;
        break;
    }

    // nothing found
    if (idx_valid == -1) return;

    // do extrapolation
    const auto steps = idx_valid;
    const double dy = (extrapolation_target - y_values[idx_valid]) / steps;
    for (int i = idx_valid-1; i >= 0; i--) {
        y_values[i] = y_values[i+1] + dy;
    }
}
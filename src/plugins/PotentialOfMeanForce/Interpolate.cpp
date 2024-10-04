#include "Interpolate.h"

Interpolate::Interpolate(double def_val, bool checks):default_value{def_val}{

}

// void Interpolate::FirstValid(){
    
//     int i=0;
//     while(!std::isfinite(y_values[i])){
//         i++;
//     }
//     if(i>=y_values.size()){
//         Log::global_log->error()<<"[Interpolate]Out of bounds"<<std::endl;
//     }
//     first_valid = i;
// }

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
        if(x_values[i]>1.0){
            y_values[i] -= v[i];
        }
    }

}

double Interpolate::CentralFiniteDifference(double r){
    // Return zero since we are at the limit and function is flat
    if(r > _simulation.getcutoffRadius()){
        Log::global_log->info()<<"Why are we checking beyond r_c="<<x_values[x_values.size()-1]<<"\t r="<<r<<std::endl;
        return 0.0;
    }

    int low, up;
    up = GetUpperLimit(r);

    
    // if(up <= first_valid && check_for_isfinite){
    //     return derivative_limit;
    // }


    low = up -1;

    double gb,ga,ra,rb;
    ga = y_values[low]; ra = x_values[low];
    gb = y_values[up]; rb = x_values[up];
    double ratio =(gb-ga)/(2.0*(rb-ra));

    // if(!std::isfinite(ratio)){
    //     return derivative_limit;
    // }

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

    for(int i= y_values.size()-1;i>=0;--i){

        if(std::isfinite(y_values[i])){
            continue;
        }

        double y_k = y_values[i+1];
        double y_k_1 = y_values[i+2];
        double x_k = x_values[i+1];
        double x_k_1 = x_values[i+2];
        double x_star = x_values[i];
        y_values[i] = y_k_1 + 1.0*(x_star-x_k_1)/(x_k-x_k_1)*(y_k-y_k_1);
        
    }

}
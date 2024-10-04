#include "Filter.h"

std::vector<double> Filter::MovingAverage(std::vector<double>& data){
    int data_size = data.size();
    std::vector<double> filtered_data;
    filtered_data.resize(data_size);

    for(int i=0;i<data_size;++i){
        double current_value=0.0;
        int count=0;
        for(int j=-2;j<=2;++j){
            if((i+j)<data_size && (i+j)>=0){
                current_value += data[i+j];
                count++;
            }
        }
        filtered_data[i] = current_value/(double)count;
    }

    return filtered_data;

}
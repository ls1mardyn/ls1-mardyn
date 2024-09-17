/*
 * Math.h
 *
 *  Created on: 11.11.2016
 *      Author: mheinen
 */

#ifndef MATH_H_
#define MATH_H_

#include <vector>
#include <numeric>

template<typename T>
struct square
{
    T operator()(const T& Left, const T& Right) const
    {
        return (Left + Right*Right);
    }
};

template<typename T1, typename T2>
inline void LinearRegression(const std::vector<T1>& x, const std::vector<T1>& y, T2& beta1, T2& beta2)
{
	// Reference: https://de.wikipedia.org/wiki/Lineare_Regression
	std::vector<T1> xy;
	typename std::vector<T1>::const_iterator it;
	for(unsigned int i=0; i<x.size(); i++)
		xy.push_back(x.at(i)*y.at(i) );

	T2 n = (T2)(x.size() );
	T2 invn = 1./n;
	T2 sumx = std::accumulate(x.begin(), x.end(), 0.0);
	T2 sumy = std::accumulate(y.begin(), y.end(), 0.0);
	T2 sumxy = std::accumulate(xy.begin(), xy.end(), 0.0);
	T2 sumx2 = std::accumulate(x.begin(), x.end(), 0.0, square<T2>() );
	T2 avgx = sumx*invn;
	T2 avgy = sumy*invn;

	T2 numerator = (n*sumxy - sumx*sumy);
	T2 dominator = (n*sumx2 - sumx*sumx);

	if(dominator > 0.)
	{
		beta2 = numerator / dominator;
		beta1 = avgy - beta2*avgx;
	}
	else
		beta2 = beta1 = 0.;
}

/**
 * @brief Used for equality checks when comparing floats.
 * 
 * Taken from AutoPas - src/autopas/utils/Math.h
 * 
 * @param a 
 * @param b 
 * @param maxRelativeDifference optional, default value 1e-9
 * @return true 
 * @return false 
 */
bool inline isNearRel(double a, double b, double maxRelativeDifference = 1e-9) 
{
	const auto greaterNumber = std::max(std::abs(a), std::abs(b));
	const auto absoluteDifference = maxRelativeDifference * greaterNumber;
	const auto diff = std::abs(a - b);
	return diff <= absoluteDifference;
}
#endif /* MATH_H_ */

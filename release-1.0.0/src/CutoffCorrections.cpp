/*
 * CutoffCorrections.h
 */

#include <cmath>

#include "CutoffCorrections.h"


double TICCu(int n,double rc,double sigma2)
{
	return -pow(rc,2*n+3) / (pow(sigma2,n)*(2*n+3));
}

double TICSu(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+3) - pow(rc-tau,2*n+3) ) * rc / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3) ) 
		+  ( pow(rc+tau,2*n+4) - pow(rc-tau,2*n+4) ) / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}

double TISSu(int n,double rc,double sigma2,double tau1,double tau2)
{
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -( pow(rc+tauPlus,2*n+4) - pow(rc+tauMinus,2*n+4) - pow(rc-tauMinus,2*n+4) + pow(rc-tauPlus,2*n+4) ) * rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) ) +  ( pow(rc+tauPlus,2*n+5) - pow(rc+tauMinus,2*n+5) - pow(rc-tauMinus,2*n+5) + pow(rc-tauPlus,2*n+5) ) / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}

double TICCv(int n,double rc,double sigma2)
{
	return 2*n * TICCu(n,rc,sigma2);
}

double TICSv(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+2) - pow(rc-tau,2*n+2) ) * rc*rc / ( 4*pow(sigma2,n)*tau*(n+1) ) - 3*TICSu(n,rc,sigma2,tau);
}

double TISSv(int n,double rc,double sigma2,double tau1,double tau2){
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -(   pow(rc+tauPlus,2*n+3) - pow(rc+tauMinus,2*n+3) - pow(rc-tauMinus,2*n+3) + pow(rc-tauPlus,2*n+3) ) * rc*rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) ) - 3*TISSu(n,rc,sigma2,tau1,tau2);
}

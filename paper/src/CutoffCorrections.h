#ifndef __CUTOFF_CORRECTIONS_H__
#define __CUTOFF_CORRECTIONS_H__


/* TODO: Comments on all the functions */
double TICCu(int n,double rc,double sigma2);
double TICSu(int n,double rc,double sigma2,double tau);
double TISSu(int n,double rc,double sigma2,double tau1,double tau2);
double TICCv(int n,double rc,double sigma2);
double TICSv(int n,double rc,double sigma2,double tau);
double TISSv(int n,double rc,double sigma2,double tau1,double tau2);

#endif /* __CUTOFF_CORRECTIONS_H__ */

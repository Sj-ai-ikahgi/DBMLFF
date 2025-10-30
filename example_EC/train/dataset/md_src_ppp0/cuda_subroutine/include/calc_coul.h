#ifndef CALC_COUL_H
#define CALC_COUL_H

#ifdef __cplusplus
extern "C" {
#endif

void fcc_calc_coul(double *rho, double *vcoul, int n1, int n2, int n3, double *AL, double pi);

#ifdef __cplusplus
}
#endif

#endif
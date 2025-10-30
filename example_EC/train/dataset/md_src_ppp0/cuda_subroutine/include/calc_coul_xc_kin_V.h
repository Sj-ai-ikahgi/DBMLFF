#ifndef CALC_COUL_XC_KIN_V_H
#define CALC_COUL_XC_KIN_V_H

#ifdef __cplusplus
extern "C" {
#endif

void fcc_calc_coul_xc_kin_V(double *rho, double *rho_z, double *vxc, double *vxc2, double *vcoul, int n1, int n2, int n3, double *AL,
                            double *E_coul, double *E_xc, double *E_kin1, double *E_kin2, double fact_kin2, double pi);

#ifdef __cplusplus
}
#endif

#endif
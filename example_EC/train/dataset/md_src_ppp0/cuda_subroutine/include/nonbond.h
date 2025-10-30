#ifndef NONBOND_H
#define NONBOND_H

#ifdef __cplusplus
extern "C" {
#endif

void fcc_nonbond_loop_82to96(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1);
void fcc_nonbond_loop_134to149(int *icent, double *w_cent, int *nat_cent, double *xatom_m, int itype_mol, int imol, double *xatom_cent, 
                        int ncent_itype_mol, int natom_mm, int nmolm);
void fcc_nonbond_loop_149to238(double *xatom_cent, double *z_cent, int *itype_cent, double *funcr2, double *AL, int ncent_itype_mol,
                        int nr, double Rm2, double *dxyz_box, double *rho, double *rho_z, int id3, int id2, int id1,
                        int n3, int n2, int n1, double Rbox2, int itype_mol, double pi, int *ion_type_cent, int *imax_ion,
                        double *r_ion, double *rho_ion, double *box, double *box2, double *box3, double vol, double *Q_type, double *z_ion,
                        int imax_nr, int imax_ntype_cent);
void fcc_nonbond_loop_400(int *itype_cent, double *xatom_cent, double * z_cent, int *icorner, double *dbox_c, double *dbox3_c, double *dxyz_box, double *AL, 
                    double *rho, double *rho_z, double *funcr2, int nr, double Rm2, int id3, int id2, int id1, int nm3, int nm2, int nm1, int n3, int n2, int n1, 
                    double Rbox2, int itype_mol, int ncent_itype_mol, double pi, int *ion_type_cent, int *imax_ion, double *r_ion, double *rho_ion, double *box, 
                    double *box2, double *box3, double vol, double *Q_type, double *z_ion, double *dbox, double *dbox2, double *dbox3, int imax_nr, int imax_ntype_cent);
void fcc_nonbond_loop_500(int *itype_cent, double *xatom_cent, double *force_cent, int *icorner, double *vxc, double *vxcm, double *vcoul, double *vcoulm, double *dbox_c,
                    double *dbox3_c, double vol_n, int id3, int id2, int id1, int n3, int n2, int n1, int nm3, int nm2, int nm1, int itype_mol, int ncent_itype_mol);
void fcc_nonbond_loop_510(int natom_m_itype_mol, int itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int n3, int n2, int n1, int nm3, int nm2, int nm1, int *icorner, double *AL, 
                    double vol_n, int id3, int id2, int id1, double *dxyz_box, double Rbox2, double *pxyz, double *dpxyz, double *vxc2, double *vxc2_m, double *vcoul, double *vcoul_m, int imol, int natom_mm,
                    int nmolm);

#ifdef __cplusplus
}
#endif

#endif
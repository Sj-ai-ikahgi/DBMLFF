#ifndef CALC_POLAR_FORCE_H
#define CALC_POLAR_FORCE_H

#ifdef __cplusplus
extern "C"{
#endif

void fcc_polar_loop_113to129(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1);
void fcc_polar_loop_500(int natom_m_itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int imol, int itype_mol, int n3, int n2, int n1, int id3, 
                        int id2, int id1, int natom_mm, int nmolm, int ntype_mm, double *AL, double *dxyz_box, double *CC_pol, double *rhop_tot, double Rbox2);
void fcc_polar_loop_312to339(int ncent_itype_mol, int *icent, int imol, int itype_mol, double *w_cent, double *xc_cent, double *xatom_cent, int *icorner, 
                            int n3, int n2, int n1, int nm3, int nm2, int nm1, double *xatom_m, int natom_mm, int nmolm, int *nat_cent);
void fcc_polar_loop_600(int ncent_itype_mol, int *itype_cent, int itype_mol, double *xatom_cent, int *icorner, double *AL, int *ion_type_cent, int *imax_ion, double *r_ion, int id3, int id2, 
                        int id1, double *box, double *box2, double *box3, double *dxyz_box, double *funcr2, double *fact_store, double *fact2_store, int n3, int n2, int n1, int nr, double Rm2, 
                        double Rbox2, double *rho_ion, double vol, int nm3, int nm2, int nm1, double *rho_m, double *Q_type, double *z_ion, int imax_nr, int imax_ntype_cent);
void fcc_polar_loop700(int natom_m_itype_mol, int *ion_type_atomp, int imol, int itype_mol, double *atom_charge_param, double *xatom_m, int *icorner, int n3, int n2, int n1,
                        int id3, int id2, int id1, double *CC_pol, double *rhop_m, double *dxyz_box, double Rbox2, int natom_mm, int nmolm, double *AL, int nm3, int nm2, int nm1);
void fcc_polar_loop800(int ncent_itype_mol, int * itype_cent, double *xatom_cent, int *icorner, double *AL, int n3, int n2, int n1, int *ion_type_cent, int *imax_ion, double *r_ion, 
                        double *dxyz_box, double *funcr2, double *box, double *box2, double *box3, double *dbox, double *dbox2, double *dbox3, double *fact_store, double *fact2_store, 
                        double Rbox2, int nr, int nm3, int nm2, int nm1, double *vrhop_tot_coul, double *vrhop_m_coul, double *dvxc_tot, double *dvxc_m, double vol_n, int id3, int id2,
                        int id1, double *force_cent, double Rm2, int itype_mol, double *rho_ion, int imax_nr, int imax_ntype_cent);
                        
#ifdef __cplusplus
}
#endif

#endif
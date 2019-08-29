#ifndef KEYFUNG_H
#define KEYFUNG_H 1

#define MXLGVN 20000 //maximum allowed langevin dipoles, should use dynamic arrays in gen_gridx. 
#define MXCENTER 50 //needed by ran_shift and elgvn_ave 
#define MXATM 984 //max amt of atoms allowed, should be dynamic 
#define MXPAIR 100000000 //max amt of pairs allowed 

/*atom types: 
1/H0, 2-He, 3-Li, 4-Be, 5-B, 6-C0, 7-C1, 8-C2, 9-N0, 10-N1,
11-N2, 12-X, 13-O0, 14-O1, 15-O2, 16-F, 17-Ne, 18-Na, 19-Mg, 20-Al,
21-Si, 22-P, 23-X, 24-S, 25-Cl, 26-Ar, 27-K, 28-Ca, 29-Sc, 30-Ti, 
31-V, 32-Cr, 33-Mn, 34-Fe, 35-Co, 36-Ni, 37-Cu, 38-Zn, 39-Ga,40-Ge,
41-As, 42-Se,43-Br, 44-Kr, 45-Rb, 46-Sr, 47-Y, 48-Zr, 49-Nb, 50-Mo, 
51-Tc, 52-Ru, 53-Rh, 54-Pd, 55-Ag, 56-Cd, 57-In, 58-Sn, 59-Sb, 60-Te, 
61-I, 62-Xe, 63-Cs, 64-Ba, 65-La, 66-Hf, 67-Ta, 68-W, 69-Re, 70-Os, 
71-Ir, 72-Pt, 73-Au, 74-Hg, 75-Tl, 76-Pb, 77-Bi, 78-Po, 79-At, 80-Rn, 
81-Fr, 82-Ra.*/
extern const double kB;
extern const double Na;
extern const double h;
extern const double pi;
extern const double e;
extern const double amu; 
extern const int iac_conv[];
extern const int mass[];
extern const double rp[];

typedef struct mol_t{
  //core topology
	char **atom; //chemical symbol
	double *zan; //atomic number 
	double *q;   //atomic charge 
	double *xw[3]; //cartesian coordinates of nucleus 
	//van der waals radii topology 
	double *rpi;
	//optional charge topology
	double *q_gas;
	double *q_mp2;
	//chemsol atom type topology
	int *iacw;
	double vdwc6[82];
  double srp;
	double ebw, elgvn, elgwa, erelax, evdw, evqdq;
	double rg_reg1, rgim; 
	double drg, rg, rg_inner, drg_inner;
	double dxp0[3];
	double rdcutl, out_cut;
	int ndxp, itl, itp; 
	double ephil1, ephil2, ephob, fsurfa[MXCENTER], evdwl[MXCENTER];
	double vdwsl, rzcut, phobsl, tdsl[MXCENTER], etds, amas, tds0;
	double clgvn, slgvn;
	double pcenter[3];
	char dumm1[8]; 
	char molname[13]; 
	char ssname[4]; 
	char fname[256]; 
	char relax[3]; 
	int n_reg1, n_inner, ngeom;
  int iter_id; 
}mol_t;

typedef struct lgvn_t{
	double esum, tdsw_a, elgvni, tds, vbornx_result, ecor;
  double temp_elgvn[MXCENTER];
  double tdsl[MXCENTER];
  double efa[3][MXLGVN];
  double efal[3][MXLGVN];
  double temp_center[3][MXLGVN];
  double xmua[3][MXLGVN];
  double **xl, *rz1, *rz_vdw, *atomfs;
  int *isd, *iz, jp, jp2;
  int n_inner, ndipole;
  int ip[0][MXLGVN], ip2[0][MXLGVN], ip3[0][MXLGVN];
  int nvol[MXCENTER];
  double **da; //double[3][MXLGVN]
  double *oshift; //double[3*MXCENTER]
  double *center_new; //double[3]
  double *vatom_result; //double[2]
  double *elgvn_ave_result; //double[4]
}lgvn_t;

int readopt(mol_t *mol, char *vdw_name);
double ran2(int seed); 
double *generate_oshift(int ndxp);
double *ran_shift(mol_t *mol, double *oshift, int i);
int lgvnx(mol_t *mol, lgvn_t *lgvn_data, int ientro);
int gen_gridx(mol_t *mol, lgvn_t *lgvn_data, int ientro);

#endif
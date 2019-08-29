#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double kB = 1.3806488e-23; // boltzmanns constant, in m^2 kg s^-2 K^-1
const double Na = 6.023e23; // avogadros constant 
const double h = 6.626e-34; // plancks constant, in joule seconds
const double pi = 3.14159; 
const double e = 2.71828; // eulers number 
const double amu = 1.66e-27; //the dalton, in kg

const int iac_conv[89] = {1,2, \
 3,4,5,6,9,13,16,17, \
 18,19,20,21,22,24,25,26, \
 27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44, \
 45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62, \
 63,64,65, \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
 66,67,68,69,70,71,72,73,74,75,76,77,78,79,80, \
 81,82};
const int mass[88] = {1,4, \
 7,9,11,12,14,16,19,20, \
 23,24,27,28,31,32,35,40, \
 39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84, \
 85,88,89,91,93,96,97,101,103,106,108,112,115,119,122,128,127,131, \
 132,137,139, \
 140,141,144,145,150,152,157,159,163,165,167,169,173,175, \
 178,181,184,186,190,192,195,197,201,204,207,208,209,210,222, \
 223,226};
const double rp[82] = {
  //  H   He  Li   Be   B    C    C1  C2   N    N1
       2.3,2.3,2.15,2.00,2.40,2.65,3.0,3.25,2.65,2.85, \
    //  N2   X   O   O1  O2    F    Ne  Na   Mg  Al
       3.2,3.0,2.32,2.65,2.8,2.46,2.5,2.58,1.82,1.70, \
    //  Si  P    X  S    Cl  Ar   K    Ca   Sc  Ti
       3.1,3.2,3.0,3.2,3.16,2.8,3.06,2.38,1.5,2.00, \
    //  V    Cr   Mn   Fe   Co   Ni  Cu1+  Zn  Ga   Ge
       1.91,1.89,1.93,1.84,1.57,1.50,1.88,1.55,2.00,2.50, \
    //  As   Se    Br   Kr   Rb  Sr    Y    Zr  Nb    Mo
       3.00,3.00,3.44,3.00,3.25,2.70,1.75,2.00,2.00,2.00, \
    //  Tc    Ru   Rh  Pd   Ag    Cd   In   Sn   Sb  Te
       2.00,2.00,2.00,2.00,2.25,1.98,2.45,2.44,3.00,3.70, \
    //  I    Xe   Cs   Ba   La   Hf   Ta   W    Re   Ir 
       3.80,3.55,3.58,2.92,2.30,3.00,3.00,3.00,3.00,3.00, \
    //  Pt   Au   Hg   Tl   Pb   Bi   Po   At   Rn   Fr
       3.00,1.77,1.94,2.00,2.55,3.00,3.00,3.00,3.00,3.00, \
    //  Ra   A
       3.50,3.00
};

int main(int argc, char **argv){
  char input_file_name[256]; 
  char vdw_name[256];
  mol_t mol;
  for(int i = 0; i<3; i+=1){
    mol.dxp0[i] = 0;
  }
	mol.rg = 26.0;
  mol.drg = 3.0;
  mol.drg_inner = 1.0;
  mol.rdcutl = 6.1;
  mol.out_cut = 16.0;
  mol.rgim = 2.0;
  mol.itl = 399;
  mol.itp = 5;
  mol.rzcut = 1.1;
  strcpy(vdw_name, argv[1]);
  int readopt_fail = readopt(&mol, vdw_name);
  FILE *input_ptr;
  if(readopt_fail){
    printf("Could not open parameter file!\n");
    exit(1);
  }
  else{
  	strcpy(input_file_name, argv[2]);
  }
  if((input_ptr = fopen(input_file_name, "r"))){
  	char temp[256];
  	fgets(mol.molname, 100, input_ptr);
  	fgets(temp, 100, input_ptr);
  	//printf("%s %s\n", molname, temp); debug 
  	mol.n_reg1 = atoi(strtok(temp, " "));
  	mol.ngeom = atoi(strtok(NULL, " "));
  	//printf("%s %d %d\n", molname, n_reg1, ngeom); debug 
  }
  else{
  	printf("error reading molecule name from file!\n");
  	exit(1);
  }
  mol.atom = malloc(mol.n_reg1*4);		  // *char[3]
  mol.zan = malloc(sizeof(double) * mol.n_reg1);
  mol.q = malloc(sizeof(double) * mol.n_reg1);
  mol.xw[0] = malloc(sizeof(double) * mol.n_reg1); // *double[3]
  mol.xw[1] = malloc(sizeof(double) * mol.n_reg1); // *double[3]
  mol.xw[2] = malloc(sizeof(double) * mol.n_reg1); // *double[3]
  mol.iacw = malloc(sizeof(int) * mol.n_reg1);
  mol.q_gas = malloc(sizeof(double) * mol.n_reg1);
  size_t len; 
  char *temp =NULL;
  mol.amas = 0; 
  printf("name %-2s # %-4s charge %-6s x %-8s y %-8s z %-5s i %-1s iacw\n", "", "", "", "", "", "", "");
  for(int i =0; i<mol.n_reg1; i+=1){
  	getline(&temp, &len, input_ptr);
  	//printf("%s", temp);
  	if(strtok(temp,"\n")==NULL || \
      strtok(temp, "\r")==NULL || \
      strtok(temp, "\r\n")==NULL){
  		i-=1;
      continue;
    }
  	mol.atom[i] = strtok(temp, " ");
  	mol.zan[i] = strtod(strtok(NULL, " "), NULL);
  	mol.q[i] = strtod(strtok(NULL, " "), NULL);
  	mol.xw[0][i] = strtod(strtok(NULL, " "), NULL);
  	mol.xw[1][i] = strtod(strtok(NULL, " "), NULL);
  	mol.xw[2][i] = strtod(strtok(NULL, " "), NULL);
  	printf("%-5s %4.1lf %10.4lf %10.4lf %10.4lf %10.4lf ", mol.atom[i], mol.zan[i], mol.q[i], mol.xw[0][i], mol.xw[1][i], mol.xw[2][i]);
  	mol.iacw[i] = (int)iac_conv[(int)mol.zan[i]];
    printf("%5d ", i);
    printf("%5d\n", mol.iacw[i]);
  	mol.amas += (double)mass[(int)mol.zan[i] - 1];
  }
  printf("===\nMass: %.1fu\nn_reg1: %d\nngeom: %d\n", mol.amas, mol.n_reg1, mol.ngeom); 

  //If charges from a PCM calculation (q_pcm) are available, they will
  //be used in the evaluation of dGlgvn. Charges can be inputed in two ways:

  //1/ q_pcm are given as the first set of charges in the
  //input file and no second set of charges is given. This scheme
  //is useful when calculating charges using the Gaussian98 code
  //to avoid extra calculation of the gas-phase charges.

  //2/ q_gas charges are given as the first set of charges and q_pcm
  //are given as the second set of charges. This is the input format 
  //used in Chemsol 1.0 - 2.0. The advantage of this scheme is that
  //the additional information about the gas-phase charge distribution 
  //can be used to evaluate EXPLICITLY the contribution of the solute polarization
  //by the solvent to the total solvation free energy (dGsolv). This 
  //contribution is called here dGrelax. Because dGrelax is implicitly
  //included in dGsolv if pcm charges are available, dGsolv results
  //are not affected by the presence/absence of the gas_phase charges 
  //in the input file.

  /* 
  			IMPLEMENT GAS LOGIC HERE 
  			(main.f90 line 155)
  */

  //check consistency of charges 
  for(int i = 0; i<mol.n_reg1; i+=1){
  	if(fabs(mol.q[i])>9){
  		printf("inconsistent charges!\n");
  		exit(1);
  	}
  }
  //the second character of the atom name is used to make a wider 
  //selection of chemsol atom types that differ in VdW radii.
  for(int i = 0; i<mol.n_reg1; i+=1){
  	if(mol.iacw[i]!=1){
  		if(mol.atom[i][1]=='1')
        mol.iacw[i]+=1;
      if(mol.atom[i][1]=='2')
        mol.iacw[i]+=2;
      if(mol.atom[i][1]=='3')
        mol.iacw[i]+=3;
  	}
    else if(mol.iacw[i]==0){
      printf("Unknown atom type on input:\n i: %d\n atom: %s\n #: %.1lf\n", i, mol.atom[i], mol.zan[i]);
      exit(1);
    }
  }
  exit(1);
  free(mol.atom); 
  free(mol.zan); 
  free(mol.q);
  free(mol.xw); 
  free(mol.iacw);
  free(mol.q_gas); 
  fclose(input_ptr);
  return 0; 

}
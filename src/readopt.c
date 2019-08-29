#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int readopt(mol_t *mol, char* vdw_name){
  FILE *vdw_ptr;
  double r2, r2min, dmax; 
  int nrp, jmin;
  if(!(vdw_ptr = fopen(vdw_name, "r")))
    return 1; 
	for(int i = 0; i<82; i+=1){
		mol->vdwc6[i] = 0; 
	}
  mol->iter_id = 1;
  mol->ndxp = 20;
  mol->dxp0[0] = 0.2;
  mol->clgvn = 0.9; 
  mol->slgvn = 0.48;
  mol->tds0 = 1.4; 
  mol->srp = 0.88;
  mol->vdwc6[0] = 1;
  mol->vdwc6[5] = 1;
  mol->vdwc6[6] = 1;
  mol->vdwc6[8] = 0.5;
  mol->vdwc6[12] = 0.5;
  mol->vdwc6[21] = 1;
  mol->vdwc6[23] = 1;
  mol->vdwsl = 0.120; 
  mol->phobsl = 0.138/10;
  mol->ephil1 = 0.002; 
  mol->ephil2 = 0.038;
  mol->rzcut = 1.5; 
  fclose(vdw_ptr); 
  //define the remaining london coefficients 
  for(int i = 0; i<mol->n_reg1; i+=1){
    if(mol->vdwc6[i] == 0){
      if(mol->iacw[i] <= 17){
        mol->vdwc6[i] = mol->vdwc6[5];
      }
      else{
        mol->vdwc6[i] = mol->vdwc6[21];
      }
    }
  }
  //calculate rp(H) as a linear function of the rp of the atom to which
  //the hydrogen is covalently bonded. Note that the coefficient of this
  //function (srp) differs for the 1st and 2nd row atoms. 
  //in addition, for inorganic oxygen (iacw=15), a separate rp(H) is used. 
  for(int i = 0; i<mol->n_reg1; i+=1){
    if(mol->iacw[i] != 1 && mol->iacw[i] != 2){
      mol->rpi[i] = rp[mol->iacw[i]];
    }
    else{
      r2min = 100;
      jmin = 21; 
      for(int j = 0; j < mol->n_reg1; j+=1){
        if(j==i)
          continue; 
        r2 = pow((mol->xw[0][i] - mol->xw[0][j]),2) + pow((mol->xw[1][i]-mol->xw[1][j]),2) \
        + pow((mol->xw[2][i]-mol->xw[2][j]),2);
        if(r2 < r2min){
          r2min = r2; 
          jmin = j; 
        }
      }
      if(mol->iacw[jmin] == 15){ //if the atom is attached to O2
        mol->iacw[i] = 2; //He and H have the same radii 
                        //iacw only ever used again to calculate vdw enthalpy 
                        //in gen_gridx; He and H are treated the same in Chemsol. 
        mol->rpi[i] = rp[mol->iacw[i]]; //use He radius (look at rp in main)
        mol->vdwc6[mol->iacw[i]] = mol->vdwc6[mol->iacw[i]-1]; //need to decouple vdwc6 from rp mapping
                                                //what happens if mol->vdwc6[mol->iacw[i]]
                                                //is used for another atom?
      }
      else{
        if(mol->iacw[jmin] < 18)
          mol->rpi[i] = mol->srp * rp[mol->iacw[jmin]]; //if the atom is attached to a second period atom
        else
          mol->rpi[i] = (mol->srp - 0.1) * rp[mol->iacw[jmin]]; //if the atom is attached to a third period or greater atom
      }
    }
  }

  //Finally, allow for the change of rp parameter of any atom without changing
  //its atom type. A new rp is specified in the end of the input file. 
  /* 
          IMPLEMENT HERE (chemsol.f90 -- line 1679)
                                                        */

  for(int i = 0; i<3; i+=1){
    mol->pcenter[i] = 0;
  }
  for(int i = 0; i < mol->n_reg1; i+=1){
    mol->pcenter[0] += mol->xw[0][i]/mol->n_reg1;
    mol->pcenter[1] += mol->xw[1][i]/mol->n_reg1;
    mol->pcenter[2] += mol->xw[2][i]/mol->n_reg1; 
  }

  //find maximal radius of the solute wrt grid center
  dmax = 0; 
  for(int i = 0; i < mol->n_reg1; i+=1){
    mol->rg_reg1 = mol->pcenter[0] - pow(mol->xw[0][i],2) \
    + mol->pcenter[1] - pow(mol->xw[1][i],2) \
    + mol->pcenter[2] - pow(mol->xw[2][i],2);
    if(mol->rg_reg1 > dmax)
      dmax = mol->rg_reg1;
  }
  mol->rg_reg1 = sqrt(dmax) + 2.40;
  //new grid radii are measured wrt solute surface. This ensures that 
  //enough space is attributed to the grid for large molecules.
  //For actual choice of grid extension, other criteria are applied: 
  //The distance from the VdW surface, (inner, 1A grid) and magnitude of the field
  //at the grid point (outer, 3A grid) - see gen_gridx function. 
  mol->rg += mol->rg_reg1;
  mol->rg_inner = mol->rg_reg1 + mol->rgim + 2.0;

  if(mol->drg < 3){
    printf("Error: DRG should be no less than 3.0A\n");
    exit(1);
  }
  else if(mol->ndxp > MXCENTER){
    printf("Error: Maximum number of grids is %d\n", MXCENTER);
    exit(1);
  }
  return 0; 
}
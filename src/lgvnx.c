#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int lgvnx(mol_t *mol, lgvn_t *lgvn_data, int ientro){
  double vdwsur[MXATM];
  double sres, tds, efn_max, gri_sp, elgvn_i, efn;
  double *vlgvn_result; // double[3]
  double fma, elgvna, ddd, dddx, dddy, dddz, epot, rx, ry, rz;
  double rqd, fs; 
  int idum, i, k;

  mol->elgvn = 0.0; //total langevin energy
  // ephil1 and ephil2 are defined in readopt
  sres = 0.0;
  // sres is surface for large fields. 
  mol->fsurfa[ientro] = 0.0; //ientro = 0 at first vlgvn call
  mol->evdwl[ientro] = 0.0; 
  lgvn_data->ndipole = 0;
  idum = 0; 
  tds = 0.0;
  gen_gridx(mol, lgvn_data, ientro);
  // calculate magnitudes of langevin dipoles and noniterative solvation 
  // energy from the electric field scaled by the distance-dependent dielectric
  lgvn_data->da = ef_ld(mol, lgvn_data, 1); // ?? 
  efn_max = -10.0;
  for(i = 0; i < mol->n_reg1; i+=1){
    vdwsur[i] = 0.0;
  }
  for(i = 0; i < lgvn_data->ndipole; i+=1){
    gri_sp = mol->drg_inner; 
    if(i > mol->n_inner)
      gri_sp = mol->drg; 
    if(i == (mol->n_inner + 1))
      elgvn_i = mol->elgvn; // save elgvn_i for inner grid langevin dipole energies report

    efn = sqrt( pow(lgvn_data->da[0][i],2) + pow(lgvn_data->da[1][i],2) + pow(lgvn_data->da[2][i],2) );
    if(efn > efn_max)
      efn_max = efn; // set a new efn_max if a greater one is found
    vlgvn_result = vlgvn_f(efn, gri_sp, mol->slgvn);
    fma = vlgvn_result[0];
    tds = vlgvn_result[1]; 
    mol->elgvn += vlgvn_result[2]; 
    lgvn_data->xmua[0][i] = fma * lgvn_data->da[0][i] / efn; 
    lgvn_data->xmua[1][i] = fma * lgvn_data->da[1][i] / efn;
    lgvn_data->xmua[2][i] = fma * lgvn_data->da[2][i] / efn;

    // calculate hydrophobic surface
    if(i < mol->n_inner){
      if(lgvn_data->rz_vdw[i] >= mol->rzcut){
        // calculate electrostatic potential at the grid point
        epot = 0.0;
        for(k = 0; k < mol->n_reg1; k+=1){
          rx = lgvn_data->xl[0][i] - mol->xw[0][k];
          ry = lgvn_data->xl[1][i] - mol->xw[1][k];
          rz = lgvn_data->xl[2][i] - mol->xw[2][k];
          rqd = sqrt( pow(rx,2) + pow(ry,2) + pow(rz,2));
          epot += mol->q[k] / rqd;
        }
        if(epot < 0.0)
          epot = -epot; // positive surfaces and the surface for the vdw term 
        if(epot <= mol->ephil1)
          fs = 1.0;
        if(epot > mol->ephil1 && epot <= mol->ephil2)
          fs = 1.0 - ((epot - mol->ephil1) / (mol->ephil2 - mol->ephil1)) * (1.0 - sres);
        else if(epot > mol->ephil2)
          fs = sres; 
        lgvn_data->atomfs[lgvn_data->iz[i]] += fs; 
        mol->fsurfa[ientro] += fs; 
        vdwsur[lgvn_data->iz[i]] += 1.0;
      }
    }
  }
  // use atom polarizabilities (vdwc6) to calculate 
  // VdW part of the solution enthalpy. 
  mol->evdwl[ientro] = 0.0;
  for(k = 0; k < mol->n_reg1; k+=1){
    mol->evdwl[ientro] += mol->vdwc6[mol->iacw[k]] * vdwsur[k];
  }
  mol->evdwl[ientro] = mol->vdwsl * mol->evdwl[ientro];
  mol->elgvn *= mol->clgvn;
  elgvn_i *= mol->clgvn;

  // calculation of the initial configuration of Langevin dipoles 
  // (0th step of the iterative calculation of dipole-dipole interactions)
  lgvn_data->da = ef_ld(mol, lgvn_data, 0);
  elgvna = 0.0; 
  for(i = 0; i < lgvn_data->ndipole; i+=1){
    gri_sp = mol->drg_inner; 
    if(i > mol->n_inner)
      gri_sp = mol->drg;
    efn = sqrt( pow(lgvn_data->da[0][i],2) + pow(lgvn_data->da[1][i],2) + pow(lgvn_data->da[2][i],2) );
    vlgvn_result = vlgvn_f(efn, gri_sp, mol->slgvn);
    fma = vlgvn_result[0];
    tds = vlgvn_result[1];
    mol->elgvn += vlgvn_result[2];

    // dipole moments of point (langevin) dipoles are oriented along
    // the field from solute for inner grid and outer surface dipoles. 
    // They are oriented randomly for odd nonsurface outer grid dipoles. 
    // (Note that this is implemented for the 0 iteration step only.)
    // (For iterative lgvn, lgvn dipoles are calculated in subroutine 
    // mu_mu_l. Only those lgvn dipoles that lie along outer surfaces are
    // constrained to be proportional to the solute field.)

    if(i > mol->n_inner || lgvn_data->isd[i] == 1 || (i % 3) == 0){
      ddd = 3.0 / (sqrt( lgvn_data->rz1[i] + 2.0 ));
      lgvn_data->xmua[0][i] = ddd * fma * lgvn_data->da[0][i] / efn;
      lgvn_data->xmua[1][i] = ddd * fma * lgvn_data->da[1][i] / efn;
      lgvn_data->xmua[2][i] = ddd * fma * lgvn_data->da[2][i] / efn;
      continue;
    }
    else{
      dddx = ran2(idum); 
      dddy = ran2(idum); 
      dddz = ran2(idum); 
      ddd = sqrt( pow(dddx,2) + pow(dddy,2) + pow(dddz,2));
      lgvn_data->xmua[0][i] = dddx * fma / ddd;
      lgvn_data->xmua[1][i] = dddy * fma / ddd;
      lgvn_data->xmua[2][i] = dddz * fma / ddd;
    }
  }
  return 0; 
}

/*call lgvnx(center_new,elgvn,ndipole,i,iterld,fsurfa,lgvn_data->da,xmua,atomfs,iz,evdwl,
  xl,drg_inner,drg,n_inner,rz_vdw,rzcut,q,ephil1,ephil2,vdwc6,iacw,vdwsl,
  clgvn,slgvn,isd,rz1,n_reg1,nvol,rg,rg_inner,rgim,rpi,xw)

subr lgvnx(center1,elgvn,ndipole,ientro,iterld,fsurfa,lgvn_data->da,xmua,atomfs,iz,evdwl,
  xl,drg_inner,drg,n_inner,rz_vdw,rzcut,q,ephil1,ephil2,vdwc6,iacw,vdwsl,
  clgvn,slgvn,isd,rz1,n_reg1,nvol,rg,rg_inner,rgim,rpi,xw)*/
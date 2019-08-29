#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gen_gridx(mol_t *mol, lgvn_t *lgvn_data, int ientro, int iflag){
  // atom fit between "vdwsl" and "evdwl"
  int ns, i0, i1, i, limit_inner, limit_outer, mid_inner, no, index;
  int ii, jj, kk;
  int icube[133][133][133];
  int igrid, jgrid, kgrid;
  int ipro, jpro, kpro;
  double rg2, drg2, drgi2, rgi2, rgim2, drg_inner_i, drg_i; //rp2[MXATM]
  double ri, rj, rk, d1, d2, c6, d6, r3, r6; 
  double xloca, yloca, zloca;
  double xp[3];
  double d2_min, efmin1, efmin2, efx, efy, efz, ddd;
  double d3, qr, efnorm, subvdw, vdw_in, vdw_out, xdrg; 
  // generate grid point for langevin dipole 
  // center1 = center of the grid (input)
  // ndipole = number of grid points (output)

  // set up parameters 
  // rg is grid radius (outer), rgi is the same for the inner grid.
  // drg is outer grid spacing (3A), rgim is a trim p mxceparameter
  // for the inner grid. 

  lgvn_data->ndipole = 0; 
  ns = 0; 
  i0 = 0; 
  i1 = mol->n_reg1; 
  rg2 = mol->rg * mol->rg; 
  drg2 = mol->drg * mol->drg; 
  drgi2 = mol->drg_inner * mol->drg_inner; 
  rgi2 = mol->rg_inner * mol->rg_inner; 
  rgim2 = mol->rgim * mol->rgim; 
  drg_inner_i = 1.0 / mol->drg_inner;
  drg_i = 1.0 / mol->drg; 
  limit_inner = (int) (mol->rg_inner * 2 / mol->drg_inner);
  limit_outer = (int) (mol->rg * 2 / mol->drg);

  // make both limits into an odd number
  if(limit_inner % 2 == 0)
    limit_inner+=1;
  if(limit_outer % 2 == 0)
    limit_outer +=1; 
  mid_inner = (int) (limit_inner / 2.0 + 0.5);
  mid_outer = (int) (limit_outer / 2.0 + 0.5);

  //===============================================
  // build inner grid, if appropriate
  //===============================================

  if(mol->rg_inner > 0.0){ //build inner grid, rg_inner should ALWAYS > 0! 
    for(ii = 0; ii < limit_inner; ii+=1){
      for(jj = 0; jj < limit_inner; jj+=1){
        for(kk = 0; kk < limit_inner; kk +=1){
          icube[ii][jj][kk] = 0;
          ri = ii - mid_inner; 
          rj = jj - mid_inner; 
          rk = kk - mid_inner; 
          d2 = (ri*ri + rj*rj + rk*rk) * drgi2; //radius from center
          if(d2 <= rgi2)
            icube[ii][jj][kk] = 1; 
        }
      }
    }

    // make cavity 
    // put solute on the grid and find the nearest grid points 
    // out of these, grid points closer than rp will bne removed. 
    // the number of points removed (nvol) is proportional to the 
    // solute volume. 

    nvol[ientro] = 0; 
    for(i = i0; i < i1; i+=1){ // loop over solute atoms
      ipro = (int) (mid_inner + (mol->xw[0][i] - center1[0] * drg_inner_i));
      jpro = (int) (mid_inner + (mol->xw[1][i] - center1[1] * drg_inner_i));
      kpro = (int) (mid_inner + (mol->xw[2][i] - center1[2] * drg_inner_i));
      if(ipro >= 0 && ipro <= limit_inner \ 
        && jpro >= 0 && jpro <= limit_inner\
        && kpro >= 0 && kpro <= limit_inner){
        // check distances to 27 closest points around the 
        // nearest grid point including itself
        for(ii = 0; ii < 9; ii+=1){
          igrid = ipro + ii - 5;
          for(jj = 0; jj < 9; jj+=1){
            jgrid = jpro + jj - 5; 
            for(kk = 0; kk < 9; kk+=1){
              kgrid = kpro + kk - 5;
              if(icube[igrid][jgrid][kgrid] == 1){
                xloca = (igrid - mid_inner) * (drg_inner + center1[0]);
                yloca = (jgrid - mid_inner) * (drg_inner + center1[1]);
                zloca = (kgrid - mid_inner) * (drg_inner + center1[2]);
                ri = xloca - mol->xw[0][i];
                rj = yloca - mol->xw[1][i];
                rk = zloca - mol->xw[2][i];
                d2 = ri*ri + rj*rj + rk*rk // distance from grid point
                if(d2 < (rpi[i]*rpi[i])){
                  icube[igrid][jgrid][kgrid] = 0;
                  nvol[ientro] += 1;
                }
              } 
            }
          }
        }
      }
    }

    // reshape grid sphere (outer region) to a solute envelope/
    for(ii = 0; ii < limit_inner; ii+=1){
      for(jj = 0; jj < limit_inner; jj+=1){
        for(kk = 0; kk < limit_inner; kk+=1){
          if(icube[ii][jj][kk] != 0){
            xp[0] = center1[0] + (ii - mid_inner) * drg_inner;
            xp[1] = center1[1] + (jj - mid_inner) * drg_inner;
            xp[2] = center1[2] + (kk - mid_inner) * drg_inner;
            d2_min = 10000.0; 
            // get distance to a closest VdW boundary 
            for(i = i0; i < i1; i+=1){
              ri = xp[0] - mol->xw[0][i];
              rj = xp[1] - mol->xw[1][i];
              rk = xp[2] - mol->xw[2][i];
              d2 = ri*ri + rj*rj + rk*rk;
              d2 = sqrt(d2);
              d2 = d2 - mol->rpi[i];
              if(d2 > 0.0){
                printf("Error: Bad grid point placement in gen_gridx!\n");
                exit(1);
              }
            }
          }
        }
      }
    }

  }
}

/*call gen_gridx(center1,ndipole,ientro,0,nvol,isd,xl,n_inner,n_reg1,rg,
  drg,drg_inner,rg_inner,rgim,rpi,xw,iacw,vdwc6,vdwsl,evdwl,rz1,iz,rz_vdw,q)

subr gen_gridx (center1,ndipole,ientro,iflag,nvol,isd,xl,n_inner,n_reg1,rg,
  drg,drg_inner,rg_inner,rgim,rpi,xw,iacw,vdwc6,vdwsl,evdwl,rz1,iz,rz_vdw,q)*/
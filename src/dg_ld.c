#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int dg_ld(mol_t *mol, lgvn_t *lgvn_data){
  int jp[MXPAIR], jp2[MXPAIR];
  lgvn_data->esum = 0; 
  mol->evqdq = 0; 
  lgvn_data->ecor = 0; 
  lgvn_data->atomfs = malloc(sizeof(double) * mol->n_reg1);
  for(int i = 0; i < mol->n_reg1; i+=1){
    lgvn_data->atomfs[i] = 0; 
  }

  // elgvn = noniterative lgvn energy (using distance-dependent dielectric)
  // elgvni = lgvn energy from the 0th iteration 
  // elgwa = iterative lgvn energy 
  // evqdq = solute relaxation energy calculated from PCM charges 
  // ecor = correction for electron correlation effects

  // generate oshift points
  lgvn_data->oshift = generate_oshift(mol->ndxp);
  for(int i = 0; i < mol->ndxp; i+=1){
    lgvn_data->center_new = ran_shift(mol, lgvn_data->oshift, i);
    lgvnx(mol, lgvn_data, i);
    lgvn_data->esum += mol->elgvn;
    sci_lgvn(mol, lgvn_data, i, jp, jp2);
    lgvn_data->temp_elgvn[i] = mol->elgwa;
    lgvn_data->tdsl[i] = lgvn_data->tds;
    lgvn_data->tdsw_a = mol->phobsl * mol->fsurfa[i];
    lgvn_data->vatom_result = vatom_f(mol, lgvn_data);
    mol->evqdq += lgvn_data->vatom_result[0];
    lgvn_data->ecor += lgvn_data->vatom_result[1];
  }
  mol->elgvn = lgvn_data->esum / (double) mol->ndxp;
  mol->evqdq = -mol->evqdq / (double) mol->ndxp; 
  lgvn_data->ecor = lgvn_data->ecor / (double) mol->ndxp; 
  lgvn_data->elgvn_ave_result = elgvn_ave_f(mol, lgvn_data);
  mol->evdw = lgvn_data->elgvn_ave_result[0];
  mol->ephob = lgvn_data->elgvn_ave_result[1];
  mol->etds = lgvn_data->elgvn_ave_result[2];
  mol->elgwa = lgvn_data->elgvn_ave_result[3];
  lgvn_data->vbornx_result = vbornx_f(mol, lgvn_data);
  mol->ebw = lgvn_data->vbornx_result;
  //remember to free(lgvn_data->atomfs)
  return 0; 
}
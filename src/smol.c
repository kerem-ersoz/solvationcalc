#include "calc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ran2(int seed){
  // returns a uniform random number between 0.0 and 1.0.
  // set seed to any negative value to initialize or reinitialize 
  // the sequence with the seed number = -seed
  const int m = 714025;
  const int ia = 1366;
  const int ic = 150889;
  const double rm = 1.0 / m;
  static int ir[97], iy, i;
  double ran2; 
  if(seed < 0){
    seed = (ic-seed) % m;
    for(i = 0; i<97; i+=1){
      seed = (ia * seed + ic) % m;
      ir[i] = seed;
    }
    seed = (ia * seed + ic) % m;
    iy = seed; 
  }
  i = 1 + (97 * iy) / m;
  if(i > 97 || i < 1){
    printf("Error: Randomizer! i = %d\n", i);
    exit(1);
  }
  iy = ir[i];
  ran2 = iy * rm; 
  seed = (ia * seed + ic) % m;
  ir[i] = seed; 
  return ran2; 
}

double *generate_oshift(int ndxp){
  int iseed, idum, dumm;
  static double oshift[MXCENTER*3];
  // initialize the random number generator and 
  // generate random origin shifts for ndxp grids.
  iseed = -931;
  idum = 1; 
  dumm = (int)ran2(iseed); //initialize ran2
  for(int i = 0; i < (3 * ndxp); i+=1){
    oshift[i] = ran2(idum);
  }
  return oshift;
}

double *ran_shift(mol_t *mol, double *oshift, int i){
  double fact, dxp[3];
  static double center2[3];
  int iseed, idum, dumm; 
  fact = mol->drg_inner; 
  if(mol->rg_inner == 0) //only works for single-atom calculations
    fact = mol->drg; 
  dxp[0] = fact * (1.0 - 2.0 * oshift[3 * i - 2]);
  dxp[1] = fact * (1.0 - 2.0 * oshift[3 * i - 1]);
  dxp[2] = fact * (1.0 - 2.0 * oshift[3 * i]);
  if(i != 1){
    center2[0] = mol->pcenter[0] + dxp[0] + mol->dxp0[0];
    center2[1] = mol->pcenter[1] + dxp[1] + mol->dxp0[1];
    center2[2] = mol->pcenter[2] + dxp[2] + mol->dxp0[2];
  }
  else{
    center2[0] = mol->pcenter[0] + dxp[0] + mol->dxp0[0];
    center2[1] = mol->pcenter[1] + dxp[1] + mol->dxp0[1];
    center2[2] = mol->pcenter[2] + dxp[2] + mol->dxp0[2]; 
  }
  return center2;
}
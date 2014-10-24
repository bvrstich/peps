/**
 * @mainpage 
 * This is an implementation of an imaginary time evolution algorithm for the optimization of PEPS.
 * @author Brecht Verstichel
 * @date 25-03-2014
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);//dimension of the lattice: LxL
   int d = atoi(argv[2]);//physical dimension
   int D = atoi(argv[3]);//virtual dimension
   int D_aux = atoi(argv[4]);//virtual dimension

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   char dir_in[200];
   sprintf(dir_in,"/home/bright/bestanden/results/peps/%dx%d/D=%d/D_aux=%d/peps",L,L,D-1,4*(D-1)*(D-1));

   PEPS<double> peps;
   peps.load(dir_in);
   peps.sD(D);

   peps.grow_bond_dimension(D,0.01);
   peps.normalize();

   global::env.calc('A',peps);
   global::env.test();

   propagate::step(true,peps,n_steps);

   return 0;

}

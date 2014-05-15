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

   //initialize the dimensions
   Global::lat.set(L,L,d);
   Environment::init();

   int D_aux = 64;

   PEPS<double> peps(D);
   peps.normalize(D_aux);
   
   Environment::calc_env('A',peps,D_aux);
   Heisenberg::energy(peps);
   
/*
   peps.init_af();

   cout << peps.dot(peps,D_aux) << endl;

   peps.sD(D);

   double tau = 0.01;

   Trotter::init(tau);
*/
}

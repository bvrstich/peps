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
   PEPS<double>::lat.set(L,L,d);

   MPS<double> mps(L,D);

   mps.canonicalize(Right);

   for(int i = 0;i < D;++i)
      for(int k = 0;k < D;++k){

         double tmp = 0.0;

         for(int j = 0;j < D;++j)
            for(int s = 0;s < d;++s)
               tmp += mps[4](i,s,j)*mps[4](k,s,j);

         cout << i << "\t" << k << "\t" << tmp << endl;

      }

}

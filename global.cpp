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

namespace global{

   int D;

   int D_aux;

   int Lx;
   int Ly;

   int d;

   Random RN;

   DArray<2> I;

   Hamiltonian ham;

   Environment env;

   Trotter trot;

   /**
    * @param D_in virtual dimension of the trial
    * @param D_aux_in auxiliary dimension for peps contraction
    * @param d_in physical dimension
    * @param Lx_in x dimension of the square lattice
    * @param Ly_in y dimension of the square lattice
    * @param tau the time step of the imaginary time evolution
    */
   void init(int D_in,int D_aux_in,int d_in,int Lx_in,int Ly_in,double tau){

      Lx = Lx_in;
      Ly = Ly_in;

      d = d_in;

      D = D_in;

      D_aux = D_aux_in;

      //set the interaction
      ham.set_heisenberg(true);

      trot = Trotter(tau);

      //identity matrix
      I.resize(d,d);

      I = 0.0;

      I(0,0) = 1.0;
      I(1,1) = 1.0;

      //initialize/allocate the environment
      env = Environment(D_in,D_aux_in);

   }

   //!function which generates random complex numbers uniformly on a square of side 2 [(-1,1):(-1,1)]
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

}

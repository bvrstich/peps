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

/**
 * construct constructs a standard PEPS object
 * @param Lx_in input x dimension of the lattice
 * @param Ly_in input y dimension of the lattice
 * @param d_in physical dimension
 * @param D_in cutoff virtual dimension
 */
template<typename T,class Q>
PEPS<T,Q>::PEPS(int Lx_in,int Ly_in,int d_in,int D_in) : vector< QSTArray<T,5,Q> >(Lx_in*Ly_in) {

   Lx = Lx_in;
   Ly = Ly_in;
   d = d_in;
   D = D_in;

}

/**
 * empty destructor
 */
template<typename T,class Q>
PEPS<T,Q>::~PEPS(){ }

template PEPS<double,SpinQuantum>::PEPS(int,int,int,int);
template PEPS<double,SpinQuantum>::~PEPS();

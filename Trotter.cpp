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

DArray<3> Trotter::LO;

DArray<3> Trotter::RO;

double Trotter::tau;

/** 
 * initialize the operators
 * @param tau timestep
 */
void Trotter::init(double tau_in) {

   tau = tau_in;

   int d = Global::lat.gd();

   //first construct S_i.S_j on a d^2 x d^2 space
   DArray<2> Sij(d*d,d*d);
   Sij = 0.0;

   //basis: |dd> , |du> , |ud> , |uu>
   Sij(0,0) = 0.25;
   Sij(1,1) = -0.25;
   Sij(2,2) = -0.25;
   Sij(3,3) = 0.25;

   Sij(1,2) = 0.5;
   Sij(2,1) = 0.5;

   DArray<1> eig(d*d);

   lapack::syev(CblasRowMajor, 'V', 'U', d*d, Sij.data(), d*d, eig.data());

   //now construct exp(- tau Sij)
   DArray<2> ts_gate(d*d,d*d);

   for(int i = 0;i < d*d;++i)
      for(int j = 0;j < d*d;++j){

         ts_gate(i,j) = 0.0;

         for(int k = 0;k < d*d;++k)
            ts_gate(i,j) += exp( -tau * eig(k) ) * Sij(i,k) * Sij(j,k);

      }

   //now singular value decomposition over two sites: first reshape to |i><i'| x |j><j'| form
   DArray<2> tmp(d*d,d*d);
   tmp = 0.0;

   //new basis |d><d| |d><u| |u><d| |u><u|
   tmp(0,0) = ts_gate(0,0);
   tmp(3,3) = ts_gate(3,3);
   tmp(0,3) = ts_gate(1,1);
   tmp(3,0) = ts_gate(2,2);

   tmp(1,2) = ts_gate(1,2);
   tmp(2,1) = ts_gate(2,1);

   DArray<2> U(d*d,d*d);
   DArray<2> V(d*d,d*d);

   //now svd of tmp
   Gesvd('S','S',tmp,eig,U,V);

   //now put them correctly into left and right operators
   LO.resize(d,d*d,d);
   RO.resize(d,d*d,d);

   LO = 0.0;
   RO = 0.0;

   for(int s = 0;s < d;++s)
      for(int s_ = 0;s_ < d;++s_){

         int i = s*d + s_;

         for(int k = 0;k < d*d;++k){

            LO(s,k,s_) = U(i,k) * sqrt( eig(k) );
            RO(s,k,s_) = sqrt( eig(k) ) * V(k,i);

         }

      }

}

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &);

   void construct_reduced_tensor(char,char,const DArray<5> &,DArray<4> &,DArray<3> &);

   void get_X(DArray<4> &,DArray<3> &);

   void invert(DArray<2> &);

   void calc_N_eff(char,int,const DArray<2> &,const DArray<4> &,const DArray<2> &,const DArray<4> &,DArray<4> &);

   void calc_N_eff(char,int,int,const DArray<3> &,const DArray<4> &,const DArray<3> &,const DArray<4> &,DArray<4> &);

   void canonicalize(DArray<3> &X,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR);

   void update(int,DArray<3> &,DArray<3> &);

}

#endif

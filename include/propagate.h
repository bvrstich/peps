#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &,int);

   void construct_reduced_tensor(char,const DArray<5> &,DArray<4> &,DArray<3> &);

   void construct_double_layer(char,const DArray<4> &,DArray<5> &);

   void get_X(DArray<4> &,DArray<3> &);

   void invert(DArray<2> &);

   void calc_N_eff(char,int,const DArray<2> &,const DArray<4> &,const DArray<2> &,const DArray<4> &,DArray<4> &);

   void canonicalize(DArray<3> &X,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR);

   void update(int,DArray<3> &,DArray<3> &);

}

#endif

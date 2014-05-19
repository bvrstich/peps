#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &,int);

   void construct_reduced_tensor(char,DArray<5> &,DArray<4> &,DArray<3> &);

   void construct_double_layer(char,DArray<4> &,DArray<5> &);

   void get_X(DArray<4> &,DArray<3> &);

}

#endif

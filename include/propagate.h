#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &,int);

   void construct_reduced_tensor(char,DArray<5> &,DArray<4> &,DArray<3> &);

}

#endif

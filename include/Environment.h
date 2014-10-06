#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

template<typename T>
class PEPS;

template<typename T>
class MPS;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      static void init();

      static void calc_env(char,const PEPS<double> &,int);

      static void calc_env(char,int,const PEPS<double> &,int D_aux);

      static void test_env();

      static void construct_double_layer(char,const DArray<5> &peps,DArray<3> &dls);

      static void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls);

      static void construct_double_layer(char,const DArray<5> &peps,DArray<4> &dlo);

      static void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo);

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      static vector< MPS<double> > l;
      static vector< MPS<double> > r;
      static vector< MPS<double> > t;
      static vector< MPS<double> > b;

   private:

};

#endif

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

      Environment();

      Environment(int,int);

      //copy constructor
      Environment(const Environment &);

      //destructor
      virtual ~Environment();

      void calc(char,const PEPS<double> &,int);

      void calc(char,int,const PEPS<double> &,int D_aux);

      void test();

      void construct_double_layer(char,const DArray<5> &peps,DArray<3> &dls);

      void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls);

      void construct_double_layer(char,const DArray<5> &peps,DArray<4> &dlo);

      void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo);

      const MPS<double> &gl(int) const;
      MPS<double> &gl(int);

      const MPS<double> &gr(int) const;
      MPS<double> &gr(int);

      const MPS<double> &gt(int) const;
      MPS<double> &gt(int);

      const MPS<double> &gb(int) const;
      MPS<double> &gb(int);

      const int gD() const;
      const int gD_aux() const;

      void sD(int);
      void sD_aux(int);

   private:

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      vector< MPS<double> > l;
      vector< MPS<double> > r;
      vector< MPS<double> > t;
      vector< MPS<double> > b;

      //!regular bond dimension of peps
      int D;

      //!Auxiliary dimension, for the contraction
      int D_aux;

};

#endif

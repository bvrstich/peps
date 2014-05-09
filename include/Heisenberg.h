#ifndef HEISENBERG_H
#define HEISENBERG_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the expectation value of a Heisenberg hamiltonian between two PEPS.
 */
class Heisenberg {

   public:

      Heisenberg();
      
      void construct_environment(const PEPS<double> &,int D_aux);

      double local(const PEPS<double> &,const DArray<2> &);

      void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls);

      void construct_double_layer(char,const DArray<5> &peps,DArray<4> &dlo);

      void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo);

      double energy(const PEPS<double> &);
   
   private:

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      vector< MPS<double> > l;
      vector< MPS<double> > r;
      vector< MPS<double> > t;
      vector< MPS<double> > b;

      //!operators!
      DArray<2> Sp;
      DArray<2> Sm;
      DArray<2> Sz;
   
};

#endif

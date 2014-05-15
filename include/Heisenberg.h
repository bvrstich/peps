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

      static void init();

      static double local(const PEPS<double> &,const DArray<2> &);

      static void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls);

      static void construct_double_layer(char,const DArray<5> &peps,DArray<4> &dlo);

      static void construct_double_layer(char,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo);

      static double energy(const PEPS<double> &);
   
   private:

      //!operators!
      static DArray<2> Sp;
      static DArray<2> Sm;
      static DArray<2> Sz;
   
};

#endif

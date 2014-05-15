#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      static void init();

      static void calc_env(char,const PEPS<double> &,int D_aux);

      static void test_env();

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      static vector< MPS<double> > l;
      static vector< MPS<double> > r;
      static vector< MPS<double> > t;
      static vector< MPS<double> > b;

   private:

};

#endif

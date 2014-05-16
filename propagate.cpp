#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

namespace propagate {

   /**
    * propagate the peps one imaginary time step
    * @param peps the PEPS to be propagated
    * @param D_aux auxiliary dimension for the contractions. Determines the accuracy of the effective environment.
    */
   void step(PEPS<double> &peps,int D_aux){

      int Lx = Global::lat.gLx();
      int Ly = Global::lat.gLy();

      int d = Global::lat.gd();

      // ##################################################### //
      // first propagate applying the gates from bottom to top //
      // ##################################################### //

      //construct the top environment:
      Environment::calc_env('T',peps,D_aux);

      //for the bottom row: construct the right operator
      vector< DArray<2> > R(Lx - 2);

      //for this we need the construct the 'bottom environment' for the bottom row:
      Environment::calc_env('B',0,peps,D_aux);

      //first the rightmost operator
      DArray<4> tmp;
      DArray<3> I;

      //tmp comes out index (t,b)
      Contract(1.0,Environment::t[0][Lx - 1],shape(1),Environment::b[0][Lx - 1],shape(1),0.0,tmp);

      //reshape tmp to a 2-index array
      R[Lx - 3] = tmp.reshape_clear(shape(Environment::t[0][Lx - 1].shape(0),Environment::b[0][Lx - 1].shape(0)));

      //now construct the rest
      for(int col = Lx - 2;col > 1;--col){

         I.clear();
         Contract(1.0,Environment::t[0][col],shape(2),R[col-1],shape(0),0.0,I);

         Contract(1.0,I,shape(1,2),Environment::b[0][col],shape(1,2),0.0,R[col-2]);

      }

      cout << R[0] << endl;

   }

}

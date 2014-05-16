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
      int D = peps.gD();

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

      //now construct the reduced tensors of the pair to propagate
      DArray<4> Q;
      DArray<3> a_L;

      construct_reduced_tensor('L',peps(0,0),Q,a_L);

      //make a 'double layer' object out of Q for contraction with environment
      DArray<5> dlQ;
      construct_double_layer(Q,dlQ);

   }

   /**
    * construct the left or right reduced tensor form of a peps element by performing QR or LQ decomposition
    * @param L == left, R == right
    */
   void construct_reduced_tensor(char option,DArray<5> &peps,DArray<4> &Q,DArray<3> &red){

      if(option == 'L'){

         DArray<5> tmp;
         Permute(peps,shape(0,1,3,2,4),tmp);

         int nrows = tmp.shape(0) * tmp.shape(1) * tmp.shape(2);
         int ncols = tmp.shape(3) * tmp.shape(4);

         int min = std::min(nrows,ncols);

         double* tau = new double [min];

         lapack::geqrf(CblasRowMajor,nrows,ncols, tmp.data(), ncols, tau);

         red.resize(shape(min,tmp.shape(3),tmp.shape(4)));
         Q.resize(tmp.shape(0),tmp.shape(1),tmp.shape(2),min);

         //r is in the upper diagonal part of tmp on exit of geqrf:
         for(int i = 0;i < min;++i)
            for(int j = i;j < ncols;++j)
               red.data()[i*ncols + j] = tmp.data()[i*ncols + j];

         //get the input for the Q construction function: lower diagonal part of the matrix
         for(int i = 0;i < nrows;++i)
            for(int j = 0;j < min;++j)
               Q.data()[i*min + j] = tmp.data()[i*ncols + j];

         //now get the Q matrix out
         if(nrows < ncols)
            lapack::orgqr(CblasRowMajor, nrows, nrows, min,Q.data(), nrows, tau);
         else
            lapack::orgqr(CblasRowMajor, nrows, ncols, min,Q.data(), ncols, tau);

         delete [] tau;

      }
      else{//LQ

         DArray<5> tmp;
         Permute(peps,shape(0,2,1,3,4),tmp);

         cout << tmp << endl;

         int nrows = tmp.shape(0) * tmp.shape(1);
         int ncols = tmp.shape(2) * tmp.shape(3) * tmp.shape(4);

         int min = std::min(nrows,ncols);

         double* tau = new double [min];

         lapack::gelqf(CblasRowMajor,nrows,ncols,tmp.data(),ncols,tau);

         red.resize(shape(tmp.shape(0),tmp.shape(1),min));
         Q.resize(min,tmp.shape(2),tmp.shape(3),tmp.shape(4));

         //l is in the lower diagonal part of tmp on exit of geqrf:
         for(int j = 0;j < ncols;++j)
            for(int i = j;i < nrows;++i)
               red.data()[i*min + j] = tmp.data()[i*ncols + j];

         //get the input for the Q construction function: upper diagonal part of the matrix
         for(int i = 0;i < min;++i)
            for(int j = 0;j < ncols;++j)
               Q.data()[i*ncols + j] = tmp.data()[i*ncols + j];

         //now get the Q matrix out
         if(nrows > ncols)
            lapack::orglq(CblasRowMajor,ncols,ncols,min,Q.data(),ncols,tau);
         else
            lapack::orglq(CblasRowMajor,nrows,ncols,min,Q.data(),ncols,tau);

         delete [] tau;

      }

   }

   /**
    * construct a double layer object out of a Q coming from a reduced tensor construction.
    * keep one leg, the one pointing to the reduced tensor, not doubled.
    * @param Q input object
    * @param dlQ output object
    */
   void construct_double_layer(DArray<4> &Q,DArray<5> &dlQ){

      //first outer product of Q
      DArray<8> tmp;
      Ger(1.0,Q,Q,tmp);

      DArray<8> reorder;

      Permute(tmp,shape(0,4,1,5,2,6,3,7),reorder);

      //and move it to dlQ
      int d0 = reorder.shape(0) * reorder.shape(1);
      int d1 = reorder.shape(2) * reorder.shape(3);
      int d2 = reorder.shape(4) * reorder.shape(5);
      int d3 = reorder.shape(6);
      int d4 = reorder.shape(7);

      dlQ = reorder.reshape_clear(shape(d0,d1,d2,d3,d4));

   }

}
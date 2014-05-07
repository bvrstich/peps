#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

/** 
 * constructor
 */
Heisenberg::Heisenberg(){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

}

/**
 * construct the enviroment mps's for the input PEPS
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Heisenberg::construct_environment(const PEPS<double> &peps,int D_aux){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   //construct bottom layer
   b[0] = MPS<double>('b',peps,peps);

   for(int i = 1;i < Ly - 1;++i){

      //i'th row as MPO
      MPO<double> mpo('H',i,peps,peps);

      MPS<double> tmp(b[i - 1]);

      //apply to form MPS with bond dimension D^4
      tmp.gemv('L',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      b[i].resize(Lx);
      b[i].compress(D_aux,tmp,5);

   }

   //then construct top layer
   t[Ly - 2] = MPS<double>('t',peps,peps);

   for(int i = Ly - 2;i > 0;--i){

      //i'th row as MPO
      MPO<double> mpo('H',i,peps,peps);

      //apply to form MPS with bond dimension D^4
      MPS<double> tmp(t[i]);

      tmp.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      t[i - 1].resize(Lx);
      t[i - 1].compress(D_aux,tmp,5);

   }

   //then left layer
   l[0] = MPS<double>('l',peps,peps);

   for(int i = 1;i < Lx - 1;++i){

      //i'th col as MPO
      MPO<double> mpo('V',i,peps,peps);

      MPS<double> tmp(l[i - 1]);

      //apply to form MPS with bond dimension D^4
      tmp.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      l[i].resize(Ly);
      l[i].compress(D_aux,tmp,5);

   }

   //finally construct right layer
   r[Lx - 2] = MPS<double>('r',peps,peps);

   for(int i = Lx - 2;i > 0;--i){

      //i'th row as MPO
      MPO<double> mpo('V',i,peps,peps);

      //apply to form MPS with bond dimension D^4
      MPS<double> tmp(r[i]);

      tmp.gemv('L',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      r[i - 1].resize(Ly);
      r[i - 1].compress(D_aux,tmp,5);

   }

   for(int i = 0;i < Ly - 1;++i)
      cout << i << "\t" << b[i].dot(t[i]) << endl;

   for(int i = 0;i < Lx - 1;++i)
      cout << i << "\t" << r[i].dot(l[i]) << endl;

}

/**
 * evaluate the expectation value of a local operator O = \sum_i O_i
 * @param O DArray<2> object of the physical dimension of a site.
 * beware, the environments have to be constructed beforehand!
 */
double Heisenberg::local(const PEPS<double> &peps,const DArray<2> &O){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   //from bottom to top: contract in mps/mpo fashion

   //bottom row: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<2> > R(Lx - 1);

   //first the rightmost operator
   DArray<4> tmp;
   DArray<3> I;

   //tmp comes out index (t,b)
   Contract(1.0,t[0][Lx - 1],shape(1),b[0][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),b[0][Lx - 1].shape(0)));

   //now construct the rest
   for(int i = Lx - 2;i > 0;--i){

      I.clear();
      Contract(1.0,t[0][i],shape(2),R[i],shape(0),0.0,I);

      Contract(1.0,I,shape(1,2),b[0][i],shape(1,2),0.0,R[i - 1]);

   }

   //now sweep from left to right to get the expectation value of the local operator on the bottom row
   double val = 0.0;

   //construct the double layer object of lowest row peps with operator O in between
   DArray<3> dl;
   construct_double_layer(peps(0,0),O,dl);

   //we will need left renormalized operators as well
   vector< DArray<2> > L(Lx - 1);

   //tmp comes out index (t,b)
   Contract(1.0,t[0][0],shape(1),dl,shape(1),0.0,tmp);

   L[0] = tmp.reshape_clear(shape(t[0][0].shape(2),dl.shape(2)));

   //first value
   val += Dot(L[0],R[0]);

   //construct left renormalized operator
   Contract(1.0,t[0][0],shape(1),b[0][0],shape(1),0.0,tmp);

   L[0] = tmp.reshape_clear(shape(t[0][0].shape(2),b[0][0].shape(2)));

   //middle of the chain:
   for(int c = 1;c < Lx-1;++c){

      construct_double_layer(peps(0,c),O,dl);

      I.clear();
      Contract(1.0,t[0][c],shape(2),R[c],shape(0),0.0,I);

      Contract(1.0,I,shape(1,2),dl,shape(1,2),0.0,R[c - 1]);
 
      cout << L[c - 1].shape() << endl;
      cout << R[c - 1].shape() << endl;
     
      val += Dot(L[c - 1],R[c - 1]);

      //construct left renormalized operator
      I.clear();
      Contract(1.0,L[c - 1],shape(0),t[0][c],shape(0),0.0,I);

      Contract(1.0,I,shape(0,1),b[0][c],shape(0,1),0.0,L[c]);

   }

   return val;

}

/**
 * construct a double layer MPS object from a PEPS taking the expectation value of an operator O between the physical indices
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dl output object
 */
void Heisenberg::construct_double_layer(const DArray<5> &peps,const DArray<2> &O,DArray<3> &dl){

   enum {i,j,k,l,m,n,o,p,q};

   DArray<5> tmp;
   Contract(1.0,peps,shape(i,j,p,k,l),O,shape(p,q),0.0,tmp,shape(i,j,q,k,l));

   DArray<8> tmp2;
   Contract(1.0,tmp,shape(i,j,k,l,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(i,n,j,o,l,p,m,q));

   int DL = tmp2.shape(0) * tmp2.shape(1);
   int d_phys = tmp2.shape(2) * tmp2.shape(3) * tmp2.shape(4) * tmp2.shape(5);
   int DR = tmp2.shape(6) * tmp2.shape(7);

   dl = tmp2.reshape_clear(shape(DL,d_phys,DR));

}

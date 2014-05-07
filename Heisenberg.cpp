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

   cout << endl;
   cout << "FROM BOTTOM TO TOP" << endl;
   cout << endl;
   for(int i = 0;i < Ly - 1;++i)
      cout << i << "\t" << b[i].dot(t[i]) << endl;

   cout << endl;
   cout << "FROM LEFT TO RIGHT" << endl;
   cout << endl;
   for(int i = 0;i < Lx - 1;++i)
      cout << i << "\t" << r[i].dot(l[i]) << endl;
   cout << endl;

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

   // -- (1) -- || bottom row: similar to overlap calculation

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
   for(int c = Lx - 2;c > 0;--c){

      I.clear();
      Contract(1.0,t[0][c],shape(2),R[c],shape(0),0.0,I);

      Contract(1.0,I,shape(1,2),b[0][c],shape(1,2),0.0,R[c - 1]);

   }

   //now sweep from left to right to get the expectation value of the local operator on the bottom row
   double val = 0.0;

   //construct the double layer object of lowest row peps with operator O in between
   DArray<3> dls;
   construct_double_layer(peps(0,0),O,dls);

   //we will need left renormalized operators as well
   DArray<2> L;

   //tmp comes out index (t,b)
   Contract(1.0,t[0][0],shape(1),dls,shape(1),0.0,tmp);

   L = tmp.reshape_clear(shape(t[0][0].shape(2),dls.shape(2)));

   //first value
   val += Dot(L,R[0]);

   //construct left renormalized operator
   Contract(1.0,t[0][0],shape(1),b[0][0],shape(1),0.0,tmp);

   L = tmp.reshape_clear(shape(t[0][0].shape(2),b[0][0].shape(2)));

   //middle of the chain:
   for(int c = 1;c < Lx-1;++c){

      construct_double_layer(peps(0,c),O,dls);

      I.clear();
      Contract(1.0,t[0][c],shape(2),R[c],shape(0),0.0,I);

      Contract(1.0,I,shape(1,2),dls,shape(1,2),0.0,R[c - 1]);
 
      val += Dot(L,R[c - 1]);

      //construct left renormalized operator
      I.clear();
      Contract(1.0,L,shape(0),t[0][c],shape(0),0.0,I);

      L.clear();
      Contract(1.0,I,shape(0,1),b[0][c],shape(0,1),0.0,L);

   }

   //last site of bottom row
   construct_double_layer(peps(0,Lx-1),O,dls);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[0][Lx - 1],shape(1),dls,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),dls.shape(0)));

   val += Dot(L,R[Lx-2]);

   // -- (2) -- || middle rows: similar to MPO/MPS expectation value
   vector< DArray<3> > RO(Lx - 1);
   DArray<3> LO;

   DArray<4> dlo;

   for(int r = 1;r < Ly - 1;++r){

      //first create right renormalized operator

      //first site make double layer object from peps
      Heisenberg::construct_double_layer(peps(r,Lx-1),dlo);

      //paste top environment on
      DArray<5> tmp5;
      Contract(1.0,t[r][Lx - 1],shape(1),dlo,shape(1),0.0,tmp5);

      //then bottom enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),b[r-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 2] = tmp6.reshape_clear(shape(t[r][Lx - 1].shape(0),dlo.shape(0),b[r-1][Lx - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int c = Lx-2;c > 0;--c){

         I4.clear();
         Contract(1.0,t[r][c],shape(2),RO[c],shape(0),0.0,I4);

         enum {i,j,k,l,m,n};

         Heisenberg::construct_double_layer(peps(r,c),dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,l),dlo,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,l));

         RO[c-1].clear();
         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(1,2),0.0,RO[c-1]);

      }

      //expectation value of operator on first site

      //first site make double layer object from peps
      Heisenberg::construct_double_layer(peps(r,0),O,dlo);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[r][0],shape(1),dlo,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[r-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LO = tmp6.reshape_clear(shape(t[r][0].shape(2),dlo.shape(3),b[r-1][0].shape(2)));

      //get expectation value for operator on this site
      val += Dot(LO,RO[0]);

      //construct left renormalized operator
      Heisenberg::construct_double_layer(peps(r,0),dlo);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[r][0],shape(1),dlo,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[r-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LO = tmp6.reshape_clear(shape(t[r][0].shape(2),dlo.shape(3),b[r-1][0].shape(2)));

      //middle sites
      for(int c = 1;c < Lx-1;++c){

         I4.clear();
         Contract(1.0,t[r][c],shape(2),RO[c],shape(0),0.0,I4);

         enum {i,j,k,l,m,n};

         Heisenberg::construct_double_layer(peps(r,c),O,dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,l),dlo,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,l));

         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(1,2),0.0,RO[c-1]);

         val += Dot(LO,RO[c - 1]);

         //construct left renormalized operator
         I4.clear();
         Contract(1.0,t[r][c],shape(0),LO,shape(0),0.0,I4);

         Heisenberg::construct_double_layer(peps(r,c),dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,l),dlo,shape(k,i,m,n),0.0,I4bis,shape(j,n,l,m));

         LO.clear();
         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(0,1),0.0,LO);

      }

      //last site: first make double layer with local operator
      Heisenberg::construct_double_layer(peps(r,Lx-1),O,dlo);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[r][Lx - 1],shape(1),dlo,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[r-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 2] = tmp6.reshape_clear(shape(t[r][Lx - 1].shape(0),dlo.shape(0),b[r-1][Lx - 1].shape(0)));

      //get expectation value
      val += Dot(LO,RO[Lx-2]);

   }//end of the middle rows!

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[Ly-2][Lx - 1],shape(1),b[Ly-2][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp.reshape_clear(shape(t[Ly-2][Lx - 1].shape(0),b[Ly-2][Lx - 1].shape(0)));

   //now construct the rest
   for(int c = Lx - 2;c > 0;--c){

      I.clear();
      Contract(1.0,t[Ly-2][c],shape(2),R[c],shape(0),0.0,I);

      R[c-1].clear();
      Contract(1.0,I,shape(1,2),b[Ly-2][c],shape(1,2),0.0,R[c - 1]);

   }

   //construct the double layer object of top row peps with operator O in between
   construct_double_layer(peps(Ly-1,0),O,dls);

   //tmp comes out index (t,b)
   Contract(1.0,dls,shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   L = tmp.reshape_clear(shape(dls.shape(2),b[Ly-2][0].shape(2)));

   //first value
   val += Dot(L,R[0]);

   //construct left renormalized operator
   Contract(1.0,t[Ly-2][0],shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   L = tmp.reshape_clear(shape(t[Ly-2][0].shape(2),b[Ly-2][0].shape(2)));

   //middle of the chain:
   for(int c = 1;c < Lx-1;++c){

      construct_double_layer(peps(Ly-1,c),O,dls);

      I.clear();
      Contract(1.0,dls,shape(2),R[c],shape(0),0.0,I);

      R[c-1].clear();
      Contract(1.0,I,shape(1,2),b[Ly-2][c],shape(1,2),0.0,R[c - 1]);
 
      val += Dot(L,R[c - 1]);

      //construct left renormalized operator
      I.clear();
      Contract(1.0,L,shape(0),t[Ly-2][c],shape(0),0.0,I);

      L.clear();
      Contract(1.0,I,shape(0,1),b[Ly-2][c],shape(0,1),0.0,L);

   }

   //last site of top row
   construct_double_layer(peps(Ly-1,Lx-1),O,dls);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,dls,shape(1),b[Ly-2][Lx-1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp.reshape_clear(shape(dls.shape(0),b[Ly-2][Lx-1].shape(0)));

   val += Dot(L,R[Lx-2]);

   return val;

}

/**
 * construct a double layer MPS object from a PEPS taking the expectation value of an operator O between the physical indices
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dls output object
 */
void Heisenberg::construct_double_layer(const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls){

   enum {i,j,k,l,m,n,o,p,q};

   DArray<5> tmp;
   Contract(1.0,peps,shape(i,j,p,k,l),O,shape(p,q),0.0,tmp,shape(i,j,q,k,l));

   DArray<8> tmp2;
   Contract(1.0,tmp,shape(i,j,k,l,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(i,n,j,o,l,p,m,q));

   int DL = tmp2.shape(0) * tmp2.shape(1);
   int d_phys = tmp2.shape(2) * tmp2.shape(3) * tmp2.shape(4) * tmp2.shape(5);
   int DR = tmp2.shape(6) * tmp2.shape(7);

   dls = tmp2.reshape_clear(shape(DL,d_phys,DR));

}

/**
 * construct a double layer MPO object from a PEPS
 * @param peps THe input PEPS elements
 * @param dlo output object
 */
void Heisenberg::construct_double_layer(const DArray<5> &peps,DArray<4> &dlo){

   enum {i,j,k,l,m,n,o,p,q};

   DArray<8> tmp;
   Contract(1.0,peps,shape(i,j,k,l,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,l,p,m,q));

   int DL = tmp.shape(0) * tmp.shape(1);
   int DU = tmp.shape(2) * tmp.shape(3);
   int DD = tmp.shape(4) * tmp.shape(5);
   int DR = tmp.shape(6) * tmp.shape(7);

   dlo = tmp.reshape_clear(shape(DL,DU,DD,DR));

}

/**
 * construct a double layer MPO object from a PEPS, using operator O in between
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dlo output object
 */
void Heisenberg::construct_double_layer(const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo){

   enum {i,j,k,l,m,n,o,p,q};

   DArray<5> tmp;
   Contract(1.0,peps,shape(i,j,p,k,l),O,shape(p,q),0.0,tmp,shape(i,j,q,k,l));

   DArray<8> tmp2;
   Contract(1.0,tmp,shape(i,j,k,l,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(i,n,j,o,l,p,m,q));

   int DL = tmp2.shape(0) * tmp2.shape(1);
   int DU = tmp2.shape(2) * tmp2.shape(3);
   int DD = tmp2.shape(4) * tmp2.shape(5);
   int DR = tmp2.shape(6) * tmp2.shape(7);

   dlo = tmp2.reshape_clear(shape(DL,DU,DD,DR));

}

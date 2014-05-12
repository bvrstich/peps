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

   int d = Global::lat.gd();

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

   //init the operators
   Sp.resize(d,d);
   Sm.resize(d,d);
   Sz.resize(d,d);

   Sp = 0.0;
   Sm = 0.0;
   Sz = 0.0;

   Sp(1,0) = 1.0;

   Sm(0,1) = 1.0;

   Sz(0,0) = -0.5;
   Sz(1,1) = 0.5;

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
      tmp.gemv('L',mpo);

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

      tmp.gemv('U',mpo);

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
 * @param peps the input PEPS 
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
   construct_double_layer('H',peps(0,0),O,dls);

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

      construct_double_layer('H',peps(0,c),O,dls);

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
   construct_double_layer('H',peps(0,Lx-1),O,dls);

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
      Heisenberg::construct_double_layer('H',peps(r,Lx-1),dlo);

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

         enum {i,j,k,m,n,o};

         Heisenberg::construct_double_layer('H',peps(r,c),dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlo,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[c-1].clear();
         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(1,2),0.0,RO[c-1]);

      }

      //expectation value of operator on first site

      //first site make double layer object from peps
      Heisenberg::construct_double_layer('H',peps(r,0),O,dlo);

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
      Heisenberg::construct_double_layer('H',peps(r,0),dlo);

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

         enum {i,j,k,o,m,n};

         Heisenberg::construct_double_layer('H',peps(r,c),O,dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlo,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(1,2),0.0,RO[c-1]);

         val += Dot(LO,RO[c - 1]);

         //construct left renormalized operator
         I4.clear();
         Contract(1.0,t[r][c],shape(0),LO,shape(0),0.0,I4);

         Heisenberg::construct_double_layer('H',peps(r,c),dlo);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlo,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LO.clear();
         Contract(1.0,I4bis,shape(2,3),b[r-1][c],shape(0,1),0.0,LO);

      }

      //last site: first make double layer with local operator
      Heisenberg::construct_double_layer('H',peps(r,Lx-1),O,dlo);

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
   construct_double_layer('H',peps(Ly-1,0),O,dls);

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

      construct_double_layer('H',peps(Ly-1,c),O,dls);

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
   construct_double_layer('H',peps(Ly-1,Lx-1),O,dls);

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
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dls output object
 */
void Heisenberg::construct_double_layer(char option,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls){

   if(option == 'H'){

      enum {i,j,k,m,n,o,p,q,s};

      DArray<5> tmp;
      Contract(1.0,peps,shape(i,j,p,k,s),O,shape(p,q),0.0,tmp,shape(i,j,q,k,s));

      DArray<8> tmp2;
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(i,n,j,o,s,p,m,q));

      int DL = tmp2.shape(0) * tmp2.shape(1);
      int d_phys = tmp2.shape(2) * tmp2.shape(3) * tmp2.shape(4) * tmp2.shape(5);
      int DR = tmp2.shape(6) * tmp2.shape(7);

      dls = tmp2.reshape_clear(shape(DL,d_phys,DR));

   }
   else{//V

      enum {i,j,k,s,m,n,o,p,q};

      DArray<5> tmp;
      Contract(1.0,peps,shape(i,j,p,k,s),O,shape(p,q),0.0,tmp,shape(i,j,q,k,s));

      DArray<8> tmp2;
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(s,p,i,n,m,q,j,o));

      int DL = tmp2.shape(0) * tmp2.shape(1);
      int d_phys = tmp2.shape(2) * tmp2.shape(3) * tmp2.shape(4) * tmp2.shape(5);
      int DR = tmp2.shape(6) * tmp2.shape(7);

      dls = tmp2.reshape_clear(shape(DL,d_phys,DR));

   }

}

/**
 * construct a double layer MPO object from a PEPS
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param dlo output object
 */
void Heisenberg::construct_double_layer(char option,const DArray<5> &peps,DArray<4> &dlo){

   if(option == 'H'){

      enum {i,j,k,s,m,n,o,p,q};

      DArray<8> tmp;
      Contract(1.0,peps,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,s,p,m,q));

      int DL = tmp.shape(0) * tmp.shape(1);
      int DU = tmp.shape(2) * tmp.shape(3);
      int DD = tmp.shape(4) * tmp.shape(5);
      int DR = tmp.shape(6) * tmp.shape(7);

      dlo = tmp.reshape_clear(shape(DL,DU,DD,DR));

   }
   else{//V

      enum {i,j,k,s,m,n,o,p,q};

      DArray<8> tmp;
      Contract(1.0,peps,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(s,p,m,q,i,n,j,o));

      int DL = tmp.shape(0) * tmp.shape(1);
      int DU = tmp.shape(2) * tmp.shape(3);
      int DD = tmp.shape(4) * tmp.shape(5);
      int DR = tmp.shape(6) * tmp.shape(7);

      dlo = tmp.reshape_clear(shape(DL,DU,DD,DR));

   }

}

/**
 * construct a double layer MPO object from a PEPS, using operator O in between
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dlo output object
 */
void Heisenberg::construct_double_layer(char option,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo){

   if(option == 'H'){

      enum {i,j,k,s,m,n,o,p,q};

      DArray<5> tmp;
      Contract(1.0,peps,shape(i,j,p,k,s),O,shape(p,q),0.0,tmp,shape(i,j,q,k,s));

      DArray<8> tmp2;
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(i,n,j,o,s,p,m,q));

      int DL = tmp2.shape(0) * tmp2.shape(1);
      int DU = tmp2.shape(2) * tmp2.shape(3);
      int DD = tmp2.shape(4) * tmp2.shape(5);
      int DR = tmp2.shape(6) * tmp2.shape(7);

      dlo = tmp2.reshape_clear(shape(DL,DU,DD,DR));

   }
   else{//V

      enum {i,j,k,s,m,n,o,p,q};

      DArray<5> tmp;
      Contract(1.0,peps,shape(i,j,p,k,s),O,shape(p,q),0.0,tmp,shape(i,j,q,k,s));

      DArray<8> tmp2;
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(s,p,m,q,i,n,j,o));

      int DL = tmp2.shape(0) * tmp2.shape(1);
      int DU = tmp2.shape(2) * tmp2.shape(3);
      int DD = tmp2.shape(4) * tmp2.shape(5);
      int DR = tmp2.shape(6) * tmp2.shape(7);

      dlo = tmp2.reshape_clear(shape(DL,DU,DD,DR));

   }

}

/**
 * evaluate the expectation value of the energy for the nn-Heisenberg model
 * @param peps the input PEPS 
 * beware, the environments have to be constructed beforehand!
 */
double Heisenberg::energy(const PEPS<double> &peps){

   // ---- || evaluate the energy in an MPO/MPS manner, first from bottom to top, then left to right || ----
   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   // #################################################################
   // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || bottom row: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<2> > R(Lx - 2);

   //first the rightmost operator
   DArray<4> tmp;
   DArray<3> I;

   //tmp comes out index (t,b)
   Contract(1.0,t[0][Lx - 1],shape(1),b[0][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),b[0][Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 1;--col){

      I.clear();
      Contract(1.0,t[0][col],shape(2),R[col-1],shape(0),0.0,I);

      Contract(1.0,I,shape(1,2),b[0][col],shape(1,2),0.0,R[col-2]);

   }

   //4 left going operators: S+, S-, Sz, and 1
   DArray<2> Lp;
   DArray<2> Lm;
   DArray<2> Lz;
   DArray<2> Lu;

   DArray<3> dlsm;
   DArray<3> dlsp;
   DArray<3> dlsz;

   //first S+
   construct_double_layer('H',peps(0,0),Sp,dlsp);

   //tmp comes out index (t,b)
   Contract(1.0,t[0][0],shape(1),dlsp,shape(1),0.0,tmp);

   Lp = tmp.reshape_clear(shape(t[0][0].shape(2),dlsp.shape(2)));

   //then S-
   construct_double_layer('H',peps(0,0),Sm,dlsm);

   //tmp comes out index (t,b)
   Contract(1.0,t[0][0],shape(1),dlsm,shape(1),0.0,tmp);

   Lm = tmp.reshape_clear(shape(t[0][0].shape(2),dlsm.shape(2)));

   //then Sz 
   construct_double_layer('H',peps(0,0),Sz,dlsz);

   //tmp comes out index (t,b)
   Contract(1.0,t[0][0],shape(1),dlsz,shape(1),0.0,tmp);

   Lz = tmp.reshape_clear(shape(t[0][0].shape(2),dlsz.shape(2)));

   //and finally unity
   Contract(1.0,t[0][0],shape(1),b[0][0],shape(1),0.0,tmp);

   Lu = tmp.reshape_clear(shape(t[0][0].shape(2),b[0][0].shape(2)));

   double val = 0.0;

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //first close down the +,- and z terms from the previous site

      //construct the right intermediate contraction (paste top to right)
      I.clear();
      Contract(1.0,t[0][col],shape(2),R[col - 1],shape(0),0.0,I);

      // 1) construct Sm double layer
      construct_double_layer('H',peps(0,col),Sm,dlsm);

      R[col-1].clear();
      Contract(1.0,I,shape(1,2),dlsm,shape(1,2),0.0,R[col - 1]);

      //contract with left S+
      val += 0.5 * Dot(Lp,R[col - 1]);

      // 2) then construct Sp double layer
      construct_double_layer('H',peps(0,col),Sp,dlsp);

      R[col-1].clear();
      Contract(1.0,I,shape(1,2),dlsp,shape(1,2),0.0,R[col - 1]);

      //contract with left S-
      val += 0.5 * Dot(Lm,R[col - 1]);

      // 3) then construct Sz double layer
      construct_double_layer('H',peps(0,col),Sz,dlsz);

      R[col-1].clear();
      Contract(1.0,I,shape(1,2),dlsz,shape(1,2),0.0,R[col - 1]);

      //contract with left Sz
      val += Dot(Lz,R[col - 1]);

      //construct left renormalized operators for next site: first paste top to Left unity
      I.clear();
      Contract(1.0,Lu,shape(0),t[0][col],shape(0),0.0,I);

      // 1) construct new Sm left operator
      Lm.clear();
      Contract(1.0,I,shape(0,1),dlsm,shape(0,1),0.0,Lm);

      // 2) construct new Sp left operator
      Lp.clear();
      Contract(1.0,I,shape(0,1),dlsp,shape(0,1),0.0,Lp);

      // 3) construct new Sz left operator
      Lz.clear();
      Contract(1.0,I,shape(0,1),dlsz,shape(0,1),0.0,Lz);

      // 4) finally construct new unity on the left
      Lu.clear();
      Contract(1.0,I,shape(0,1),b[0][col],shape(0,1),0.0,Lu);

   }

   //last site of bottom row:close down the left +,- and z

   //1) Sm to close down Lp
   construct_double_layer('H',peps(0,Lx-1),Sm,dlsm);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[0][Lx - 1],shape(1),dlsm,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),dlsm.shape(0)));

   val += 0.5 * Dot(Lp,R[Lx-3]);

   //2) Sp to close down Lm
   construct_double_layer('H',peps(0,Lx-1),Sp,dlsp);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[0][Lx - 1],shape(1),dlsp,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),dlsp.shape(0)));

   val += 0.5 * Dot(Lm,R[Lx-3]);

   //3) Sz to close down Lz
   construct_double_layer('H',peps(0,Lx-1),Sz,dlsz);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[0][Lx - 1],shape(1),dlsz,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(t[0][Lx - 1].shape(0),dlsz.shape(0)));

   val += Dot(Lz,R[Lx-3]);

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< DArray<3> > RO(Lx - 2);

   //4 left renormalized operators needed
   DArray<3> LOp;
   DArray<3> LOm;
   DArray<3> LOz;
   DArray<3> LOu;

   //double layer objects
   DArray<4> dlop;
   DArray<4> dlom;
   DArray<4> dloz;
   DArray<4> dlou;

   for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operator

      //first site make double layer object from peps
      Heisenberg::construct_double_layer('H',peps(row,Lx-1),dlou);

      //paste top environment on
      DArray<5> tmp5;
      Contract(1.0,t[row][Lx - 1],shape(1),dlou,shape(1),0.0,tmp5);

      //then bottom enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(t[row][Lx - 1].shape(0),dlou.shape(0),b[row-1][Lx - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int col = Lx-2;col > 1;--col){

         I4.clear();
         Contract(1.0,t[row][col],shape(2),RO[col-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         Heisenberg::construct_double_layer('H',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-2].clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(1,2),0.0,RO[col-2]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) S+ -- make double layer object from peps with Sp
      Heisenberg::construct_double_layer('H',peps(row,0),Sp,dlop);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][0],shape(1),dlop,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOp = tmp6.reshape_clear(shape(t[row][0].shape(2),dlop.shape(3),b[row-1][0].shape(2)));

      // 2) S- -- make double layer object from peps with Sm
      Heisenberg::construct_double_layer('H',peps(row,0),Sm,dlom);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][0],shape(1),dlom,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOm = tmp6.reshape_clear(shape(t[row][0].shape(2),dlom.shape(3),b[row-1][0].shape(2)));

      // 3) Sz -- make double layer object from peps with Sz
      Heisenberg::construct_double_layer('H',peps(row,0),Sz,dloz);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][0],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOz = tmp6.reshape_clear(shape(t[row][0].shape(2),dlom.shape(3),b[row-1][0].shape(2)));

      // 4) 1 -- finally construct left renormalized operator with unity
      Heisenberg::construct_double_layer('H',peps(row,0),dlou);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][0],shape(1),dlou,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOu = tmp6.reshape_clear(shape(t[row][0].shape(2),dlou.shape(3),b[row-1][0].shape(2)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int col = 1;col < Lx - 1;++col){

         //first add top to the right side, put it in I4
         I4.clear();
         Contract(1.0,t[row][col],shape(2),RO[col-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         //1) close down LOp with Sm
         Heisenberg::construct_double_layer('H',peps(row,col),Sm,dlom);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += 0.5 * Dot(LOp,RO[col-1]);

         //2) close down LOm with Sp
         Heisenberg::construct_double_layer('H',peps(row,col),Sp,dlop);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += 0.5 * Dot(LOm,RO[col-1]);

         //3) finally close down LOz with Sz
         Heisenberg::construct_double_layer('H',peps(row,col),Sz,dloz);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += Dot(LOz,RO[col-1]);

         // now construct the new left going renormalized operators
         //first attach top to left unity
         I4.clear();
         Contract(1.0,t[row][col],shape(0),LOu,shape(0),0.0,I4);

         // 1) construct left Sp operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOp.clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(0,1),0.0,LOp);

         // 2) construct left Sm operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOm.clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(0,1),0.0,LOm);

         // 3) construct left Sz operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOz.clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(0,1),0.0,LOz);

         // 4) finally construct new left unity
         Heisenberg::construct_double_layer('H',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOu.clear();
         Contract(1.0,I4bis,shape(2,3),b[row-1][col],shape(0,1),0.0,LOu);

      }

      //last site on the right: close down on the incomings

      //1) first Lp with Sm
      Heisenberg::construct_double_layer('H',peps(row,Lx-1),Sm,dlom);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][Lx - 1],shape(1),dlom,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(t[row][Lx - 1].shape(0),dlom.shape(0),b[row-1][Lx - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOp,RO[Lx - 3]);

      //2) then Lm with Sp
      Heisenberg::construct_double_layer('H',peps(row,Lx-1),Sp,dlop);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][Lx - 1],shape(1),dlop,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(t[row][Lx - 1].shape(0),dlop.shape(0),b[row-1][Lx - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOm,RO[Lx - 3]);

      //3) then Lz with Sz
      Heisenberg::construct_double_layer('H',peps(row,Lx-1),Sz,dloz);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,t[row][Lx - 1],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(t[row][Lx - 1].shape(0),dloz.shape(0),b[row-1][Lx - 1].shape(0)));

      //add to value
      val += Dot(LOz,RO[Lx - 3]);

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,t[Ly-2][Lx - 1],shape(1),b[Ly-2][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(t[Ly-2][Lx - 1].shape(0),b[Ly-2][Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 1;--col){

      I.clear();
      Contract(1.0,t[Ly-2][col],shape(2),R[col-1],shape(0),0.0,I);

      R[col-2].clear();
      Contract(1.0,I,shape(1,2),b[Ly-2][col],shape(1,2),0.0,R[col-2]);

   }

   //construct the left going operators on the first top site

   //first S+
   construct_double_layer('H',peps(Ly-1,0),Sp,dlsp);

   //tmp comes out index (t,b)
   Contract(1.0,dlsp,shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   Lp = tmp.reshape_clear(shape(dlsp.shape(2),b[Ly-2][0].shape(2)));

   //then S-
   construct_double_layer('H',peps(Ly-1,0),Sm,dlsm);

   //tmp comes out index (t,b)
   Contract(1.0,dlsm,shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   Lm = tmp.reshape_clear(shape(dlsm.shape(2),b[Ly-2][0].shape(2)));

   //then Sz 
   construct_double_layer('H',peps(Ly-1,0),Sz,dlsz);

   //tmp comes out index (t,b)
   Contract(1.0,dlsz,shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   Lz = tmp.reshape_clear(shape(dlsz.shape(2),b[Ly-2][0].shape(2)));

   //and finally unity
   Contract(1.0,t[Ly-2][0],shape(1),b[Ly-2][0],shape(1),0.0,tmp);

   Lu = tmp.reshape_clear(shape(t[Ly-2][0].shape(2),b[Ly-2][0].shape(2)));

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the +,- and z terms from the previous site

      //construct the right intermediate contraction (paste bottom to right)
      I.clear();
      Contract(1.0,b[Ly-2][col],shape(2),R[col - 1],shape(1),0.0,I);

      // 1) construct Sm double layer
      construct_double_layer('H',peps(Ly-1,col),Sm,dlsm);

      R[col-1].clear();
      Contract(1.0,dlsm,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

      //contract with left S+
      val += 0.5 * Dot(Lp,R[col - 1]);

      // 2) construct Sp double layer
      construct_double_layer('H',peps(Ly-1,col),Sp,dlsp);

      R[col-1].clear();
      Contract(1.0,dlsp,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

      //contract with left S-
      val += 0.5 * Dot(Lm,R[col - 1]);

      // 3) construct Sz double layer
      construct_double_layer('H',peps(Ly-1,col),Sz,dlsz);

      R[col-1].clear();
      Contract(1.0,dlsz,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

      //contract with left Sz
      val += Dot(Lz,R[col - 1]);

      //construct left renormalized operators for next site: first paste bottom to Left unity
      I.clear();
      Contract(1.0,Lu,shape(1),b[Ly-2][col],shape(0),0.0,I);

      // 1) construct new Sm left operator
      Lm.clear();
      Contract(1.0,dlsm,shape(0,1),I,shape(0,1),0.0,Lm);

      // 2) construct new Sp left operator
      Lp.clear();
      Contract(1.0,dlsp,shape(0,1),I,shape(0,1),0.0,Lp);

      // 3) construct new Sz left operator
      Lz.clear();
      Contract(1.0,dlsz,shape(0,1),I,shape(0,1),0.0,Lz);

      // 4) finally construct new unity on the left
      Lu.clear();
      Contract(1.0,t[Ly-2][col],shape(0,1),I,shape(0,1),0.0,Lu);

   }

   //finally close down on last top site

   //1) Sm to close down Lp
   construct_double_layer('H',peps(Ly-1,Lx-1),Sm,dlsm);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,dlsm,shape(1),b[Ly-2][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(dlsm.shape(0),b[Ly-2][Lx - 1].shape(0)));

   val += 0.5 * Dot(Lp,R[Lx-3]);

   //2) Sp to close down Lm
   construct_double_layer('H',peps(Ly-1,Lx-1),Sp,dlsp);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,dlsp,shape(1),b[Ly-2][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(dlsp.shape(0),b[Ly-2][Lx - 1].shape(0)));

   val += 0.5 * Dot(Lm,R[Lx-3]);

   //3) Sz to close down Lz
   construct_double_layer('H',peps(Ly-1,Lx-1),Sz,dlsz);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,dlsz,shape(1),b[Ly-2][Lx - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(dlsz.shape(0),b[Ly-2][Lx - 1].shape(0)));

   val += Dot(Lz,R[Lx-3]);

   // #################################################################
   // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || left column: similar to overlap calculation

   //first construct the right renormalized operators

   //first the rightmost operator

   //tmp comes out index (r,l)
   Contract(1.0,r[0][Ly - 1],shape(1),l[0][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(r[0][Ly - 1].shape(0),l[0][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 1;--row){

      I.clear();
      Contract(1.0,r[0][row],shape(2),R[row-1],shape(0),0.0,I);

      R[row-2].clear();
      Contract(1.0,I,shape(1,2),l[0][row],shape(1,2),0.0,R[row-2]);

   }

   //4 left going operators: S+, S-, Sz, and 1

   //first S+
   construct_double_layer('V',peps(0,0),Sp,dlsp);

   //tmp comes out index (r,l)
   Contract(1.0,r[0][0],shape(1),dlsp,shape(1),0.0,tmp);

   Lp = tmp.reshape_clear(shape(r[0][0].shape(2),dlsp.shape(2)));

   //then S-
   construct_double_layer('V',peps(0,0),Sm,dlsm);

   //tmp comes out index (r,l)
   Contract(1.0,r[0][0],shape(1),dlsm,shape(1),0.0,tmp);

   Lm = tmp.reshape_clear(shape(r[0][0].shape(2),dlsm.shape(2)));

   //then Sz 
   construct_double_layer('V',peps(0,0),Sz,dlsz);

   //tmp comes out index (r,l)
   Contract(1.0,r[0][0],shape(1),dlsz,shape(1),0.0,tmp);

   Lz = tmp.reshape_clear(shape(r[0][0].shape(2),dlsz.shape(2)));

   //and finally unity
   Contract(1.0,r[0][0],shape(1),l[0][0],shape(1),0.0,tmp);

   Lu = tmp.reshape_clear(shape(r[0][0].shape(2),l[0][0].shape(2)));

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //first close down the +,- and z terms from the previous site

      //construct the right intermediate contraction (paste 'right' to R)
      I.clear();
      Contract(1.0,r[0][row],shape(2),R[row - 1],shape(0),0.0,I);

      // 1) construct Sm double layer
      construct_double_layer('V',peps(row,0),Sm,dlsm);

      R[row-1].clear();
      Contract(1.0,I,shape(1,2),dlsm,shape(1,2),0.0,R[row - 1]);

      //contract with left S+
      val += 0.5 * Dot(Lp,R[row - 1]);

      // 2) then construct Sp double layer
      construct_double_layer('V',peps(row,0),Sp,dlsp);

      R[row-1].clear();
      Contract(1.0,I,shape(1,2),dlsp,shape(1,2),0.0,R[row - 1]);

      //contract with left S-
      val += 0.5 * Dot(Lm,R[row - 1]);

      // 3) then construct Sz double layer
      construct_double_layer('V',peps(row,0),Sz,dlsz);

      R[row-1].clear();
      Contract(1.0,I,shape(1,2),dlsz,shape(1,2),0.0,R[row - 1]);

      //contract with left Sz
      val += Dot(Lz,R[row - 1]);

      //construct left renormalized operators for next site: first paste top to Left unity
      I.clear();
      Contract(1.0,Lu,shape(0),r[0][row],shape(0),0.0,I);

      // 1) construct new Sm left operator
      Lm.clear();
      Contract(1.0,I,shape(0,1),dlsm,shape(0,1),0.0,Lm);

      // 2) construct new Sp left operator
      Lp.clear();
      Contract(1.0,I,shape(0,1),dlsp,shape(0,1),0.0,Lp);

      // 3) construct new Sz left operator
      Lz.clear();
      Contract(1.0,I,shape(0,1),dlsz,shape(0,1),0.0,Lz);

      // 4) finally construct new unity on the left
      Lu.clear();
      Contract(1.0,I,shape(0,1),l[0][row],shape(0,1),0.0,Lu);

   }

   //last site of left column: close down the left +,- and z

   //1) Sm to close down Lp
   construct_double_layer('V',peps(Ly-1,0),Sm,dlsm);

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,r[0][Ly - 1],shape(1),dlsm,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(r[0][Ly - 1].shape(0),dlsm.shape(0)));

   val += 0.5 * Dot(Lp,R[Ly-3]);

   //2) Sp to close down Lm
   construct_double_layer('V',peps(Ly-1,0),Sp,dlsp);

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,r[0][Ly - 1],shape(1),dlsp,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(r[0][Ly - 1].shape(0),dlsp.shape(0)));

   val += 0.5 * Dot(Lm,R[Ly-3]);

   //3) Sz to close down Lz
   construct_double_layer('V',peps(Ly-1,0),Sz,dlsz);

   //tmp comes out index (t,b)
   tmp.clear();
   Contract(1.0,r[0][Ly - 1],shape(1),dlsz,shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(r[0][Ly - 1].shape(0),dlsz.shape(0)));

   val += Dot(Lz,R[Ly-3]);

   // -- (2) -- now move from left to right calculating everything like an MPO/MPS expectation value
   for(int col = 1;col < Lx - 1;++col){

      //first create right renormalized operator

      //first site make double layer object from peps
      Heisenberg::construct_double_layer('V',peps(Ly-1,col),dlou);

      //paste right environment on
      DArray<5> tmp5;
      Contract(1.0,r[col][Ly - 1],shape(1),dlou,shape(1),0.0,tmp5);

      //then left enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),l[col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(r[col][Ly - 1].shape(0),dlou.shape(0),l[col-1][Ly - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int row = Ly-2;row > 1;--row){

         I4.clear();
         Contract(1.0,r[col][row],shape(2),RO[row-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         Heisenberg::construct_double_layer('V',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-2].clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(1,2),0.0,RO[row-2]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) S+ -- make double layer object from peps with Sp
      Heisenberg::construct_double_layer('V',peps(0,col),Sp,dlop);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][0],shape(1),dlop,shape(1),0.0,tmp5);

      //then left enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (right-env,peps-col,left-env)
      LOp = tmp6.reshape_clear(shape(r[col][0].shape(2),dlop.shape(3),l[col-1][0].shape(2)));

      // 2) S- -- make double layer object from peps with Sm
      Heisenberg::construct_double_layer('V',peps(0,col),Sm,dlom);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][0],shape(1),dlom,shape(1),0.0,tmp5);

      //then left enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: 
      LOm = tmp6.reshape_clear(shape(r[col][0].shape(2),dlom.shape(3),l[col-1][0].shape(2)));

      // 3) Sz -- make double layer object from peps with Sz
      Heisenberg::construct_double_layer('V',peps(0,col),Sz,dloz);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][0],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOz = tmp6.reshape_clear(shape(r[col][0].shape(2),dlom.shape(3),l[col-1][0].shape(2)));

      // 4) 1 -- finally construct left renormalized operator with unity
      Heisenberg::construct_double_layer('V',peps(0,col),dlou);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][0],shape(1),dlou,shape(1),0.0,tmp5);

      //then left enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: 
      LOu = tmp6.reshape_clear(shape(r[col][0].shape(2),dlou.shape(3),l[col-1][0].shape(2)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int row = 1;row < Ly - 1;++row){

         //first add right to the right side, put it in I4
         I4.clear();
         Contract(1.0,r[col][row],shape(2),RO[row-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         //1) close down LOp with Sm
         Heisenberg::construct_double_layer('V',peps(row,col),Sm,dlom);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-1].clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         val += 0.5 * Dot(LOp,RO[row-1]);

         //2) close down LOm with Sp
         Heisenberg::construct_double_layer('V',peps(row,col),Sp,dlop);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-1].clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         val += 0.5 * Dot(LOm,RO[row-1]);

         //3) finally close down LOz with Sz
         Heisenberg::construct_double_layer('V',peps(row,col),Sz,dloz);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-1].clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         val += Dot(LOz,RO[row-1]);

         // now construct the new left going renormalized operators
         //first attach top to left unity
         I4.clear();
         Contract(1.0,r[col][row],shape(0),LOu,shape(0),0.0,I4);

         // 1) construct left Sp operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOp.clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(0,1),0.0,LOp);

         // 2) construct left Sm operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOm.clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(0,1),0.0,LOm);

         // 3) construct left Sz operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOz.clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(0,1),0.0,LOz);

         // 4) finally construct new left unity
         Heisenberg::construct_double_layer('V',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOu.clear();
         Contract(1.0,I4bis,shape(2,3),l[col-1][row],shape(0,1),0.0,LOu);

      }

      //last site on the right: close down on the incomings

      //1) first Lp with Sm
      Heisenberg::construct_double_layer('V',peps(Ly-1,col),Sm,dlom);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][Ly - 1],shape(1),dlom,shape(1),0.0,tmp5);

      //then left enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 3] = tmp6.reshape_clear(shape(r[col][Ly - 1].shape(0),dlom.shape(0),l[col-1][Ly - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOp,RO[Ly - 3]);

      //2) then Lm with Sp
      Heisenberg::construct_double_layer('V',peps(Ly-1,col),Sp,dlop);

      //paste right environment on
      tmp5.clear();
      Contract(1.0,r[col][Ly - 1],shape(1),dlop,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 3] = tmp6.reshape_clear(shape(r[col][Ly - 1].shape(0),dlop.shape(0),l[col-1][Ly - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOm,RO[Ly - 3]);

      //3) then Lz with Sz
      Heisenberg::construct_double_layer('V',peps(Ly-1,col),Sz,dloz);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,r[col][Ly - 1],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),l[col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 3] = tmp6.reshape_clear(shape(r[col][Ly - 1].shape(0),dloz.shape(0),l[col-1][Ly - 1].shape(0)));

      //add to value
      val += Dot(LOz,RO[Ly - 3]);

   }

   // -- (3) -- || right column = Lx-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,r[Lx-2][Ly - 1],shape(1),l[Lx-2][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Lx - 3] = tmp.reshape_clear(shape(r[Lx-2][Ly - 1].shape(0),l[Lx-2][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 1;--row){

      I.clear();
      Contract(1.0,r[Lx-2][row],shape(2),R[row-1],shape(0),0.0,I);

      R[row-2].clear();
      Contract(1.0,I,shape(1,2),l[Lx-2][row],shape(1,2),0.0,R[row-2]);

   }

   //construct the left going operators on the first top site

   //first S+
   construct_double_layer('V',peps(0,Lx-1),Sp,dlsp);

   //tmp comes out index (r,l)
   Contract(1.0,dlsp,shape(1),l[Lx-2][0],shape(1),0.0,tmp);

   Lp = tmp.reshape_clear(shape(dlsp.shape(2),l[Lx-2][0].shape(2)));

   //then S-
   construct_double_layer('V',peps(0,Lx-1),Sm,dlsm);

   //tmp comes out index (r,l)
   Contract(1.0,dlsm,shape(1),l[Lx-2][0],shape(1),0.0,tmp);

   Lm = tmp.reshape_clear(shape(dlsm.shape(2),l[Lx-2][0].shape(2)));

   //then Sz 
   construct_double_layer('V',peps(0,Lx-1),Sz,dlsz);

   //tmp comes out index (r,l)
   Contract(1.0,dlsz,shape(1),l[Lx-2][0],shape(1),0.0,tmp);

   Lz = tmp.reshape_clear(shape(dlsz.shape(2),l[Lx-2][0].shape(2)));

   //and finally unity
   Contract(1.0,r[Lx-2][0],shape(1),l[Lx-2][0],shape(1),0.0,tmp);

   Lu = tmp.reshape_clear(shape(r[Lx-2][0].shape(2),l[Lx-2][0].shape(2)));

   //middle of the chain:
   for(int row = 1;row < Ly-1;++row){

      //first close down the +,- and z terms from the previous site

      //construct the right intermediate contraction (paste left to 'right')
      I.clear();
      Contract(1.0,l[Lx-2][row],shape(2),R[row - 1],shape(1),0.0,I);

      // 1) construct Sm double layer
      construct_double_layer('V',peps(row,Lx-1),Sm,dlsm);

      R[row-1].clear();
      Contract(1.0,dlsm,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

      //contract with left S+
      val += 0.5 * Dot(Lp,R[row - 1]);

      // 2) construct Sp double layer
      construct_double_layer('V',peps(row,Lx-1),Sp,dlsp);

      R[row-1].clear();
      Contract(1.0,dlsp,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

      //contract with left S-
      val += 0.5 * Dot(Lm,R[row - 1]);

      // 3) construct Sz double layer
      construct_double_layer('V',peps(row,Lx-1),Sz,dlsz);

      R[row-1].clear();
      Contract(1.0,dlsz,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

      //contract with left Sz
      val += Dot(Lz,R[row - 1]);

      //construct left renormalized operators for next site: first paste bottom to Left unity
      I.clear();
      Contract(1.0,Lu,shape(1),l[Lx-2][row],shape(0),0.0,I);

      // 1) construct new Sm left operator
      Lm.clear();
      Contract(1.0,dlsm,shape(0,1),I,shape(0,1),0.0,Lm);

      // 2) construct new Sp left operator
      Lp.clear();
      Contract(1.0,dlsp,shape(0,1),I,shape(0,1),0.0,Lp);

      // 3) construct new Sz left operator
      Lz.clear();
      Contract(1.0,dlsz,shape(0,1),I,shape(0,1),0.0,Lz);

      // 4) finally construct new unity on the left
      Lu.clear();
      Contract(1.0,r[Lx-2][row],shape(0,1),I,shape(0,1),0.0,Lu);

   }

   //finally close down on last 'right' site

   //1) Sm to close down Lp
   construct_double_layer('V',peps(Ly-1,Lx-1),Sm,dlsm);

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,dlsm,shape(1),l[Lx-2][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(dlsm.shape(0),l[Lx-2][Ly - 1].shape(0)));

   val += 0.5 * Dot(Lp,R[Ly-3]);

   //2) Sp to close down Lm
   construct_double_layer('V',peps(Ly-1,Lx-1),Sp,dlsp);

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,dlsp,shape(1),l[Lx-2][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(dlsp.shape(0),l[Lx-2][Ly - 1].shape(0)));

   val += 0.5 * Dot(Lm,R[Ly-3]);

   //3) Sz to close down Lz
   construct_double_layer('V',peps(Ly-1,Lx-1),Sz,dlsz);

   //tmp comes out index (r,l)
   tmp.clear();
   Contract(1.0,dlsz,shape(1),l[Lx-2][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(dlsz.shape(0),l[Lx-2][Ly - 1].shape(0)));

   val += Dot(Lz,R[Ly-3]);

   return val;

}

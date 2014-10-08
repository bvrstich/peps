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

using namespace global;

/** 
 * empty constructor
 */
Environment::Environment(){ }

/** 
 * constructor with allocation
 * @param D_in bond dimension of peps state
 * @param D_aux_in contraction bond dimension
 */
Environment::Environment(int D_in,int D_aux_in){

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

   D = D_in;
   D_aux = D_aux_in;

   //allocate the memory
   
   //bottom
   int tmp = D*D;

   for(int i = 0;i < Ly - 1;++i){

      if(tmp < D_aux)
         b[i] = MPO<double>(Lx,D,tmp);
      else
         b[i] = MPO<double>(Lx,D,D_aux);

      tmp *= D*D;

   }
   
   //top
   tmp = D*D;

   for(int i = Ly - 2;i >= 0;--i){

      if(tmp < D_aux)
         t[i] = MPO<double>(Lx,D,tmp);
      else
         t[i] = MPO<double>(Lx,D,D_aux);

      tmp *= D*D;

   }

   //left
   tmp = D*D;

   for(int i = 0;i < Lx - 1;++i){

      if(tmp < D_aux)
         l[i] = MPO<double>(Ly,D,tmp);
      else
         l[i] = MPO<double>(Ly,D,D_aux);

      tmp *= D*D;

   }
   
   //finally right
   tmp = D*D;

   for(int i = Lx - 2;i >= 0;--i){

      if(tmp < D_aux)
         r[i] = MPO<double>(Ly,D,tmp);
      else
         r[i] = MPO<double>(Ly,D,D_aux);

      tmp *= D*D;

   }

}

/** 
 * copy constructor with allocation
 */
Environment::Environment(const Environment &env_copy){

   t = env_copy.gt();
   b = env_copy.gb();

   r = env_copy.gr();
   l = env_copy.gl();

   D = env_copy.gD();
   D_aux = env_copy.gD_aux();

}

/**
 * empty destructor
 */
Environment::~Environment(){ }

/**
 * construct the enviroment mps's for the input PEPS
 * @param option if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 *               if 'A' construct all environments
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc(const char option,const PEPS<double> &peps){

   if(option == 'B' || option == 'A'){

      //construct bottom layer
      b[0].fill('b',peps);

      //for(int i = 1;i < Ly - 1;++i)
      int i = 1;
         this->add_layer('b',i,peps,5);

   }

   if(option == 'T' || option == 'A'){

      t[Ly - 2].fill('t',peps);

      for(int i = Ly - 3;i >= 0;--i)
         this->add_layer('t',i,peps,5);

   }
/*
   if(option == 'L' || option == 'A'){

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

   }

   if(option == 'R' || option == 'A'){

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

   }
   */
}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test(){

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
 * construct the enviroment for a specific row/column of the input PEPS, it is assumed that 
 * all prerequisites for the construction of the environment are there!
 * @param option if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 * @param rc row or column index for the L,R,T or B environment to be constructed
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc(char option,int rc,const PEPS<double> &peps,int D_aux){
   /* 
      if(option == 'B'){

   //construct bottom layer
   if(rc == 0)
   b[0] = MPS<double>('b',peps,peps);
   else{

   //i'th row as MPO
   MPO<double> mpo('H',rc,peps,peps);

   MPS<double> tmp(b[rc - 1]);

   //apply to form MPS with bond dimension D^4
   tmp.gemv('L',mpo);

   //reduce the dimensions of the edge states using thin svd
   tmp.cut_edges();

   //compress in sweeping fashion
   b[rc].resize(Lx);
   b[rc].compress(D_aux,tmp,5);

   }

   }
   else if(option == 'T'){

   //then construct top layer
   if(rc == Ly-1)
   t[Ly - 2] = MPS<double>('t',peps,peps);
   else{

   //i'th row as MPO
   MPO<double> mpo('H',rc,peps,peps);

   //apply to form MPS with bond dimension D^4
   MPS<double> tmp(t[rc]);

   tmp.gemv('U',mpo);

   //reduce the dimensions of the edge states using thin svd
   tmp.cut_edges();

   //compress in sweeping fashion
   t[rc - 1].resize(Lx);
   t[rc - 1].compress(D_aux,tmp,5);

   }

   }
   else if(option == 'L'){

   //then left layer
   if(rc == 0)
   l[0] = MPS<double>('l',peps,peps);
   else{

   //i'th col as MPO
   MPO<double> mpo('V',rc,peps,peps);

   MPS<double> tmp(l[rc - 1]);

   //apply to form MPS with bond dimension D^4
   tmp.gemv('L',mpo);

   //reduce the dimensions of the edge states using thin svd
   tmp.cut_edges();

   //compress in sweeping fashion
   l[rc].resize(Ly);
   l[rc].compress(D_aux,tmp,5);

}

}
else{//option == R

   //finally construct right layer
   if(rc == Lx - 1)
      r[Lx - 2] = MPS<double>('r',peps,peps);
   else{

      //i'th row as MPO
      MPO<double> mpo('V',rc,peps,peps);

      //apply to form MPS with bond dimension D^4
      MPS<double> tmp(r[rc]);

      tmp.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      r[rc - 1].resize(Ly);
      r[rc - 1].compress(D_aux,tmp,5);

   }

}
*/
}

/**
 * construct a double layer MPO object from a PEPS
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param dlo output object
 */
void Environment::construct_double_layer(char option,const DArray<5> &peps,DArray<4> &dlo){

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
void Environment::construct_double_layer(char option,const DArray<5> &peps,const DArray<2> &O,DArray<4> &dlo){

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
 * construct a double layer MPS object from a PEPS taking the expectation value of an operator O between the physical indices
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param O single particle operator
 * @param dls output object
 */
void Environment::construct_double_layer(char option,const DArray<5> &peps,const DArray<2> &O,DArray<3> &dls){

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
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(s,p,m,q,i,n,j,o));

      int DL = tmp2.shape(0) * tmp2.shape(1);
      int d_phys = tmp2.shape(2) * tmp2.shape(3) * tmp2.shape(4) * tmp2.shape(5);
      int DR = tmp2.shape(6) * tmp2.shape(7);

      dls = tmp2.reshape_clear(shape(DL,d_phys,DR));

   }

}

/**
 * construct a double layer MPS object from a PEPS 
 * @param option == 'H' horizontal if 'V' vertical
 * @param peps THe input PEPS elements
 * @param dls output object
 */
void Environment::construct_double_layer(char option,const DArray<5> &peps,DArray<3> &dls){

   if(option == 'H'){

      enum {i,j,k,m,n,o,p,q,s};

      DArray<8> tmp;
      Contract(1.0,peps,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,s,p,m,q));

      int DL = tmp.shape(0) * tmp.shape(1);
      int d_phys = tmp.shape(2) * tmp.shape(3) * tmp.shape(4) * tmp.shape(5);
      int DR = tmp.shape(6) * tmp.shape(7);

      dls = tmp.reshape_clear(shape(DL,d_phys,DR));

   }
   else{//V

      enum {i,j,k,s,m,n,o,p,q};

      DArray<8> tmp;
      Contract(1.0,peps,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(s,p,m,q,i,n,j,o));

      int DL = tmp.shape(0) * tmp.shape(1);
      int d_phys = tmp.shape(2) * tmp.shape(3) * tmp.shape(4) * tmp.shape(5);
      int DR = tmp.shape(6) * tmp.shape(7);

      dls = tmp.reshape_clear(shape(DL,d_phys,DR));

   }

}

/**
 * const version
 * @param col the column index
 * @return the right boundary 'MPO' environment on column col
 */
const MPO<double> &Environment::gr(int col) const {

   return r[col];

}

/**
 * @param col the column index: access version
 * @return the right boundary 'MPO' environment on column col
 */
MPO<double> &Environment::gr(int col) {

   return r[col];

}

/**
 * const version
 * @param col the column index
 * @return the left boundary 'MPO' environment on column col
 */
const MPO<double> &Environment::gl(int col) const {

   return l[col];

}

/**
 * @param col the column index: access version
 * @return the left boundary 'MPO' environment on column col
 */
MPO<double> &Environment::gl(int col) {

   return l[col];

}

/**
 * const version
 * @param row the row index
 * @return the top boundary 'MPO' environment on row 'row'
 */
const MPO<double> &Environment::gt(int row) const {

   return t[row];

}

/**
 * access version
 * @param row the row index
 * @return the top boundary 'MPO' environment on row 'row'
 */
MPO<double> &Environment::gt(int row) {

   return t[row];

}

/**
 * const version
 * @param row the row index
 * @return the bottom boundary 'MPO' environment on row 'row'
 */
const MPO<double> &Environment::gb(int row) const {

   return b[row];

}

/**
 * access version
 * @param row the row index
 * @return the bottom boundary 'MPO' environment on row 'row'
 */
MPO<double> &Environment::gb(int row) {

   return b[row];

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
const int Environment::gD_aux() const {

   return D_aux;

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
const int Environment::gD() const {

   return D;

}

/**
 * set a new bond dimension
 */
void Environment::sD(int D_in) {

   D = D_in;

}

/**
 * set a new auxiliary bond dimension
 */
void Environment::sD_aux(int D_aux_in) {

   D_aux = D_aux_in;

}

/**
 * @return the full bottom boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gb() const {

   return b;

}

/**
 * @return the full top boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gt() const {

   return t;

}

/**
 * @return the full left boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gl() const {

   return l;

}

/**
 * @return the full right boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gr() const {

   return r;

}

/**
 * construct the (t,b,l or r) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPO
 * @param option 't'op, 'b'ottom, 'l'eft or 'r'ight environment
 * @param rc row or column index
 * @param peps the input PEPS<double> object 
 * @param n_iter number of iterations of compression sweep
 */
void Environment::add_layer(const char option,int rc,const PEPS<double> &peps,int n_iter){

   if(option == 'b'){

      b[rc].fill_Random();

      vector< DArray<4> > R(Lx - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,b[rc - 1][Lx - 1],shape(1),peps(rc,Lx - 1),shape(3),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,5),peps(rc,Lx - 1),shape(3,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),b[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

      R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(3,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,6),peps(rc,i),shape(3,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),b[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }
/*
      int iter = 0;

      while(iter < n_iter){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,b[rc - 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(rc,0),shape(1,4),tmp6,shape(1,3),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(rc,0),shape(1,2,4),tmp7,shape(4,1,5),0.0,tmp6);

         b[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(b[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,b[rc-1][0],shape(1),peps(rc,0),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc,0),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(rc,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,b[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(b[rc][i],tmp2);

            //construct left operator

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),b[rc - 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(rc,i),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,4),peps(rc,i),shape(0,1,2),0.0,tmp6);

            Contract(1.0,tmp6,shape(0,2,4),b[rc][i],shape(0,1,2),0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),b[rc - 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(rc,Lx - 1),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,4),peps(rc,Lx - 1),shape(0,1,2),0.0,tmp6);

         b[rc][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,b[rc][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,b[rc - 1][Lx - 1],shape(1),peps(rc,Lx - 1),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc,Lx - 1),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(rc,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,b[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,b[rc][i]);

            //construct right operator
            tmp6.clear();
            Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps(rc,i),shape(1,2,4),0.0,tmp6);

            Contract(1.0,tmp6,shape(3,5,1),b[rc][i],shape(1,2,3),0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,b[rc][0],tmp2,0.0,tmp4);

         b[rc][0] = std::move(tmp4);

         ++iter;

      }
*/
   }
   else if(option == 't'){

      t[rc].fill_Random();

      vector< DArray<4> > R(Lx - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,t[rc + 1][Lx - 1],shape(1),peps(rc+1,Lx - 1),shape(1),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(rc+1,Lx - 1),shape(1,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),t[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

      R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps(rc+1,i),shape(1,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),t[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < n_iter){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,t[rc + 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(rc+1,0),shape(1,4),tmp6,shape(1,3),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(rc+1,0),shape(1,2,4),tmp7,shape(4,1,5),0.0,tmp6);

         t[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(t[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,t[rc+1][0],shape(1),peps(rc+1,0),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc+1,0),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(rc+1,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,t[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(t[rc][i],tmp2);

            //construct left operator
            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),t[rc + 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(rc+1,i),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,4),peps(rc+1,i),shape(0,1,2),0.0,tmp6);

            Contract(1.0,tmp6,shape(0,2,4),t[rc][i],shape(0,1,2),0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),t[rc + 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(rc+1,Lx - 1),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,4),peps(rc+1,Lx - 1),shape(0,1,2),0.0,tmp6);

         t[rc][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,t[rc][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,t[rc + 1][Lx - 1],shape(1),peps(rc+1,Lx - 1),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc+1,Lx - 1),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(rc+1,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,t[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,t[rc][i]);

            //construct right operator
            tmp6.clear();
            Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps(rc+1,i),shape(1,4,2),0.0,tmp6);

            Contract(1.0,tmp6,shape(3,5,1),t[rc][i],shape(1,2,3),0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc][0],tmp2,0.0,tmp4);

         t[rc][0] = std::move(tmp4);

         ++iter;

      }


   }
   else if(option == 'l'){

   }
   else{

   }

}

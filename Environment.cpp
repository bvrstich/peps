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

//statics
vector< MPS<double> > Environment::l;
vector< MPS<double> > Environment::r;
vector< MPS<double> > Environment::t;
vector< MPS<double> > Environment::b;

/** 
 * constructor
 */
void Environment::init(){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

}

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
void Environment::calc_env(char option,const PEPS<double> &peps,int D_aux){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   int d = Global::lat.gd();

   if(option == 'B' || option == 'A'){

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

   }

   if(option == 'T' || option == 'A'){

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

   }

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

}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test_env(){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

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
void Environment::calc_env(char option,int rc,const PEPS<double> &peps,int D_aux){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   int d = Global::lat.gd();

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
/*
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
*/
   }
   else{//option == R
/*
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
*/
   }

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
      Contract(1.0,tmp,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp2,shape(s,p,i,n,m,q,j,o));

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
      Contract(1.0,peps,shape(i,j,k,s,m),peps,shape(n,o,k,p,q),0.0,tmp,shape(s,p,i,n,m,q,j,o));

      int DL = tmp.shape(0) * tmp.shape(1);
      int d_phys = tmp.shape(2) * tmp.shape(3) * tmp.shape(4) * tmp.shape(5);
      int DR = tmp.shape(6) * tmp.shape(7);

      dls = tmp.reshape_clear(shape(DL,d_phys,DR));

   }

}

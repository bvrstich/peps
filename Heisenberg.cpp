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
   cout << "top to bottom" << endl;
   cout << endl;

   //test
   for(int i = 0;i < Ly - 1;++i)
      cout << t[i].dot(b[i]) << endl;
   cout << endl;
   cout << "left to right" << endl;
   cout << endl;

   for(int i = 0;i < Lx - 1;++i)
      cout << r[i].dot(l[i]) << endl;

}

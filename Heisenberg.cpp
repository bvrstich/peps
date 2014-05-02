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

   //make it left canonicalized
   b[0].canonicalize(Left,false);

   for(int i = 1;i < Ly - 1;++i){

      //i'th row as MPO
      MPO<double> mpo(i,peps,peps);

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

   //make it left canonicalized
   t[Ly - 2].canonicalize(Left,false);

   for(int i = Ly - 2;i > 0;--i){

      //i'th row as MPO
      MPO<double> mpo(i,peps,peps);

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
/*
   //make it left canonicalized
   t[].canonicalize(Left,false);

   for(int i = Ly - 2;i > 0;--i){

      //i'th row as MPO
      MPO<double> mpo(i,peps,peps);

      //apply to form MPS with bond dimension D^4
      MPS<double> tmp(t[i]);

      tmp.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      tmp.cut_edges();

      //compress in sweeping fashion
      t[i - 1].resize(Lx);
      t[i - 1].compress(D_aux,tmp,5);

   }
*/

   //test
   for(int i = 0;i < Ly - 1;++i)
      cout << t[i].dot(b[i]) << endl;

}

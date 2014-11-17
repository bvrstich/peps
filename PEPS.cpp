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
 * construct an empty PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 */
template<typename T>
PEPS<T>::PEPS() : vector< TArray<T,5> >(Lx * Ly) { }

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(Lx * Ly) {

   D = D_in;
   
   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(1,D,d,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx-1 ].resize(D,D,d,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(1,1,d,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx-1 ].resize(D,1,d,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(D,D,d,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(D,1,d,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(1,D,d,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(D,D,d,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(D,D,d,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].generate(rgen<T>);

         Normalize((*this)[ r*Lx + c ]);
         Scal((T)D,(*this)[ r*Lx + c ]);

      }

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,5> >(peps_copy) {

   D = peps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
PEPS<T>::~PEPS(){ }

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,5> &PEPS<T>::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[r*Lx + c];

}

/**
 * @return the cutoff virtual dimension
 */
template<typename T>
int PEPS<T>::gD() const {

   return D;

}

/**
 * @param D_in value to the D to
 */
template<typename T>
void PEPS<T>::sD(int D_in) {

   this->D = D_in;

}

/**
 * initialize the peps to the direct sum of two antiferromagnetic D=1 structures
 * @param f jastrow factor
 */
template<>
void PEPS<double>::initialize_jastrow(double f) {

   enum {i,j,k,l,m,n,o,p,q,r,s};

   D = 2;

   //bottom row, first site
   (*this)[0].resize(1,D,d,1,D);
   (*this)[0] = 0.0;

   (*this)[0](0,0,0,0,0) = 1.0;
   (*this)[0](0,1,1,0,1) = 1.0;

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[col].resize(D,D,d,1,D);
      (*this)[col] = 0.0;

      (*this)[col](0,0,0,0,0) = f;
      (*this)[col](1,0,0,0,0) = 1.0;
      (*this)[col](0,1,1,0,1) = 1.0;
      (*this)[col](1,1,1,0,1) = f;

   }

   //bottom row, last site
   (*this)[Lx-1].resize(D,D,d,1,1);
   (*this)[Lx-1] = 0.0;

   (*this)[Lx-1](0,0,0,0,0) = f;
   (*this)[Lx-1](1,0,0,0,0) = 1.0;

   (*this)[Lx-1](0,1,1,0,0) = 1.0;
   (*this)[Lx-1](1,1,1,0,0) = f;

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      (*this)[row*Lx].resize(1,D,d,D,D);
      (*this)[row*Lx] = 0.0;

      (*this)[row*Lx](0,0,0,0,0) = f;
      (*this)[row*Lx](0,0,0,1,0) = 1.0;
      (*this)[row*Lx](0,1,1,0,1) = 1.0;
      (*this)[row*Lx](0,1,1,1,1) = f;

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         (*this)[row*Lx + col].resize(D,D,d,D,D);
         (*this)[row*Lx + col] = 0.0;

         (*this)[row*Lx + col](0,0,0,0,0) = f*f;
         (*this)[row*Lx + col](0,0,0,1,0) = f;
         (*this)[row*Lx + col](1,0,0,0,0) = f;
         (*this)[row*Lx + col](1,0,0,1,0) = 1.0;

         (*this)[row*Lx + col](0,1,1,0,1) = 1.0;
         (*this)[row*Lx + col](0,1,1,1,1) = f;
         (*this)[row*Lx + col](1,1,1,0,1) = f;
         (*this)[row*Lx + col](1,1,1,1,1) = f*f;

      }

      //rightmost site on row 'row'
      (*this)[row*Lx + Lx - 1].resize(D,D,d,D,1);
      (*this)[row*Lx + Lx - 1] = 0.0;

      (*this)[row*Lx + Lx - 1](0,0,0,0,0) = f*f;
      (*this)[row*Lx + Lx - 1](0,0,0,1,0) = f;
      (*this)[row*Lx + Lx - 1](1,0,0,0,0) = f;
      (*this)[row*Lx + Lx - 1](1,0,0,1,0) = 1.0;

      (*this)[row*Lx + Lx - 1](0,1,1,0,0) = 1.0;
      (*this)[row*Lx + Lx - 1](1,1,1,0,0) = f;
      (*this)[row*Lx + Lx - 1](0,1,1,1,0) = f;
      (*this)[row*Lx + Lx - 1](1,1,1,1,0) = f*f;

   }

   //top row
   //leftmost site
   (*this)[(Ly - 1)*Lx].resize(1,1,d,D,D);
   (*this)[(Ly - 1)*Lx] = 0.0;

   (*this)[(Ly - 1)*Lx](0,0,0,0,0) = f;
   (*this)[(Ly - 1)*Lx](0,0,0,1,0) = 1.0;
   (*this)[(Ly - 1)*Lx](0,0,1,0,1) = 1.0;
   (*this)[(Ly - 1)*Lx](0,0,1,1,1) = f;

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[(Ly - 1)*Lx + col].resize(D,1,d,D,D);
      (*this)[(Ly - 1)*Lx + col] = 0.0;

      (*this)[(Ly - 1)*Lx + col](0,0,0,0,0) = f*f;
      (*this)[(Ly - 1)*Lx + col](0,0,0,1,0) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,0,0,0) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,0,1,0) = 1.0;

      (*this)[(Ly - 1)*Lx + col](0,0,1,0,1) = 1.0;
      (*this)[(Ly - 1)*Lx + col](0,0,1,1,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,1,0,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,1,1,1) = f*f;

   }

   //top row rightmost site
   (*this)[(Ly - 1)*Lx + Lx - 1].resize(D,1,d,D,1);
   (*this)[(Ly - 1)*Lx + Lx - 1] = 0.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,0,0) = f*f;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,1,0) = 1.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,1,0,0) = 1.0;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,1,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,1,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,1,1,0) = f*f;

}

/**
 * increase the bond dimension by one
 * @param D_in bond dimension to grow to
 * @param noise level of noise added to the initial state
 */
template<>
void PEPS<double>::grow_bond_dimension(int D_in,double noise) {

   D = D_in;

   DArray<5> tmp;

   //bottom row, first site
   tmp.resize(1,D,d,1,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[0].shape(0);++i)
      for(int j = 0;j < (*this)[0].shape(1);++j)
         for(int k = 0;k < (*this)[0].shape(2);++k)
            for(int l = 0;l < (*this)[0].shape(3);++l)
               for(int m = 0;m < (*this)[0].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[0](i,j,k,l,m);

   (*this)[0] = std::move(tmp);

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,D,d,1,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[col].shape(0);++i)
         for(int j = 0;j < (*this)[col].shape(1);++j)
            for(int k = 0;k < (*this)[col].shape(2);++k)
               for(int l = 0;l < (*this)[col].shape(3);++l)
                  for(int m = 0;m < (*this)[col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[col](i,j,k,l,m);

      (*this)[col] = std::move(tmp);

   }

   //bottom row, last site
   tmp.clear();
   tmp.resize(D,D,d,1,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[Lx-1].shape(0);++i)
      for(int j = 0;j < (*this)[Lx-1].shape(1);++j)
         for(int k = 0;k < (*this)[Lx-1].shape(2);++k)
            for(int l = 0;l < (*this)[Lx-1].shape(3);++l)
               for(int m = 0;m < (*this)[Lx-1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[Lx - 1](i,j,k,l,m);

   (*this)[Lx-1] = std::move(tmp);

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      tmp.clear();
      tmp.resize(1,D,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx](i,j,k,l,m);

      (*this)[row*Lx] = std::move(tmp);

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         tmp.clear();
         tmp.resize(D,D,d,D,D);

         tmp.generate(rgen<double>);
         Scal(noise,tmp);

         for(int i = 0;i < (*this)[row*Lx + col].shape(0);++i)
            for(int j = 0;j < (*this)[row*Lx + col].shape(1);++j)
               for(int k = 0;k < (*this)[row*Lx + col].shape(2);++k)
                  for(int l = 0;l < (*this)[row*Lx + col].shape(3);++l)
                     for(int m = 0;m < (*this)[row*Lx + col].shape(4);++m)
                        tmp(i,j,k,l,m) += (*this)[row*Lx + col](i,j,k,l,m);

         (*this)[row*Lx + col] = std::move(tmp);

      }

      //rightmost site on row 'row'
      tmp.clear();
      tmp.resize(D,D,d,D,1);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx + Lx - 1].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx + Lx - 1].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx + Lx - 1].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx + Lx - 1].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx + Lx - 1].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx + Lx - 1](i,j,k,l,m);

      (*this)[row*Lx + Lx - 1] = std::move(tmp);

   }

   //top row
   //leftmost site
   tmp.clear();
   tmp.resize(1,1,d,D,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx] = std::move(tmp);

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,1,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[(Ly - 1)*Lx + col].shape(0);++i)
         for(int j = 0;j < (*this)[(Ly - 1)*Lx + col].shape(1);++j)
            for(int k = 0;k < (*this)[(Ly - 1)*Lx + col].shape(2);++k)
               for(int l = 0;l < (*this)[(Ly - 1)*Lx + col].shape(3);++l)
                  for(int m = 0;m < (*this)[(Ly - 1)*Lx + col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + col](i,j,k,l,m);

      (*this)[(Ly - 1)*Lx + col] = std::move(tmp);

   }

   //top row rightmost site
   tmp.clear();
   tmp.resize(D,1,d,D,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx + Lx - 1].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx + Lx - 1].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx + Lx - 1].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx + Lx - 1].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx + Lx - 1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + Lx - 1](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx + Lx - 1] = std::move(tmp);

}

/**
 * @param peps_i peps to take the overlap with
 * @param init boolean, if true the top-bottom environment has already been calculated and we just need to do the MPO-MPO contraction
 * @return the inner product of two PEPS <psi1|psi2> 
 */
template<>
double PEPS<double>::dot(const PEPS<double> &peps_i,bool init) const {

   int half = Ly/2;

   if(!init){

      //construct bottom environment until half
      env.gb(0).fill('b',peps_i);

      for(int i = 1;i <= half;++i)
         env.add_layer('b',i,peps_i);

      env.gt(Ly - 2).fill('t',peps_i);

      for(int i = Ly - 3;i >= half;--i)
         env.add_layer('t',i,peps_i);

   }

   return env.gb(half).dot(env.gt(half));

}

/** 
 * normalize the peps approximately, using a contraction with auxiliary dimension
 * @param init if true the environment has already been initialized and does not need to be calculated
 */
template<>
void PEPS<double>::normalize(bool init){

   double val = sqrt(this->dot(*this,init));
   val = pow(val,1.0/(double)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(1.0/val,(*this)[ r*Lx + c ]);

}

/**
 * scale the peps with a number
 * @param val scalar to be multiplied with the peps
 */
template<typename T>
void PEPS<T>::scal(T val){

   val = pow(val,(T)1.0/(T)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(val,(*this)[ r*Lx + c ]);

}

/**
 * @param mpx will be written to file
 * @param filename name of the file
 * save the MPX object to a file in binary format.
 */

template<typename T>
void PEPS<T>::save(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ofstream fout(name);
         fout.precision(16);

         int Da = (*this)(row,col).shape(0);
         int Db = (*this)(row,col).shape(1);
         int Dc = (*this)(row,col).shape(2);
         int Dd = (*this)(row,col).shape(3);
         int De = (*this)(row,col).shape(4);

         fout << Da << "\t" << Db << "\t" << Dc << "\t" << Dd << "\t" << De << endl;

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fout << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\t" << e_ << "\t" << (*this)(row,col)(a_,b_,c_,d_,e_) << endl;

      }

}

/**
 * @param mpx will be constructed from file
 * @param filename name of the file
 * load the MPX object from a file in binary format.
 */
template<typename T>
void PEPS<T>::load(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ifstream fin(name);

         int Da,Db,Dc,Dd,De;

         fin >> Da >> Db >> Dc >> Dd >> De;

         (*this)(row,col).resize(Da,Db,Dc,Dd,De);

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fin >> a_ >> b_ >> c_ >> d_ >> e_ >> (*this)(row,col)(a_,b_,c_,d_,e_);

      }

}

/**
 * evaluate the expectation value of the energy for the nn-Heisenberg model
 * beware, the environments have to be constructed beforehand!
 */
template<>
double PEPS<double>::energy(){

   // ---- || evaluate the energy in an MPO/MPS manner, first from bottom to top, then left to right || ----
   int delta = ham.gdelta();

   // #################################################################
   // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || bottom row: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<3> > R(Lx - 1);

   contractions::init_ro('b',*this,R); 

   //left going operators: Li
   std::vector< DArray<3> > Li( delta ); 

   //left unit operator
   DArray<3> Lu;

   //peps contracted with a local operator
   std::vector< DArray<5> > peps_i( delta );

   DArray<4> tmp4;
   DArray<5> tmp5;
   DArray<6> tmp6;
   DArray<6> tmp6bis;
   DArray<7> tmp7;
   DArray<8> tmp8;
   DArray<8> tmp8bis;
   DArray<10> tmp10;

   //construct the left operator with two open physical bonds
   Contract(1.0,env.gt(0)[0],shape(1),(*this)(0,0),shape(1),0.0,tmp7);

   Contract(1.0,tmp7,shape(1),(*this)(0,0),shape(1),0.0,tmp10);

   DArray<10> tmp10bis;
   Permute(tmp10,shape(0,2,6,4,8,1,5,9,3,7),tmp10bis);

   int m = tmp10.shape(1) * tmp10.shape(5) * tmp10.shape(9);
   int n = tmp10.shape(3) * tmp10.shape(7);

   //add left hamiltonian operators
   for(int i = 0;i < delta;++i){

      Li[i].resize(shape(tmp10.shape(1),tmp10.shape(5),tmp10.shape(9)));
      blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, ham.gL(i).data(), 1, 0.0, Li[i].data(), 1);

   }

   //and finally unity
   Lu.resize(shape(tmp10.shape(1),tmp10.shape(5),tmp10.shape(9)));
   blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, I.data(), 1, 0.0, Lu.data(), 1);

   double val = 0.0;

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //close down the left renormalized operators from the previous site

      //construct the intermediate contraction (paste top right)
      tmp5.clear();
      Contract(1.0,env.gt(0)[col],shape(3),R[col-1],shape(0),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,tmp5,shape(1,3),(*this)(0,col),shape(1,4),0.0,tmp6);

      for(int i = 0;i < delta;++i){

         //first contract with the peps with right operators
         Contract(1.0,ham.gR(i),shape(1),(*this)(0,col),shape(2),0.0,peps_i[i]);

         // paste the operator peps' to intermediary
         tmp5.clear();
         Contract(1.0,tmp6,shape(1,2,4),peps_i[i],shape(2,4,0),0.0,tmp5);

         //contract with left 
         val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp5.data(),1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,Lu,shape(0),env.gt(0)[col],shape(0),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,tmp5,shape(0,2),(*this)(0,col),shape(0,1),0.0,tmp6);

      // construct new left operators
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(0,col),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,tmp6,shape(3,0,1),peps_i[i],shape(0,1,2),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape(env.gt(0)[col].shape(3),(*this)(0,col).shape(4),(*this)(0,col).shape(4)));

      }

      // 4) finally construct new unity on the left
      Contract(1.0,tmp6,shape(0,1,3),(*this)(0,col),shape(0,1,2),0.0,tmp5);
      Lu = tmp5.reshape_clear(shape(env.gt(0)[col].shape(3),(*this)(0,col).shape(4),(*this)(0,col).shape(4)));

   }

   //last site of bottom row: close down the left +,- and z

   //contract with left operators for energy contribution
   for(int i = 0;i < delta;++i){

      //first contract with the peps with the right hamiltonian operators
      peps_i[i].clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(0,Lx-1),shape(2),0.0,peps_i[i]);

      //add top to left operator
      tmp5.clear();
      Contract(1.0,Li[i],shape(0),env.gt(0)[Lx - 1],shape(0),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,(*this)(0,Lx-1),shape(0,1),tmp5,shape(0,2),0.0,tmp6);

      val += ham.gcoef(i) * blas::dot(tmp6.size(), tmp6.data(), 1, peps_i[i].data(), 1);

   }

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< DArray<4> > RO(Lx - 2);

   //array of delta left renormalized operators needed
   std::vector< DArray<4> > LOi(delta);

   //left unit renormalized operator
   DArray<4> LOu;

   for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operator
      contractions::init_ro('H',row,*this,RO);

      // --- now move from left to right to get the expecation value of the interactions ---

      // --- First construct the left going operators for the first site -----

      //paste top environment on
      tmp7.clear();
      Contract(1.0,env.gt(row)[0],shape(1),(*this)(row,0),shape(1),0.0,tmp7);

      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gL(i),shape(1),(*this)(row,0),shape(2),0.0,peps_i[i]);

         // add peps_i's to intermediate
         tmp8.clear();
         Contract(1.0,tmp7,shape(4,1),peps_i[i],shape(0,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),env.gb(row-1)[0],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
         LOi[i] = tmp8bis.reshape_clear(shape(env.gt(row)[0].shape(3),(*this)(row,0).shape(4),(*this)(row,0).shape(4),env.gb(row-1)[0].shape(3)));

      }

      // 4) 1 -- finally construct left renormalized operator with unity
      Contract(1.0,tmp7,shape(1,4),(*this)(row,0),shape(1,2),0.0,tmp8);
      Contract(1.0,tmp8,shape(3,6),env.gb(row-1)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
      LOu = tmp8bis.reshape_clear(shape(env.gt(row)[0].shape(3),(*this)(row,0).shape(4),(*this)(row,0).shape(4),env.gb(row-1)[0].shape(3)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int col = 1;col < Lx - 1;++col){

         //add top and one peps to the right side
         tmp6.clear();
         Contract(1.0,env.gt(row)[col],shape(3),RO[col-1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),(*this)(row,col),shape(1,4),0.0,tmp7);

         for(int i = 0;i < delta;++i){

            //contract the peps with right operators: peps_i's
            peps_i[i].clear();
            Contract(1.0,ham.gR(i),shape(1),(*this)(row,col),shape(2),0.0,peps_i[i]);

            //close down the LOi's with peps_i's
            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps_i[i],shape(2,4,0),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[col - 1].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,RO[col - 1]);

            //expectation value:
            val += ham.gcoef(i) * Dot(LOi[i],RO[col-1]);

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity
         tmp6.clear();
         Contract(1.0,env.gt(row)[col],shape(0),LOu,shape(0),0.0,tmp6);

         //add peps to it, intermediary
         tmp7.clear();
         Contract(1.0,tmp6,shape(3,0),(*this)(row,col),shape(0,1),0.0,tmp7);

         //paste left operators to intermediary to make new left operators
         for(int i = 0;i < delta;++i){

            Contract(1.0,ham.gL(i),shape(1),(*this)(row,col),shape(2),0.0,peps_i[i]);

            tmp6.clear();
            Contract(1.0,tmp7,shape(4,2,0),peps_i[i],shape(0,1,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            LOi[i].clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LOi[i]);

         }

         //finally construct new left unity
         Contract(1.0,tmp7,shape(2,0,4),(*this)(row,col),shape(0,1,2),0.0,tmp6);
         Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

         LOu.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LOu);

      }

      //last site on the right: close down on the incomings
      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,Lx - 1),shape(2),0.0,peps_i[i]);

         tmp6.clear();
         Contract(1.0,env.gt(row)[Lx - 1],shape(0),LOi[i],shape(0),0.0,tmp6);

         //add peps to it
         tmp7.clear();
         Contract(1.0,tmp6,shape(3,0),(*this)(row,Lx - 1),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(4,2,0),peps_i[i],shape(0,1,2),0.0,tmp6);

         val += ham.gcoef(i) * blas::dot(tmp6.size(), tmp6.data(), 1, env.gb(row-1)[Lx-1].data(), 1);

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro('t',*this,R);

   //construct the left operator with two open physical bonds
   tmp7.clear();
   Contract(1.0,env.gb(Ly-2)[0],shape(1),(*this)(Ly-1,0),shape(3),0.0,tmp7);

   tmp10.clear();
   Contract(1.0,tmp7,shape(1),(*this)(Ly-1,0),shape(3),0.0,tmp10);

   tmp10bis.clear();
   Permute(tmp10,shape(0,6,2,3,7,5,9,1,8,4),tmp10bis);

   //multiply intermediate with Left hamiltonian operators
   for(int i = 0;i < delta;++i){

      Li[i].resize( shape( (*this)(Ly-1,0).shape(4),(*this)(Ly-1,0).shape(4),env.gb(Ly-2)[0].shape(3) ) );

      m = tmp10.shape(1) * tmp10.shape(5) * tmp10.shape(9);
      n = tmp10.shape(4) * tmp10.shape(8);

      blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, ham.gL(i).data(), 1, 0.0, Li[i].data(), 1);

   }

   //and finally unity
   Lu.resize( shape( (*this)(Ly-1,0).shape(4),(*this)(Ly-1,0).shape(4),env.gb(Ly-2)[0].shape(3) ) );
   blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, I.data(), 1, 0.0, Lu.data(), 1);

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //close down the left renormalized terms from the previous site

      //construct the intermediate contraction (paste top right)
      tmp5.clear();
      Contract(1.0,env.gb(Ly-2)[col],shape(3),R[col-1],shape(2),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(3,4),tmp5,shape(2,4),0.0,tmp6);

      //paste right hamiltonian operators to intermediary an contract with left renormalized operator
      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(Ly - 1,col),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,peps_i[i],shape(0,3,4),tmp6,shape(2,4,5),0.0,tmp5);

         //contract with left hamiltonian operator
         val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp5.data(), 1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,env.gb(Ly-2)[col],shape(0),Lu,shape(2),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(3,0),tmp5,shape(1,4),0.0,tmp6);

      //now paste on left hamiltonian operators 
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(Ly - 1,col),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,peps_i[i],shape(0,3,1),tmp6,shape(1,3,5),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape((*this)(Ly-1,col).shape(4),(*this)(Ly-1,col).shape(4),env.gb(Ly-2)[col].shape(3)));

      }
      //finally construct new unity on the left
      tmp5.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(2,3,0),tmp6,shape(1,3,5),0.0,tmp5);

      Lu = tmp5.reshape_clear(shape((*this)(Ly-1,col).shape(4),(*this)(Ly-1,col).shape(4),env.gb(Ly-2)[col].shape(3)));

   }

   //finally close down on last top site

   //construct intermediate
   tmp7.clear();
   Contract(1.0,(*this)(Ly-1,Lx-1),shape(3),env.gb(Ly-2)[Lx - 1],shape(2),0.0,tmp7);

   //) close down Li[i] right hamiltonian operators
   for(int i = 0;i < delta;++i){

      peps_i[i].clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly - 1,Lx - 1),shape(2),0.0,peps_i[i]);

      tmp8.clear();
      Contract(1.0,peps_i[i],shape(0,3),tmp7,shape(2,5),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp8.data(),1);

   }

   // #################################################################
   // ###   ---- from right left : contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || right column: similar to overlap calculation

   //first construct the right renormalized operators
   R.resize(Ly - 2);

   contractions::init_ro('r',*this,R);

   //construct the left operator with two open physical bonds
   tmp7.clear();
   Contract(1.0,env.gl(Ly - 2)[0],shape(1),(*this)(0,Lx-1),shape(0),0.0,tmp7);

   tmp10.clear();
   Contract(1.0,tmp7,shape(1),(*this)(0,Lx-1),shape(0),0.0,tmp10);

   tmp10bis.clear();
   Permute(tmp10,shape(0,4,8,5,9,1,2,6,3,7),tmp10bis);

   //paste the left hamiltonian operators to this intermediate
   for(int i = 0;i < delta;++i){

      Li[i].resize(shape(tmp10.shape(1),tmp10.shape(2),tmp10.shape(6)));

      int m = tmp10.shape(1) * tmp10.shape(2) * tmp10.shape(6);
      int n = tmp10.shape(3) * tmp10.shape(7);

      blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, ham.gL(i).data(), 1, 0.0, Li[i].data(), 1);

   }

   //construct left unit operator
   Lu.resize(shape(tmp10.shape(1),tmp10.shape(2),tmp10.shape(6)));
   blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, I.data(), 1, 0.0, Lu.data(), 1);

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //first close down the interaction terms from the previous site

      //construct the intermediate contraction (paste top right)
      tmp5.clear();
      Contract(1.0,env.gl(Lx - 2)[row],shape(3),R[row-1],shape(0),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,tmp5,shape(1,3),(*this)(row,Lx-1),shape(0,1),0.0,tmp6);

      //paste right hamiltonian operator to intermediary and constact with left renormalized operator
      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,Lx-1),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,tmp6,shape(3,1,2),peps_i[i],shape(0,1,2),0.0,tmp5);

         //contract with left S+
         val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp5.data(),1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,Lu,shape(0),env.gl(Lx - 2)[row],shape(0),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,tmp5,shape(0,2),(*this)(row,Lx - 1),shape(3,0),0.0,tmp6);

      // then contract with left Hamiltonian operator
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(row,Lx-1),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,tmp6,shape(0,1,4),peps_i[i],shape(3,1,0),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape(env.gl(Lx - 2)[row].shape(3),(*this)(row,Lx - 1).shape(1),(*this)(row,Lx - 1).shape(1)));

      }

      // finally construct new unity on the left
      Contract(1.0,tmp6,shape(0,1,4),(*this)(row,Lx-1),shape(3,0,2),0.0,tmp5);
      Lu = tmp5.reshape_clear(shape(env.gl(Lx - 2)[row].shape(3),(*this)(row,Lx - 1).shape(1),(*this)(row,Lx - 1).shape(1)));

   }

   //last site: close down the left operators with right hamiltonian operators

   //construct intermediate first
   tmp7.clear();
   Contract(1.0,env.gl(Ly-2)[Lx - 1],shape(1),(*this)(Ly-1,Lx-1),shape(0),0.0,tmp7);

   //contract with right operator for energy
   for(int i = 0;i < delta;++i){

      peps_i[i].clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,Lx-1),shape(2),0.0,peps_i[i]);

      tmp8.clear();
      Contract(1.0,tmp7,shape(4,1),peps_i[i],shape(0,1),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp8.data(), 1);

   }

   // -- (2) -- now move from right to left calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   RO.resize(Ly - 2);

   //loop over the columns
   for(int col = Lx - 2;col > 0;--col){

      //first create right renormalized operator
      contractions::init_ro('V',col,*this,RO);

      // --- now move from left to right to get the expecation value of the interactions ---

      // --- First construct the left going operators for the first site -----

      //paste left environment on
      tmp7.clear();
      Contract(1.0,env.gl(col - 1)[0],shape(1),(*this)(0,col),shape(0),0.0,tmp7);

      // add left hamiltonian operators to intermediate
      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gL(i),shape(1),(*this)(0,col),shape(2),0.0,peps_i[i]);

         tmp8.clear();
         Contract(1.0,tmp7,shape(4,1),peps_i[i],shape(0,1),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(4,7),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
         LOi[i] = tmp8bis.reshape_clear(shape(env.gl(col - 1)[0].shape(3),(*this)(0,col).shape(1),(*this)(0,col).shape(1),env.gr(col)[0].shape(3)));

      }

      // construct left renormalized operator with unity
      Contract(1.0,tmp7,shape(1,4),(*this)(0,col),shape(0,2),0.0,tmp8);
      Contract(1.0,tmp8,shape(4,7),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
      LOu = tmp8bis.reshape_clear(shape(env.gl(col - 1)[0].shape(3),(*this)(0,col).shape(1),(*this)(0,col).shape(1),env.gr(col)[0].shape(3)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int row = 1;row < Ly - 1;++row){

         //first add 'left' and one peps to the right side to construct intermediate
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[row],shape(3),RO[row-1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),(*this)(row,col),shape(0,1),0.0,tmp7);

         // close down LOp with right hamiltonian operators pasted on intermediate for energy contribution
         for(int i = 0;i < delta;++i){

            peps_i[i].clear();
            Contract(1.0,ham.gR(i),shape(1),(*this)(row,col),shape(2),0.0,peps_i[i]);

            tmp6.clear();
            Contract(1.0,tmp7,shape(4,1,2),peps_i[i],shape(0,1,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[row - 1].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(col)[row],0.0,RO[row - 1]);

            //expectation value:
            val += ham.gcoef(i) * Dot(LOi[i],RO[row-1]);

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity, make intermediary
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[row],shape(0),LOu,shape(0),0.0,tmp6);

         //add peps to it
         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),(*this)(row,col),shape(0,3),0.0,tmp7);

         // add left hamiltonian operator to intermediary to construct next left renormalized operator
         for(int i = 0;i < delta;++i){

            Contract(1.0,ham.gL(i),shape(1),(*this)(row,col),shape(2),0.0,peps_i[i]);

            tmp6.clear();
            Contract(1.0,tmp7,shape(5,0,2),peps_i[i],shape(0,1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            LOi[i].clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LOi[i]);

         }

         // finally construct new left unity
         Contract(1.0,tmp7,shape(0,5,2),(*this)(row,col),shape(0,2,3),0.0,tmp6);
         Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

         LOu.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LOu);

      }

      //last site on the right: close down on the incomings
      for(int i = 0;i < delta;++i){

         //add right hamiiltonian operator to peps
         peps_i[i].clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,col),shape(2),0.0,peps_i[i]);

         //make intermediary
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[Ly - 1],shape(0),LOi[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),(*this)(Ly - 1,col),shape(0,3),0.0,tmp7);

         //contract with right ham operator
         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps_i[i],shape(1,3,0),0.0,tmp6);

         val += ham.gcoef(i) * blas::dot(tmp6.size(), tmp6.data(), 1, env.gr(col)[Ly-1].data(), 1);

      }

   }

   // -- (3) -- || left column = 0: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro('l',*this,R);

   //construct the left operator with two open physical bonds
   tmp7.clear();
   Contract(1.0,env.gr(0)[0],shape(1),(*this)(0,0),shape(4),0.0,tmp7);

   tmp10.clear();
   Contract(1.0,tmp7,shape(1),(*this)(0,0),shape(4),0.0,tmp10);

   tmp10bis.clear();
   Permute(tmp10,shape(0,5,9,2,6,3,7,1,8,4),tmp10bis);

   //contract with left hamiltonian operators to construct left renormalized operators
   for(int i = 0;i < delta;++i){

      Li[i].resize( shape( (*this)(0,0).shape(1),(*this)(0,0).shape(1),env.gr(0)[0].shape(3) ) );

      int m = tmp10.shape(1) * tmp10.shape(3) * tmp10.shape(7);
      int n = tmp10.shape(4) * tmp10.shape(8);

      blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, ham.gL(i).data(), 1, 0.0, Li[i].data(), 1);

   }

   //and contract with peps to get left unit operator
   Lu.resize( shape( (*this)(0,0).shape(1),(*this)(0,0).shape(1),env.gr(0)[0].shape(3) ) );
   blas::gemv(CblasRowMajor,CblasNoTrans, m, n, 1.0, tmp10bis.data(), n, I.data(), 1, 0.0, Lu.data(), 1);

   //middle of the chain:
   for(int row = 1;row < Ly-1;++row){

      //close down the left renormalized operators from the previous site

      //construct the intermediate contraction (paste right)
      tmp5.clear();
      Contract(1.0,env.gr(0)[row],shape(3),R[row-1],shape(2),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,(*this)(row,0),shape(4,1),tmp5,shape(2,4),0.0,tmp6);

      //paste right hamiltonian operators to it for energy 
      for(int i = 0;i < delta;++i){

         peps_i[i].clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,0),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,peps_i[i],shape(0,4,2),tmp6,shape(1,4,5),0.0,tmp5);

         //contract with left S+
         val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp5.data(), 1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,env.gr(0)[row],shape(0),Lu,shape(2),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,(*this)(row,0),shape(3,4),tmp5,shape(4,1),0.0,tmp6);

      // paste left hamiltonian operator to it
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(row,0),shape(2),0.0,peps_i[i]);

         tmp5.clear();
         Contract(1.0,peps_i[i],shape(0,3,4),tmp6,shape(2,5,3),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape((*this)(row,0).shape(1),(*this)(row,0).shape(1),env.gr(0)[row].shape(3)));

      }

      // construct new unity on the left
      Contract(1.0,(*this)(row,0),shape(2,3,4),tmp6,shape(2,5,3),0.0,tmp5);

      Lu = tmp5.reshape_clear(shape((*this)(row,0).shape(1),(*this)(row,0).shape(1),env.gr(0)[row].shape(3)));

   }

   //finally close down on last left site

   //first construct intermediate
   tmp7.clear();
   Contract(1.0,(*this)(Ly-1,0),shape(4),env.gr(0)[Ly - 1],shape(2),0.0,tmp7);

   // paste on right hamiltonian operator to close down left renormalized operator
   for(int i = 0;i < delta;++i){

      peps_i[i].clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,0),shape(2),0.0,peps_i[i]);

      tmp8.clear();
      Contract(1.0,peps_i[i],shape(0,4),tmp7,shape(2,5),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp8.data(),1);

   }

   return val;

}

//forward declarations for types to be used!
template PEPS<double>::PEPS();
template PEPS< complex<double> >::PEPS();

template PEPS<double>::PEPS(int);
template PEPS< complex<double> >::PEPS(int);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,5> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,5> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::scal(double val);
template void PEPS< complex<double> >::scal(complex<double> val);

template void PEPS<double>::load(const char *filename);
template void PEPS< complex<double> >::load(const char *filename);

template void PEPS<double>::save(const char *filename);
template void PEPS< complex<double> >::save(const char *filename);

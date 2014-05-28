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
 * construct an empty PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 */
template<typename T>
PEPS<T>::PEPS() : vector< TArray<T,5> >(Global::lat.gLx() * Global::lat.gLy()) { }

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(Global::lat.gLx() * Global::lat.gLy()) {

   D = D_in;

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();
   int d = Global::lat.gd();

   //corners first

   //r == 0 : c == 0
   (*this)[ Global::lat.grc2i(0,0) ].resize(1,D,d,1,D);

   //r == 0 : c == L - 1
   (*this)[ Global::lat.grc2i(0,Lx-1) ].resize(D,D,d,1,1);

   //r == L - 1 : c == 0
   (*this)[ Global::lat.grc2i(Ly-1,0) ].resize(1,1,d,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ Global::lat.grc2i(Ly-1,Lx-1) ].resize(D,1,d,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ Global::lat.grc2i(0,c) ].resize(D,D,d,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ Global::lat.grc2i(Ly-1,c) ].resize(D,1,d,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ Global::lat.grc2i(r,0) ].resize(1,D,d,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ Global::lat.grc2i(r,Lx - 1) ].resize(D,D,d,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ Global::lat.grc2i(r,c) ].resize(D,D,d,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ Global::lat.grc2i(r,c) ].generate(Global::rgen<T>);

         Normalize((*this)[ Global::lat.grc2i(r,c) ]);
         Scal((T)D,(*this)[ Global::lat.grc2i(r,c) ]);

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

   return (*this)[Global::lat.grc2i(r,c)];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[Global::lat.grc2i(r,c)];

}

/**
 * @return the cutoff virutal dimension
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
 * initialize the peps to an antiferromagnetic D=1 structure, where one time step has been acted on, and compressed to dimension D if necessary
 * @param D compressed dimension of the state
 * @param noise level of noise added to the initial state
 */
template<>
void PEPS<double>::initialize_state(int D_in,double noise) {

   this->D = D_in;

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();
   int d = Global::lat.gd();

   //resize to D=1
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ Global::lat.grc2i(r,c) ].resize(1,1,d,1,1);

         if( (r + c)%2 == 0){

            (*this)[ Global::lat.grc2i(r,c) ](0,0,0,0,0) = 0.0;
            (*this)[ Global::lat.grc2i(r,c) ](0,0,1,0,0) = 1.0;

         }
         else{

            (*this)[ Global::lat.grc2i(r,c) ](0,0,0,0,0) = 1.0;
            (*this)[ Global::lat.grc2i(r,c) ](0,0,1,0,0) = 0.0;

         }

      }

   //now act the trotter gates on perform simple svd

   //first on the rows, i.e. horizontal
   for(int r = 0;r < Ly;++r){

      enum {i,j,k,l,m,n,p};
      IVector<5> pshape;

      //first the even bonds, (0,1)-(2,3),...
      for(int c = 0;c < Lx-1;c+=2){

         //left
         pshape = (*this)[Global::lat.grc2i(r,c)].shape();

         DArray<6> tmp;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(i,j,k,l,m),Trotter::LO,shape(n,p,k),0.0,tmp,shape(i,j,n,l,m,p));

         (*this)[Global::lat.grc2i(r,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1],d,pshape[3],pshape[4]*Trotter::LO.shape(1)));

         //right
         pshape = (*this)[Global::lat.grc2i(r,c+1)].shape();

         Contract(1.0,Trotter::RO,shape(i,j,k),(*this)[Global::lat.grc2i(r,c+1)],shape(l,m,k,n,p),0.0,tmp,shape(l,j,m,i,n,p));

         (*this)[Global::lat.grc2i(r,c+1)] = tmp.reshape_clear(shape(pshape[0]*Trotter::RO.shape(1),pshape[1],d,pshape[3],pshape[4]));

         //now create 'two-site' object
         DArray<8> ts;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(4),(*this)[Global::lat.grc2i(r,c+1)],shape(0),0.0,ts);

         //svd the fucker
         DArray<1> S;
         Gesvd ('S','S', ts, S,(*this)[Global::lat.grc2i(r,c)],(*this)[Global::lat.grc2i(r,c+1)],D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,(*this)[Global::lat.grc2i(r,c+1)]);
         Dimm((*this)[Global::lat.grc2i(r,c)],S);

      }

      //then the odd bonds, (1,2)-(3,4),...
      for(int c = 1;c < Lx-1;c+=2){

         //left
         pshape = (*this)[Global::lat.grc2i(r,c)].shape();

         DArray<6> tmp;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(i,j,k,l,m),Trotter::LO,shape(n,p,k),0.0,tmp,shape(i,j,n,l,m,p));

         (*this)[Global::lat.grc2i(r,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1],d,pshape[3],pshape[4]*Trotter::LO.shape(1)));

         //right
         pshape = (*this)[Global::lat.grc2i(r,c+1)].shape();

         Contract(1.0,Trotter::RO,shape(i,j,k),(*this)[Global::lat.grc2i(r,c+1)],shape(l,m,k,n,p),0.0,tmp,shape(l,j,m,i,n,p));

         (*this)[Global::lat.grc2i(r,c+1)] = tmp.reshape_clear(shape(pshape[0]*Trotter::RO.shape(1),pshape[1],d,pshape[3],pshape[4]));

         //now create 'two-site' object
         DArray<8> ts;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(4),(*this)[Global::lat.grc2i(r,c+1)],shape(0),0.0,ts);

         //svd the fucker
         DArray<1> S;
         Gesvd ('S','S', ts, S,(*this)[Global::lat.grc2i(r,c)],(*this)[Global::lat.grc2i(r,c+1)],D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,(*this)[Global::lat.grc2i(r,c+1)]);
         Dimm((*this)[Global::lat.grc2i(r,c)],S);

      }

   }

   //then on the columns, i.e. vertical bonds
   for(int c = 0;c < Lx;++c){

      enum {i,j,k,l,m,n,p};
      IVector<5> pshape;

      //first the even bonds, (0,1)-(2,3),...
      for(int r = 0;r < Ly-1;r+=2){

         //left
         pshape = (*this)[Global::lat.grc2i(r,c)].shape();
         
         DArray<6> tmp;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(i,j,k,l,m),Trotter::LO,shape(n,p,k),0.0,tmp,shape(i,j,p,n,l,m));

         (*this)[Global::lat.grc2i(r,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1]*Trotter::LO.shape(1),d,pshape[3],pshape[4]));

         //right
         pshape = (*this)[Global::lat.grc2i(r+1,c)].shape();

         Contract(1.0,Trotter::RO,shape(i,j,k),(*this)[Global::lat.grc2i(r+1,c)],shape(l,m,k,n,p),0.0,tmp,shape(l,m,i,n,j,p));

         (*this)[Global::lat.grc2i(r+1,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1],d,pshape[3]*Trotter::RO.shape(1),pshape[4]));

         //now create 'two-site' object
         DArray<8> ts;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(1),(*this)[Global::lat.grc2i(r+1,c)],shape(3),0.0,ts);

         //svd the fucker
         DArray<1> S;
         DArray<5> U;
         DArray<5> VT;

         Gesvd ('S','S', ts, S,U,VT,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(U,S);
         Dimm(S,VT);

         //permute the memory the way it should be
         Permute(U,shape(0,4,1,2,3),(*this)[Global::lat.grc2i(r,c)]);
         Permute(VT,shape(1,2,3,0,4),(*this)[Global::lat.grc2i(r+1,c)]);

      }

      //then the odd bonds, (1,2)-(3,4),...
      for(int r = 1;r < Ly-1;r+=2){

         //left
         pshape = (*this)[Global::lat.grc2i(r,c)].shape();
         
         DArray<6> tmp;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(i,j,k,l,m),Trotter::LO,shape(n,p,k),0.0,tmp,shape(i,j,p,n,l,m));

         (*this)[Global::lat.grc2i(r,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1]*Trotter::LO.shape(1),d,pshape[3],pshape[4]));

         //right
         pshape = (*this)[Global::lat.grc2i(r+1,c)].shape();

         Contract(1.0,Trotter::RO,shape(i,j,k),(*this)[Global::lat.grc2i(r+1,c)],shape(l,m,k,n,p),0.0,tmp,shape(l,m,i,n,j,p));

         (*this)[Global::lat.grc2i(r+1,c)] = tmp.reshape_clear(shape(pshape[0],pshape[1],d,pshape[3]*Trotter::RO.shape(1),pshape[4]));

         //now create 'two-site' object
         DArray<8> ts;
         Contract(1.0,(*this)[Global::lat.grc2i(r,c)],shape(1),(*this)[Global::lat.grc2i(r+1,c)],shape(3),0.0,ts);

         //svd the fucker
         DArray<1> S;
         DArray<5> U;
         DArray<5> VT;

         Gesvd ('S','S', ts, S,U,VT,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(U,S);
         Dimm(S,VT);

         //permute the memory the way it should be
         Permute(U,shape(0,4,1,2,3),(*this)[Global::lat.grc2i(r,c)]);
         Permute(VT,shape(1,2,3,0,4),(*this)[Global::lat.grc2i(r+1,c)]);

      }

   }

   //make some noise!
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         for(int index = 0;index < (*this)(r,c).size();++index)
            (*this)(r,c).data()[index] += noise * Global::rgen<double>();

}

/**
 * @param peps_i peps to take the overlap with
 * @param D_aux auxiliary dimension of the contraction (determines the accuracy of the contraction)
 * @return the inner product of two PEPS <psi1|psi2> 
 */
template<typename T>
T PEPS<T>::dot(const PEPS<T> &peps_i,int D_aux) const {

   //start from bottom
   MPS<T> mps_b('b',*this,peps_i);

   for(int i = 1;i < Global::lat.gLy()/2;++i){

      //i'th row as MPO
      MPO<T> mpo('H',i,*this,peps_i);

      //apply to form MPS with bond dimension D^4
      mps_b.gemv('L',mpo);

      //reduce the dimensions of the edge states using thin svd
      mps_b.cut_edges();

      if(mps_b.gD() > D_aux){

         MPS<T> mps_c(mps_b.size());

         //compress in sweeping fashion
         mps_c.compress(D_aux,mps_b,1);

         mps_b = std::move(mps_c);

      }

   }

   //then from top 
   MPS<T> mps_t('t',*this,peps_i);

   for(int i = Global::lat.gLy() - 2;i >= Global::lat.gLy()/2;--i){

      //i'th row as MPO
      MPO<T> mpo('H',i,*this,peps_i);

      //apply to form MPS with bond dimension D^4
      mps_t.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      mps_t.cut_edges();

      MPS<T> mps_c(mps_t.size());

      //compress in sweeping fashion
      mps_c.compress(D_aux,mps_t,5);

      mps_t = std::move(mps_c);

   }

   return mps_b.dot(mps_t);

}

/** 
 * normalize the peps approximately, using a contraction with auxiliary dimension
 * @param D_aux the auxiliary dimension
 */
template<typename T>
void PEPS<T>::normalize(int D_aux){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   T val = sqrt(this->dot(*this,D_aux));
   val = pow(val,(T)1.0/(T)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(1.0/val,(*this)[ Global::lat.grc2i(r,c) ]);

}

/**
 * scale the peps with a number
 * @param val scalar to be multiplied with the peps
 */
template<typename T>
void PEPS<T>::scal(T val){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   val = pow(val,(T)1.0/(T)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(val,(*this)[ Global::lat.grc2i(r,c) ]);

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

template double PEPS<double>::dot(const PEPS<double> &peps_i,int D_aux) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::normalize(int D_aux);
template void PEPS< complex<double> >::normalize(int D_aux);

template void PEPS<double>::scal(double val);
template void PEPS< complex<double> >::scal(complex<double> val);

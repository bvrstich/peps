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
 * constructor: just sets the length of the vector, nothing is allocates or initialized
 * @param L_in length of the chain
 */
template<typename T>
MPS<T>::MPS(int L_in) : vector< TArray<T,3> >(L_in) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
template<typename T>
MPS<T>::MPS(int L_in,int D_in) : vector< TArray<T,3> >(L_in) {

   D = D_in;

   int d = PEPS<T>::lat.gd();

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i){

      (*this)[i].resize(vdim[i],d,vdim[i+1]);
      (*this)[i].generate(PEPS<T>::rgen);

   }

}

/**
 * construct constructs a standard MPS object, by creating a double layer peps object from the top row
 */
template<typename T>
MPS<T>::MPS(const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,3> >(PEPS<T>::lat.gLx()) {

   D = peps_1.gD() * peps_2.gD();

   //c == 0
   (*this)[0].reshape(shape(1,D,D));

   for(int c = 1;c < this->size() - 1;++c)
      (*this)[c].reshape(shape(D,D,D));

   (*this)[this->size() - 1].reshape(shape(D,D,1));

   enum {i,j,k,l,m,n,o,p,q};

   TArray<T,8> tmp;

   for(int c = 0;c < this->size();++c){

      tmp.clear();

      Contract((T)1.0,peps_1(0,c),shape(i,j,k,l,m),peps_2(0,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

      tmp.move((*this)[c]);

   }

}

/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   D = mps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
MPS<T>::~MPS(){ }

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gD() const {

   return D;

}

/**
 * act with an MPO on this MPS, resulting MPS is returned as *this object
 * @param uplo if == 'U' contract with the upper physical index of the MPO, if == 'L', contract with the lower
 * @param mpo the MPO
 */
template<typename T>
void MPS<T>::gemv(char uplo,const MPO<T> &mpo){

   int L = this->size();

   if(uplo == 'U'){

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[0].reshape(shape(1,D,D*D));

      tmp.move((*this)[0]);

      //middle sites
      for(int c = 1;c < L - 1;++c){

         tmp.clear();

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

         (*this)[c].reshape(shape(D*D,D,D*D));

         tmp.move((*this)[c]);

      }

      //last site
      tmp.clear();

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[L - 1].reshape(shape(D*D,D,1));

      tmp.move((*this)[L - 1]);

   }
   else{//L

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[0].reshape(shape(1,D,D*D));

      tmp.move((*this)[0]);

      //middle sites
      for(int c = 1;c < L - 1;++c){

         tmp.clear();

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

         (*this)[c].reshape(shape(D*D,D,D*D));

         tmp.move((*this)[c]);

      }

      //last site
      tmp.clear();

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[L - 1].reshape(shape(D*D,D,1));

      tmp.move((*this)[L - 1]);

   }

   D *= mpo.gD();

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 */
template<typename T>
void MPS<T>::canonicalize(const BTAS_SIDE &dir){

   if(dir == Left){//QR

      TArray<T,2> R;
      TArray<T,3> tmp;

      for(int i = 0;i < this->size() - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,R,shape(1),(*this)[i + 1],shape(0),(T)0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

   }
   else{//LQ

      TArray<T,2> L;
      TArray<T,3> tmp;

      for(int i = this->size() - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)[i - 1],shape(2),L,shape(0),(T)0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

   }

}

/**
 * find an approximate form of the state 'mps' compressed to a bond dimension 'Dc' by performing an SVD on an non-canonical state.
 * @param dir Left or Right - going compression
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
template<typename T>
void MPS<T>::guess(const BTAS_SIDE &dir,int Dc,const MPS<T> &mps){

   int L = mps.size();

   if(dir == Left){

      TArray<T,3> U;
      TArray<T,2> V;
      TArray<typename remove_complex<T>::type,1> S;

      Gesvd('S','S',mps[0],S,U,V,Dc);

      (*this)[0] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      //and contract V with mps on next site
      Contract((T)1.0,V,shape(1),mps[1],shape(0),(T)0.0,(*this)[1]);

      for(int i = 1;i < L - 1;++i){

         T nrm = sqrt(Dotc((*this)[i],(*this)[i]));
         Scal(1.0/nrm,(*this)[i]);

         this->scal(nrm);

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(U);

         //multiply S to V
         Dimm(S,V);

         //and contract V with mps on next site
         Contract((T)1.0,V,shape(1),mps[i + 1],shape(0),(T)0.0,(*this)[i + 1]);

      }

      T nrm = sqrt(Dotc((*this)[L - 1],(*this)[L - 1]));
      Scal(1.0/nrm,(*this)[L - 1]);

      this->scal(nrm);

   }
   else{

      TArray<T,2> U;
      TArray<T,3> V;
      TArray<typename remove_complex<T>::type,1> S;

      Gesvd('S','S',mps[L - 1],S,U,V,Dc);

      (*this)[L - 1] = std::move(V);

      //multiply U and S
      Dimm(U,S);

      //and contract U with mps on previous site
      Contract((T)1.0,mps[L - 2],shape(2),U,shape(0),(T)0.0,(*this)[L - 2]);

      for(int i = L - 2;i > 0;--i){

         T nrm = sqrt(Dotc((*this)[i],(*this)[i]));
         Scal(1.0/nrm,(*this)[i]);

         this->scal(nrm);

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(V);

         //multiply S to V
         Dimm(U,S);

         //and contract V with mps on next site
         Contract((T)1.0,mps[i - 1],shape(2),U,shape(0),(T)0.0,(*this)[i - 1]);

      }

      T nrm = sqrt(Dotc((*this)[0],(*this)[0]));
      Scal(1.0/nrm,(*this)[0]);

      this->scal(nrm);

   }

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS<double>::scal(double alpha){

   int sign;

   if(alpha > 0)
      sign = 1;
   else
      sign = -1;

   alpha = pow(fabs(alpha),1.0/(double)this->size());

   Scal(sign * alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS< complex<double> >::scal(complex<double> alpha){

   alpha = pow(fabs(alpha),1.0/(complex<double>)this->size());

   Scal(alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

/**
 * find the best compression of the state 'mps' a bond dimension 'Dc' by optimizing the tensor in a sweeping fashion
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
template<typename T>
void MPS<T>::compress(int Dc,const MPS<T> &mps){

   //initial guess by performing svd compression of uncanonicalized state
   this->guess(Right,Dc,mps);

}

template MPS<double>::MPS(const PEPS<double> &,const PEPS<double> &);
template MPS< complex<double> >::MPS(const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template MPS<double>::MPS(int,int);
template MPS< complex<double> >::MPS(int,int);

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::gemv(char uplo,const MPO<double> &mpo);
template void MPS< complex<double> >::gemv(char uplo,const MPO< complex<double> > &mpo);

template void MPS<double>::canonicalize(const BTAS_SIDE &dir);
template void MPS< complex<double> >::canonicalize(const BTAS_SIDE &dir);

template void MPS<double>::guess(const BTAS_SIDE &dir,int Dc,const MPS<double> &mps);
template void MPS< complex<double> >::guess(const BTAS_SIDE &dir,int Dc,const MPS< complex<double> > &mps);

template void MPS<double>::compress(int Dc,const MPS<double> &mps);
template void MPS< complex<double> >::compress(int Dc,const MPS< complex<double> > &mps);

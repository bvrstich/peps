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
 * empty constructor: just sets the length of the vector
 */
template<typename T>
MPS<T>::MPS() : vector< TArray<T,3> >() { }

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

   int d = Global::lat.gd();

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
      (*this)[i].generate(Global::rgen<T>);

   }

}

/**
 * construct constructs a standard MPS object, by creating a double layer peps object
 * @param option construct mps from bottom layer (r == 0) if option == 'b', from top layer (r = Ly-1) if option == 't'
 * from left 'l' c = 0, or right 'r' c = Lx-1
 */
template<typename T>
MPS<T>::MPS(char option,const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,3> >() {

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   if(option == 'b'){

      this->resize(Lx);

      D = peps_1.gD() * peps_2.gD();

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;
 
      for(int c = 0;c < Lx;++c){

         Contract((T)1.0,peps_1(0,c),shape(i,j,k,l,m),peps_2(0,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

         int DL = peps_1(0,c).shape(0) * peps_2(0,c).shape(0);
         int phys_d = peps_1(0,c).shape(1) * peps_2(0,c).shape(1);
         int DR = peps_1(0,c).shape(4) * peps_2(0,c).shape(4);

         (*this)[c] = tmp.reshape_clear(shape(DL,phys_d,DR));

      }

   }
   else if(option == 't'){//top

      this->resize(Lx);

      D = peps_1.gD() * peps_2.gD();

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      for(int c = 0;c < Lx;++c){

         int DL = peps_1(Ly-1,c).shape(0) * peps_2(Ly-1,c).shape(0);
         int phys_d = peps_1(Ly-1,c).shape(3) * peps_2(Ly-1,c).shape(3);
         int DR = peps_1(Ly-1,c).shape(4) * peps_2(Ly-1,c).shape(4);

         Contract((T)1.0,peps_1(Ly-1,c),shape(i,j,k,l,m),peps_2(Ly-1,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

         (*this)[c] = tmp.reshape_clear(shape(DL,phys_d,DR));

      }

   }
   else if(option == 'l'){//left

      this->resize(Ly);

      D = peps_1.gD() * peps_2.gD();

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      for(int r = 0;r < Ly;++r){

         int DL = peps_1(r,0).shape(3) * peps_2(r,0).shape(3);
         int phys_d = peps_1(r,0).shape(4) * peps_2(r,0).shape(4);
         int DR = peps_1(r,0).shape(1) * peps_2(r,0).shape(1);

         Contract((T)1.0,peps_1(r,0),shape(i,j,k,l,m),peps_2(r,0),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,m,q,i,n,j,o));

         (*this)[r] = tmp.reshape_clear(shape(DL,phys_d,DR));

      }

   }
   else{//finally right

      this->resize(Ly);

      D = peps_1.gD() * peps_2.gD();

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      for(int r = 0;r < Ly;++r){

         int DL = peps_1(r,Lx-1).shape(3) * peps_2(r,Lx-1).shape(3);
         int phys_d = peps_1(r,Lx-1).shape(0) * peps_2(r,Lx-1).shape(0);
         int DR = peps_1(r,Lx-1).shape(1) * peps_2(r,Lx-1).shape(1);

         Contract((T)1.0,peps_1(r,Lx-1),shape(i,j,k,l,m),peps_2(r,Lx-1),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,i,n,m,q,j,o));

         (*this)[r] = tmp.reshape_clear(shape(DL,phys_d,DR));

      }

   }

}

/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   this->D = mps_copy.gD();

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

   int DO = mpo.gD();

   int L = this->size();

   if(uplo == 'U'){

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int d_phys = (*this)[0].shape(1);
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }
   else{//L

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int d_phys = (*this)[0].shape(1);
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }

   //vdim is increased
   D *= DO;

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
template<typename T>
void MPS<T>::canonicalize(const BTAS_SIDE &dir,bool norm){

   int length = this->size();

   if(dir == Left){//QR

      TArray<T,2> R;
      TArray<T,3> tmp;

      for(int i = 0;i < length - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,R,shape(1),(*this)[i + 1],shape(0),(T)0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt(Dotc((*this)[length-1],(*this)[length-1]));
         Scal(1.0/nrm,(*this)[length-1]);

      }

   }
   else{//LQ

      TArray<T,2> L;
      TArray<T,3> tmp;

      for(int i = length - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)[i - 1],shape(2),L,shape(0),(T)0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt(Dotc((*this)[0],(*this)[0]));
         Scal(1.0/nrm,(*this)[0]);

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
      (*this)[1].clear();

      Contract((T)1.0,V,shape(1),mps[1],shape(0),(T)0.0,(*this)[1]);

      for(int i = 1;i < L - 1;++i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(U);

         //multiply S to V
         Dimm(S,V);

         //and contract V with mps on next site
         (*this)[i + 1].clear();

         Contract((T)1.0,V,shape(1),mps[i + 1],shape(0),(T)0.0,(*this)[i + 1]);

      }

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
      (*this)[L - 2].clear();

      Contract((T)1.0,mps[L - 2],shape(2),U,shape(0),(T)0.0,(*this)[L - 2]);

      for(int i = L - 2;i > 0;--i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(V);

         //multiply S to V
         Dimm(U,S);

         //and contract V with mps on next site
         (*this)[i - 1].clear();

         Contract((T)1.0,mps[i - 1],shape(2),U,shape(0),(T)0.0,(*this)[i - 1]);

      }

   }

   if(Dc < mps.gD())
      this->D = Dc;
   else
      Dc = mps.gD();

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
void MPS<T>::compress(int Dc,const MPS<T> &mps,int n_iter){

   int L = mps.size();

   //initial guess by performing svd compression of uncanonicalized state: output is right-canonicalized state
   guess(Right,Dc,mps);

   //construct renormalized operators
   std::vector< TArray<T,2> > RO(L - 1);
   std::vector< TArray<T,2> > LO(L - 1);

   compress::init_ro(Right,RO,mps,*this);

   int iter = 0;

   while(iter < n_iter){

      //first site
      (*this)[0].clear();

      Contract((T)1.0,mps[0],shape(2),RO[0],shape(1),(T)0.0,(*this)[0]);

      //QR
      Geqrf((*this)[0],RO[0]);

      //paste to next matrix
      TArray<T,3> tmp;

      Contract((T)1.0,RO[0],shape(1),(*this)[1],shape(0),(T)0.0,tmp);

      (*this)[1] = std::move(tmp);

      compress::update_L(0,LO,mps,*this);

      for(int i = 1;i < L - 1;++i){

         TArray<T,3> I;

         Contract((T)1.0,mps[i],shape(2),RO[i],shape(1),(T)0.0,I);

         (*this)[i].clear();

         Contract((T)1.0,LO[i - 1],shape(1),I,shape(0),(T)0.0,(*this)[i]);

         Geqrf((*this)[i],RO[i]);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,RO[i],shape(1),(*this)[i + 1],shape(0),(T)0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

         compress::update_L(i,LO,mps,*this);

      }

      //and backward!
      (*this)[L - 1].clear();

      Contract((T)1.0,LO[L - 2],shape(1),mps[L - 1],shape(0),(T)0.0,(*this)[L - 1]);

      //LQ
      Gelqf(LO[L - 2],(*this)[L - 1]);

      //paste to next matrix
      tmp.clear();

      Contract((T)1.0,(*this)[L - 2],shape(2),LO[L -  2],shape(0),(T)0.0,tmp);

      (*this)[L - 2] = std::move(tmp);

      compress::update_R(L-1,RO,mps,*this);

      for(int i = L - 2;i > 0;--i){

         TArray<T,3> I;

         Contract((T)1.0,mps[i],shape(2),RO[i],shape(1),(T)0.0,I);

         (*this)[i].clear();

         Contract((T)1.0,LO[i - 1],shape(1),I,shape(0),(T)0.0,(*this)[i]);

         Gelqf(LO[i],(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)[i - 1],shape(2),LO[i],shape(0),(T)0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

         compress::update_R(i,RO,mps,*this);

      }

      ++iter;

   }

   this->D = Dc;

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPS's, with *this being the ket
 */
template<typename T>
T MPS<T>::dot(const MPS<T> &bra) const {

   TArray<T,2> E;

   Contract((T)1.0,bra[0],shape(0,1),(*this)[0],shape(0,1),(T)0.0,E);

   TArray<T,3> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract((T)1.0,bra[i],shape(0),E,shape(0),(T)0.0,I);

      E.clear();

      Contract((T)1.0,I,shape(2,0),(*this)[i],shape(0,1),(T)0.0,E);

   }

   return E(0,0);

}

/**
 * normalize the mps
 * @return the norm before
 */
template<typename T>
T MPS<T>::normalize(){

   T nrm = std::sqrt(this->dot(*this));

   this->scal(1.0/nrm);

   return nrm;

}


/**
 * reduce the dimension of the edge states after MPO action using thin svd.
 */
template<typename T>
void MPS<T>::cut_edges() {

   int L = this->size();

   //Left
   TArray<T,3> U;
   TArray<T,2> V;
   TArray<typename remove_complex<T>::type,1> S;

   int i = 0;

   //easy compression
   while( (*this)[i].shape(0)*(*this)[i].shape(1) < (*this)[i].shape(2) ){

      U.clear();
      S.clear();
      V.clear();

      Gesvd('S','S',(*this)[i],S,U,V);

      (*this)[i] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      U.clear();

      //and contract V with mps on next site
      Contract((T)1.0,V,shape(1),(*this)[i+1],shape(0),(T)0.0,U);

      (*this)[i+1] = std::move(U);

      ++i;

   }

   i = L - 1;

   while( (*this)[i].shape(0) > (*this)[L - 1].shape(1)*(*this)[i].shape(2) ){

      //Right
      U.clear();
      V.clear();
      S.clear();

      Gesvd('S','S',(*this)[i],S,V,U);

      (*this)[i] = std::move(U);

      //multiply U and S
      Dimm(V,S);

      //and contract U with mps on previous site
      U.clear();

      Contract((T)1.0,(*this)[i-1],shape(2),V,shape(0),(T)0.0,U);
      (*this)[i-1] = std::move(U);

      --i;

   }

}

template MPS<double>::MPS(char,const PEPS<double> &,const PEPS<double> &);
template MPS< complex<double> >::MPS(char,const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template MPS<double>::MPS(int,int);
template MPS< complex<double> >::MPS(int,int);

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS();
template MPS< complex<double> >::MPS();

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::gemv(char uplo,const MPO<double> &mpo);
template void MPS< complex<double> >::gemv(char uplo,const MPO< complex<double> > &mpo);

template void MPS<double>::canonicalize(const BTAS_SIDE &dir,bool);
template void MPS< complex<double> >::canonicalize(const BTAS_SIDE &dir,bool);

template void MPS<double>::guess(const BTAS_SIDE &dir,int Dc,const MPS<double> &mps);
template void MPS< complex<double> >::guess(const BTAS_SIDE &dir,int Dc,const MPS< complex<double> > &mps);

template void MPS<double>::compress(int Dc,const MPS<double> &mps,int);
template void MPS< complex<double> >::compress(int Dc,const MPS< complex<double> > &mps,int);

template double MPS<double>::dot(const MPS<double> &bra) const;
template  complex<double>  MPS< complex<double> >::dot(const MPS< complex<double> > &bra) const;

template double MPS<double>::normalize();
template complex<double> MPS< complex<double> >::normalize();

template void MPS<double>::cut_edges();
template void MPS< complex<double> >::cut_edges();

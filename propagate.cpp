#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

namespace propagate {

   /**
    * propagate the peps one imaginary time step
    * @param peps the PEPS to be propagated
    * @param D_aux auxiliary dimension for the contractions. Determines the accuracy of the effective environment.
    */
   void step(PEPS<double> &peps,int D_aux){

      int Lx = Global::lat.gLx();
      int Ly = Global::lat.gLy();

      int d = Global::lat.gd();
      int D = peps.gD();

      enum {i,j,k,m,n,o,p,q};

      // ########################################################### //
      // ########################################################### //
      // ## first propagate applying the gates from bottom to top ## //
      // ########################################################### //
      // ########################################################### //


      // --------------------------------------//
      // --- !!! (1) the bottom row (1) !!! ---// 
      // --------------------------------------//

      //construct the top environment:
      Environment::calc_env('T',peps,D_aux);

      //for the bottom row: construct the right renormalized operator
      vector< DArray<2> > R(Lx - 2);
      DArray<2> L;

      //for this we need the construct the 'bottom environment' for the bottom row:
      Environment::calc_env('B',0,peps,D_aux);

      //first the rightmost operator
      DArray<4> tmp4;
      DArray<3> tmp3;

      //tmp comes out index (t,b)
      Contract(1.0,Environment::t[0][Lx - 1],shape(1),Environment::b[0][Lx - 1],shape(1),0.0,tmp4);

      //reshape tmp to a 2-index array
      R[Lx - 3] = tmp4.reshape_clear(shape(Environment::t[0][Lx - 1].shape(0),Environment::b[0][Lx - 1].shape(0)));

      //now construct the rest
      for(int col = Lx - 2;col > 1;--col){

         tmp3.clear();
         Contract(1.0,Environment::t[0][col],shape(2),R[col-1],shape(0),0.0,tmp3);

         Contract(1.0,tmp3,shape(1,2),Environment::b[0][col],shape(1,2),0.0,R[col-2]);

      }

      //now construct the reduced tensors of the first pair to propagate
      DArray<4> QL;
      DArray<3> a_L;

      //Left
      construct_reduced_tensor('L',peps(0,0),QL,a_L);

      //Right
      DArray<4> QR;
      DArray<3> a_R;

      construct_reduced_tensor('R',peps(0,1),QR,a_R);

      DArray<4> N_eff;
      calc_N_eff('b',0,L,QL,R[0],QR,N_eff);

      //now get the 'X' matrix:
      DArray<3> X;
      get_X(N_eff,X);
      
      //make the environment as 'unitary' as possible
      canonicalize(X,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update(D,a_L,a_R);

      //now expand updated reduced tensors back to the full tensors
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,0),shape(i,j,m,k,n));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(0,1),shape(i,o,j,m,n));

      //finally construct the 'Left renormalized operator' for use on next site

      //first construct a double layer object for the newly updated bottom left site
      Environment::construct_double_layer('H',peps(0,0),Environment::b[0][0]);

      Contract(1.0,Environment::t[0][0],shape(1),Environment::b[0][0],shape(1),0.0,tmp4);

      L = tmp4.reshape_clear(shape(Environment::t[0][0].shape(2),Environment::b[0][0].shape(2)));

      //middle sites of the bottom row:
      for(int col = 1;col < Lx-2;++col){

         //first construct the reduced tensors of the first pair to propagate

         construct_reduced_tensor('L',peps(0,col),QL,a_L);
         construct_reduced_tensor('R',peps(0,col+1),QR,a_R);

         //calculate the effective environment N_eff
         calc_N_eff('b',col,L,QL,R[col],QR,N_eff);

         //extract positive appromixant
         get_X(N_eff,X);

         //make environment close to unitary before the update
         canonicalize(X,a_L,QL,a_R,QR);

         //now do the update! Apply the gates!
         update(D,a_L,a_R);

         //and expand back to the full tensors
         Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,col),shape(i,j,m,k,n));
         Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(0,col+1),shape(i,o,j,m,n));

         //first construct a double layer object for the newly updated bottom 
         Environment::construct_double_layer('H',peps(0,col),Environment::b[0][col]);

         //update the left renormalized operator:
         tmp3.clear();
         Contract(1.0,L,shape(0),Environment::t[0][col],shape(0),0.0,tmp3);

         L.clear();
         Contract(1.0,tmp3,shape(0,1),Environment::b[0][col],shape(0,1),0.0,L);

      }

      //right bottom pair update

      //get the reduced tensors
      construct_reduced_tensor('L',peps(0,Lx-2),QL,a_L);
      construct_reduced_tensor('R',peps(0,Lx-1),QR,a_R);

      //calculate effective environment
      calc_N_eff('b',Lx-2,L,QL,R[Lx-3],QR,N_eff);

      //get positive approximant
      get_X(N_eff,X);

      //make environment close to unitary before the update
      canonicalize(X,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update(D,a_L,a_R);

      //and expand back to the full tensors
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,Lx-2),shape(i,j,m,k,n));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(0,Lx-1),shape(i,o,j,m,n));

      //finally construct a double layer objects for the two new tensors on the bottom
      Environment::construct_double_layer('H',peps(0,Lx-2),Environment::b[0][Lx-2]);
      Environment::construct_double_layer('H',peps(0,Lx-1),Environment::b[0][Lx-1]);
/*
      // ---------------------------------------------------//
      // --- !!! (2) the middle rows (1 -> Ly-2) (2) !!! ---// 
      // ---------------------------------------------------//

      //Right renormalized operators
      vector< DArray<3> > RO(Lx - 2);

      DArray<4> dlo;
      DArray<4> tmp4bis;

      //for(int row = 1;row < Ly-1;++r){
      int row = 1;

      //first create right renormalized operator

      //last site make double layer object from peps
      Environment::construct_double_layer('H',peps(row,Lx-1),dlo);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),dlo,shape(1),0.0,tmp5);

      //then bottom enviroment
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),dlo.shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      //now construct the middle operators
      for(int col = Lx-2;col > 1;--col){

         tmp4.clear();
         Contract(1.0,Environment::t[row][col],shape(2),RO[col-1],shape(0),0.0,tmp4);

         enum {i,j,k,o,m,n};

         Environment::construct_double_layer('H',peps(row,col),dlo);

         tmp4bis.clear();
         Contract(1.0,tmp4,shape(i,j,k,o),dlo,shape(m,j,n,k),0.0,tmp4bis,shape(i,m,n,o));

         RO[col-2].clear();
         Contract(1.0,tmp4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-2]);

      }

      //now construct the reduced tensor and effective environment!

      //Left
      construct_reduced_tensor('L',peps(row,0),QL,a_L);

      //make a 'double layer' object out of Q for contraction with environment
      construct_double_layer('L',QL,tmp5);

      //construct the 'Left' eff environment
      Contract(1.0,Environment::t[row][0],shape(1),tmp5,shape(1),0.0,tmp6);

      DArray<7> tmp7;
      Contract(1.0,tmp6,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp7);

      DArray<4> LO_env = tmp7.reshape_clear(shape(Environment::t[row][0].shape(2),tmp5.shape(3),tmp5.shape(4),Environment::b[row-1][0].shape(2)));

      //Right
      construct_reduced_tensor('R',peps(row,1),QR,a_R);

      //make a 'double layer' object out of Q for contraction with environment
      construct_double_layer('R',QR,tmp5);

      //contract with right renormalized operator:
      tmp4.clear();
      Contract(1.0,Environment::t[row][1],shape(2),RO[0],shape(0),0.0,tmp4);

      //to construct the R_environment
      DArray<5> tmp5bis;
      Contract(1.0,tmp4,shape(i,j,k,m),tmp5,shape(n,o,j,p,k),0.0,tmp5bis,shape(i,n,o,p,m));

      //construct the 'Right' eff environment
      DArray<4> RO_env;
      Contract(1.0,tmp5bis,shape(i,j,k,m,n),Environment::b[row-1][1],shape(p,m,n),0.0,RO_env,shape(i,j,k,p));

      //construct effective environment
      N_eff.clear();
      Contract(1.0,LO_env,shape(i,j,k,m),RO_env,shape(i,n,o,m),0.0,N_eff,shape(j,n,k,o));

      get_X(N_eff,X);

      //now QR and LQ the X matrix and paste it on the aR and aL
      X_copy = X;

      //QR
      Geqrf(X_copy,tmp2);

      //first paste it on the reduced tensor: a_R * R
      tmp3.clear();
      Contract(1.0,a_R,shape(2),tmp2,shape(1),0.0,tmp3);

      a_R = std::move(tmp3);

      //paste the inverse to the environment tensor: R^{-1} * QR
      invert(tmp2);

      tmp4.clear();
      Contract(1.0,tmp2,shape(0),QR,shape(0),0.0,tmp4);

      QR = std::move(tmp4);

      //LQ
      tmp2.clear();
      Gelqf(tmp2,X);

      //first paste it on the reduced tensor: L * a_L
      tmp3.clear();
      Contract(1.0,tmp2,shape(0),a_L,shape(0),0.0,tmp3);

      a_L = std::move(tmp3);

      //paste the inverse to the environment tensor: QL * L^{-1}
      invert(tmp2);

      tmp4.clear();
      Contract(1.0,QL,shape(3),tmp2,shape(1),0.0,tmp4);

      QL = std::move(tmp4);
*/
      //}

   }

   /**
    * construct the left or right reduced tensor form of a peps element by performing QR or LQ decomposition
    * @param L == left, R == right
    */
   void construct_reduced_tensor(char option,const DArray<5> &peps,DArray<4> &Q,DArray<3> &red){

      if(option == 'L'){

         DArray<5> tmp;
         Permute(peps,shape(0,1,3,2,4),tmp);

         int nrows = tmp.shape(0) * tmp.shape(1) * tmp.shape(2);
         int ncols = tmp.shape(3) * tmp.shape(4);

         int min = std::min(nrows,ncols);

         double* tau = new double [min];

         lapack::geqrf(CblasRowMajor,nrows,ncols, tmp.data(), ncols, tau);

         red.resize(shape(min,tmp.shape(3),tmp.shape(4)));
         Q.resize(tmp.shape(0),tmp.shape(1),tmp.shape(2),min);

         //r is in the upper diagonal part of tmp on exit of geqrf:
         for(int i = 0;i < min;++i)
            for(int j = i;j < ncols;++j)
               red.data()[i*ncols + j] = tmp.data()[i*ncols + j];

         //get the input for the Q construction function: lower diagonal part of the matrix
         for(int i = 0;i < nrows;++i)
            for(int j = 0;j < min;++j)
               Q.data()[i*min + j] = tmp.data()[i*ncols + j];

         //now get the Q matrix out
         if(nrows < ncols)
            lapack::orgqr(CblasRowMajor, nrows, nrows, min,Q.data(), nrows, tau);
         else
            lapack::orgqr(CblasRowMajor, nrows, ncols, min,Q.data(), ncols, tau);

         delete [] tau;

      }
      else{//LQ

         DArray<5> tmp;
         Permute(peps,shape(0,2,1,3,4),tmp);

         int nrows = tmp.shape(0) * tmp.shape(1);
         int ncols = tmp.shape(2) * tmp.shape(3) * tmp.shape(4);

         int min = std::min(nrows,ncols);

         double* tau = new double [min];

         lapack::gelqf(CblasRowMajor,nrows,ncols,tmp.data(),ncols,tau);

         red.resize(shape(tmp.shape(0),tmp.shape(1),min));
         Q.resize(min,tmp.shape(2),tmp.shape(3),tmp.shape(4));

         //l is in the lower diagonal part of tmp on exit of geqrf:
         for(int j = 0;j < ncols;++j)
            for(int i = j;i < nrows;++i)
               red.data()[i*min + j] = tmp.data()[i*ncols + j];

         //get the input for the Q construction function: upper diagonal part of the matrix
         for(int i = 0;i < min;++i)
            for(int j = 0;j < ncols;++j)
               Q.data()[i*ncols + j] = tmp.data()[i*ncols + j];

         //now get the Q matrix out
         if(nrows > ncols)
            lapack::orglq(CblasRowMajor,ncols,ncols,min,Q.data(),ncols,tau);
         else
            lapack::orglq(CblasRowMajor,nrows,ncols,min,Q.data(),ncols,tau);

         delete [] tau;

      }

   }

   /**
    * construct a double layer object out of a Q coming from a reduced tensor construction.
    * keep one leg, the one pointing to the reduced tensor, not doubled.
    * @param option 'L'eft or 'R'ight
    * @param Q input object
    * @param dlQ output object
    */
   void construct_double_layer(char option,const DArray<4> &Q,DArray<5> &dlQ){

      //first outer product of Q
      DArray<8> tmp;
      Ger(1.0,Q,Q,tmp);

      if(option == 'L'){

         DArray<8> reorder;
         Permute(tmp,shape(0,4,1,5,2,6,3,7),reorder);

         //and move it to dlQ
         int d0 = reorder.shape(0) * reorder.shape(1);
         int d1 = reorder.shape(2) * reorder.shape(3);
         int d2 = reorder.shape(4) * reorder.shape(5);
         int d3 = reorder.shape(6);
         int d4 = reorder.shape(7);

         dlQ = reorder.reshape_clear(shape(d0,d1,d2,d3,d4));

      }
      else{

         DArray<8> reorder;
         Permute(tmp,shape(0,4,1,5,2,6,3,7),reorder);

         //and move it to dlQ
         int d0 = reorder.shape(0);
         int d1 = reorder.shape(1);
         int d2 = reorder.shape(2) * reorder.shape(3);
         int d3 = reorder.shape(4) * reorder.shape(5);
         int d4 = reorder.shape(6) * reorder.shape(7);

         dlQ = reorder.reshape_clear(shape(d0,d1,d2,d3,d4));

      }

   }

   /**
    * get the X matrix which is the square root of the closest positive approximation to the effective environment
    * @param N_eff input effective environment
    * @param X output DArray<3> is X
    */
   void get_X(DArray<4> &N_eff,DArray<3> &X){

      int matdim = N_eff.shape(0)*N_eff.shape(1);

      //symmetrize
      for(int i = 0;i < matdim;++i)
         for(int j = i + 1;j < matdim;++j){

            N_eff.data()[i*matdim + j] = 0.5 * (N_eff.data()[i*matdim + j]  + N_eff.data()[j*matdim + i]);
            N_eff.data()[j*matdim + i] = N_eff.data()[i*matdim + j];

         }

      DArray<1> eig(matdim);

      //diagonalize
      lapack::syev(CblasRowMajor, 'V','U', matdim, N_eff.data(), matdim, eig.data());

      X.resize(N_eff.shape(0),matdim,N_eff.shape(1));

      //get the square root of the positive approximant:
      for(int iL = 0;iL < N_eff.shape(0);++iL)
         for(int iR = 0;iR < N_eff.shape(1);++iR)
            for(int kL = 0;kL < N_eff.shape(0);++kL)
               for(int kR = 0;kR < N_eff.shape(1);++kR){

                  if(eig(kL*N_eff.shape(1) + kR) > 0.0)
                     X(iL,kL*N_eff.shape(1) + kR,iR) = sqrt( eig(kL*N_eff.shape(1) + kR) ) * N_eff(iL,iR,kL,kR);

               }

   }

   /** 
    * wrapper function invert square general matrix DArray<2>.
    * @param A both input as output matrix: on input A, on output A^{-1}
    */
   void invert(DArray<2> &A){

      int *ipiv = new int [A.shape(0)];

      lapack::getrf(CblasRowMajor,A.shape(0),A.shape(1), A.data(), A.shape(1), ipiv);

      lapack::getri(CblasRowMajor,A.shape(0), A.data(), A.shape(1), ipiv);

      delete [] ipiv;

   }

   /**
    * calculate the effective environment of a pair with the left tensor on site (row,col)
    * @param option 't'op or 'b'ottom row
    * @param col column index of the left tensor
    * @param L left environment matrix
    * @param QL left unitary matrix coming out of the reduced tensor construction
    * @param R left environment matrix
    * @param QR right unitary matrix coming out of the reduced tensor construction
    * @param N_eff output DArray<4> object containing the effective norm environment on exit
    */
   void calc_N_eff(char option,int col,const DArray<2> &L,const DArray<4> &QL,const DArray<2> &R, const DArray<4> &QR,DArray<4> &N_eff){

      int Lx = Global::lat.gLx();
      int Ly = Global::lat.gLy();

      if(option == 'b'){

         if(col == 0){//left edge

            //make a 'double layer' object out of Q for contraction with environment
            DArray<5> tmp5;
            construct_double_layer('L',QL,tmp5);

            //for this one only top contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,Environment::t[0][0],shape(1),tmp5,shape(1),0.0,tmp6);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp6.reshape_clear(shape(Environment::t[0][0].shape(2),tmp5.shape(3),tmp5.shape(4)));

            //make a 'double layer' object out of Q for contraction with environment
            construct_double_layer('R',QR,tmp5);

            //contract with right renormalized operator:
            DArray<3> tmp3;
            Contract(1.0,Environment::t[0][1],shape(2),R,shape(0),0.0,tmp3);

            //to construct the R_environment
            DArray<4> tmp4;
            Contract(1.0,tmp3,shape(1,2),tmp5,shape(2,4),0.0,tmp4);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp4.reshape_clear(shape(tmp4.shape(0),tmp4.shape(1),tmp4.shape(2)));

            //now contract left and right environment to form N_eff
            enum {i,j,k,m,n};
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,n,m),0.0,N_eff,shape(j,n,k,m));

         }
         else if(col == Lx - 2){//right edge

            //Left

            //make a 'double layer' object out of Q for contraction with environment
            DArray<5> tmp5;
            construct_double_layer('L',QL,tmp5);

            //contraction with left renormalized operator
            DArray<3> tmp3;
            Contract(1.0,L,shape(0),Environment::t[0][Lx-2],shape(0),0.0,tmp3);

            DArray<4> tmp4;
            Contract(1.0,tmp3,shape(0,1),tmp5,shape(0,1),0.0,tmp4);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp4.reshape_clear(shape(Environment::t[0][Lx-2].shape(2),tmp5.shape(3),tmp5.shape(4)));

            //Right

            //make a 'double layer' object out of Q for contraction with environment
            construct_double_layer('R',QR,tmp5);

            //only attach to top to construct R_env
            DArray<6> tmp6;
            Contract(1.0,Environment::t[0][Lx-1],shape(1),tmp5,shape(2),0.0,tmp6);

            DArray<3> R_env = tmp6.reshape_clear(shape(Environment::t[0][Lx-1].shape(0),tmp5.shape(0),tmp5.shape(1)));

            //construct effective environment
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,m,n),0.0,N_eff,shape(j,m,k,n));

         }
         else{//middle

            enum {i,j,k,m,n};

            //make a 'double layer' object out of Q for contraction with environment
            DArray<5> tmp5;
            construct_double_layer('L',QL,tmp5);

            //contraction with left renormalized operator
            DArray<3> tmp3;
            Contract(1.0,L,shape(0),Environment::t[0][col],shape(0),0.0,tmp3);

            DArray<4> tmp4;
            Contract(1.0,tmp3,shape(0,1),tmp5,shape(0,1),0.0,tmp4);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp4.reshape_clear(shape(Environment::t[0][col].shape(2),tmp5.shape(3),tmp5.shape(4)));

            //make a 'double layer' object out of Q for contraction with environment
            tmp5.clear();
            construct_double_layer('R',QR,tmp5);

            //contract with right renormalized operator:
            tmp3.clear();
            Contract(1.0,Environment::t[0][col+1],shape(2),R,shape(0),0.0,tmp3);

            //to construct the R_environment
            Contract(1.0,tmp3,shape(1,2),tmp5,shape(2,4),0.0,tmp4);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp4.reshape_clear(shape(tmp4.shape(0),tmp4.shape(1),tmp4.shape(2)));

            //now contract left and right environment to form N_eff
            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,m,n),0.0,N_eff,shape(j,m,k,n));

         }

      }

   }

   /**
    * make the environment as 'canonical' as possible so that the svd for the pair update is as optimal as possible.
    * @param X (XX^T) is the best positive approximation to the environment of the pair
    * @param a_L left reduced tensor, will be multiplied with the L of the environment
    * @param QL unitary part of the left tensor reduction, will be multiplied with the inverse of the L of the environment
    * @param a_R right reduced tensor, will be multiplied with the R of the environment
    * @param QR unitary part of the right tensor reduction, will be multiplied with the inverse of the R of the environment
    */
   void canonicalize(DArray<3> &X,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR){

      //now QR and LQ the X matrix and paste it on the aR and aL
      DArray<3> X_copy(X);

      //QR
      DArray<2> tmp2;
      Geqrf(X_copy,tmp2);

      //first paste it on the reduced tensor: a_R * R
      DArray<3> tmp3;
      Contract(1.0,a_R,shape(2),tmp2,shape(1),0.0,tmp3);

      a_R = std::move(tmp3);

      //paste the inverse to the environment tensor: R^{-1} * QR
      invert(tmp2);

      DArray<4> tmp4;
      Contract(1.0,tmp2,shape(0),QR,shape(0),0.0,tmp4);

      QR = std::move(tmp4);

      //LQ
      Gelqf(tmp2,X);

      //first paste it on the reduced tensor: L * a_L
      tmp3.clear();
      Contract(1.0,tmp2,shape(0),a_L,shape(0),0.0,tmp3);

      a_L = std::move(tmp3);

      //paste the inverse to the environment tensor: QL * L^{-1}
      invert(tmp2);

      tmp4.clear();
      Contract(1.0,QL,shape(3),tmp2,shape(1),0.0,tmp4);

      QL = std::move(tmp4);

   }

   /**
    * update a tensor pair by applying a trotter gate over their mutual bond
    * after which a svd is performed over the bond and the dimensions are set back to D
    */
   void update(int D,DArray<3> &a_L,DArray<3> &a_R){

      enum {i,j,k,m,n};

      //left
      DArray<4> tmp4;
      Contract(1.0,a_L,shape(i,j,k),Trotter::LO,shape(n,m,j),0.0,tmp4,shape(i,n,k,m));

      a_L = tmp4.reshape_clear(shape(a_L.shape(0),a_L.shape(1),a_L.shape(2)*Trotter::LO.shape(1)));

      //right
      tmp4.clear();
      Contract(1.0,Trotter::RO,shape(i,j,k),a_R,shape(n,k,m),0.0,tmp4,shape(n,j,i,m));

      a_R = tmp4.reshape_clear(shape(a_R.shape(0)*Trotter::RO.shape(1),a_R.shape(1),a_R.shape(2)));

      //now create 'two-site' object
      tmp4.clear();
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      //svd the fucker
      DArray<1> S;
      Gesvd ('S','S', tmp4, S,a_L,a_R,D);

      //take the square root of the sv's
      for(int i = 0;i < S.size();++i)
         S(i) = sqrt(S(i));

      //and multiply it left and right to the tensors
      Dimm(S,a_R);
      Dimm(a_L,S);

   }

}

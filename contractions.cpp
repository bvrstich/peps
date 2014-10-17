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
using namespace global;

namespace contractions {

  /**
    * update left renormalized operator on site col 
    * @param option == 't'op ,'b'ottom, 'l'eft or 'r'ight
    * @param rc is row or column index, col for t,b row for r,l
    */
   void update_L(char option,int rc,DArray<3> &L){

      if(option == 'b'){//bottom

         if(rc == 0){

            DArray<4> tmp4;
            Contract(1.0,env.gt(0)[0],shape(1,2),env.gb(0)[0],shape(1,2),0.0,tmp4);

            L = tmp4.reshape_clear(shape(env.gt(0)[0].shape(3),D,D));

         }
         else if(rc < Lx - 1){


            //update the left renormalized operator:
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gt(0)[rc],shape(0),0.0,tmp5);

            DArray<4> tmp4 = tmp5.reshape_clear( shape(tmp5.shape(0)*tmp5.shape(1),tmp5.shape(2),tmp5.shape(3),tmp5.shape(4)) );

            DArray<2> tmp2;
            Contract(1.0,tmp4,shape(0,1,2),env.gb(0)[rc],shape(0,1,2),0.0,tmp2);

            L = tmp2.reshape_clear(shape(env.gt(0)[rc].shape(3),D,D));

         }

      }
      else if(option == 't'){//top
         /*
            if(rc == 0){

            DArray<4> tmp4;
            Contract(1.0,env.gt(Ly-2)[0],shape(1),env.gb(Ly-2)[0],shape(1),0.0,tmp4);

            L = tmp4.reshape_clear(shape(env.gt(Ly-2)[0].shape(2),env.gb(Ly-2)[0].shape(2)));

            }
            else{

         //update the left renormalized operator:
         DArray<3> tmp3;
         Contract(1.0,L,shape(0),env.gt(Ly-2)[rc],shape(0),0.0,tmp3);

         L.clear();
         Contract(1.0,tmp3,shape(0,1),env.gb(Ly-2)[rc],shape(0,1),0.0,L);

         }
         */
      }
      else if(option == 'l'){//left
         /*
            if(rc == 0){

            DArray<4> tmp4;
            Contract(1.0,env.gr(0)[0],shape(1),env.gl(0)[0],shape(1),0.0,tmp4);

            L = tmp4.reshape_clear(shape(env.gr(0)[0].shape(2),env.gl(0)[0].shape(2)));

            }
            else{

         //update the left renormalized operator:
         DArray<3> tmp3;
         Contract(1.0,L,shape(0),env.gr(0)[rc],shape(0),0.0,tmp3);

         L.clear();
         Contract(1.0,tmp3,shape(0,1),env.gl(0)[rc],shape(0,1),0.0,L);

         }
         */
      }
      else{//right
         /*
            if(rc == 0){

            DArray<4> tmp4;
            Contract(1.0,env.gr(Lx-2)[0],shape(1),env.gl(Lx-2)[0],shape(1),0.0,tmp4);

            L = tmp4.reshape_clear(shape(env.gr(Lx-2)[0].shape(2),env.gl(Lx-2)[0].shape(2)));

            }
            else{

         //update the left renormalized operator:
         DArray<3> tmp3;
         Contract(1.0,L,shape(0),env.gr(Lx-2)[rc],shape(0),0.0,tmp3);

         L.clear();
         Contract(1.0,tmp3,shape(0,1),env.gl(Lx-2)[rc],shape(0,1),0.0,L);

         }
         */
      }

   }

   /** 
    * init the right renormalized operator for the middle rows: 
    * @param option 'H'orizontal or 'V'ertical
    * @param rc 'row' index for Horizontal, 'col' index for Vertical
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   void init_ro(char option,int rc,const PEPS<double> &peps,vector< DArray<4> > &RO){

      if(option == 'H'){

         DArray<6> tmp6;
         DArray<6> tmp6bis;
         DArray<7> tmp7;
         DArray<8> tmp8;
         DArray<8> tmp8bis;

         //paste top peps 'operators'
         Contract(1.0,env.gt(rc)[Lx - 1],shape(1),peps(rc,Lx-1),shape(1),0.0,tmp7);

         Contract(1.0,tmp7,shape(1,4),peps(rc,Lx-1),shape(1,2),0.0,tmp8);

         Contract(1.0,tmp8,shape(3,6),env.gb(rc-1)[Lx-1],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object
         RO[Lx - 3] = tmp8bis.reshape_clear(shape(env.gt(rc)[Lx - 1].shape(0),peps(rc,Lx-1).shape(0),peps(rc,Lx-1).shape(0),env.gb(rc-1)[Lx - 1].shape(0)));

         //now construct the middle operators
         for(int col = Lx - 2;col > 1;--col){

            tmp6.clear();
            Contract(1.0,env.gt(rc)[col],shape(3),RO[col-1],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,col),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps(rc,col),shape(1,4,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[col-2].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(rc - 1)[col],0.0,RO[col - 2]);

         }

      }
      else{//vertical, columns

         DArray<6> tmp6;
         DArray<6> tmp6bis;
         DArray<7> tmp7;
         DArray<8> tmp8;
         DArray<8> tmp8bis;

         Contract(1.0,env.gl(rc - 1)[Ly - 1],shape(1),peps(Ly-1,rc),shape(0),0.0,tmp7);

         Contract(1.0,tmp7,shape(1,4),peps(Ly-1,rc),shape(0,2),0.0,tmp8);

         Contract(1.0,tmp8,shape(4,7),env.gr(rc)[Ly-1],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object
         RO[Ly - 3] = tmp8bis.reshape_clear(shape(env.gl(rc - 1)[Ly - 1].shape(0),peps(Ly-1,rc).shape(3),peps(Ly-1,rc).shape(3),env.gr(rc)[Ly - 1].shape(0)));

         //now construct the middle operators
         for(int row = Ly - 2;row > 1;--row){

            tmp6.clear();
            Contract(1.0,env.gl(rc - 1)[row],shape(3),RO[row-1],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(row,rc),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,4),peps(row,rc),shape(0,1,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[row-2].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(rc)[row],0.0,RO[row - 2]);

         }

      }

   }

   /** 
    * init the right renormalized operator for the top or bottom row, or left or right column
    * @param option == 'l'eft 'r'ight 'top' or 'b'ottom
    * @param R vector containing the right operators on exit
    */
   void init_ro(char option,const PEPS<double> &peps,vector< DArray<3> > &R){

      if(option == 'b'){

         //first the rightmost operator
         DArray<4> tmp4;
         DArray<5> tmp5;

         //tmp comes out index (t,b)
         Contract(1.0,env.gt(0)[Lx - 1],shape(1,2),env.gb(0)[Lx - 1],shape(1,2),0.0,tmp4);

         //reshape tmp to a 2-index array
         R[Lx - 3] = tmp4.reshape_clear(shape(env.gt(0)[Lx - 1].shape(0),peps(0,Lx-1).shape(0),peps(0,Lx-1).shape(0)));

         //now construct the rest
         for(int col = Lx - 2;col > 1;--col){

            tmp5.clear();
            Contract(1.0,env.gt(0)[col],shape(3),R[col-1],shape(0),0.0,tmp5);

            R[col - 2].resize(env.gt(0)[col].shape(0),peps(0,col).shape(0),peps(0,col).shape(0));

            int m = tmp5.shape(0);//rows of op(A)
            int n = env.gb(0)[col].shape(0);//col of op(B)
            int k = tmp5.shape(1) * tmp5.shape(2) * tmp5.shape(3) * tmp5.shape(4);

            blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans,m,n,k,1.0, tmp5.data(),k, env.gb(0)[col].data(), k,0.0, R[col - 2].data(), n);

         }

      }
      else if(option == 't'){

         //first the rightmost operator
         DArray<4> tmp4;
         DArray<5> tmp5;

         //tmp comes out index (t,b)
         Contract(1.0,env.gt(Ly-2)[Lx - 1],shape(1,2),env.gb(Ly-2)[Lx - 1],shape(1,2),0.0,tmp4);

         //reshape tmp to a 2-index array
         R[Lx - 3] = tmp4.reshape_clear(shape(peps(Ly-1,Lx-1).shape(0),peps(Ly-1,Lx-1).shape(0),env.gb(Ly-2)[Lx-1].shape(0)));

         //now construct the rest
         for(int col = Lx - 2;col > 1;--col){

            int m = env.gt(Ly - 2)[col].shape(0) * env.gt(Ly - 2)[col].shape(1) * env.gt(Ly - 2)[col].shape(2);//rows of op(A)
            int n = R[col - 1].shape(2);//col of op(B)
            int k = R[col - 1].shape(0) * R[col - 1].shape(1);

            tmp5.resize( shape( peps(Ly-1,col).shape(0),peps(Ly-1,col).shape(0),env.gt(Ly-2)[col].shape(1),env.gt(Ly-2)[col].shape(2),R[col - 1].shape(2) ) );

            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,m,n,k,1.0, env.gt(Ly-2)[col].data(),k, R[col - 1].data(), n,0.0,tmp5.data(), n);

            R[col - 2].clear();
            Contract(1.0,tmp5,shape(2,3,4),env.gb(Ly-2)[col],shape(1,2,3),0.0,R[col - 2]);

         }

      }
      else if(option == 'l'){

         //first the rightmost operator
         DArray<4> tmp4;
         DArray<5> tmp5;

         //tmp comes out index (l,r)
         Contract(1.0,env.gl(0)[Ly - 1],shape(1,2),env.gr(0)[Ly - 1],shape(1,2),0.0,tmp4);

         //reshape tmp to a 2-index array
         R[Ly - 3] = tmp4.reshape_clear( shape(peps(Ly-1,0).shape(3),peps(Ly-1,0).shape(3),env.gr(0)[Ly-1].shape(0)));

         //now construct the rest
         for(int row = Ly - 2;row > 1;--row){

            int m = env.gl(0)[row].shape(0) * env.gl(0)[row].shape(1) * env.gl(0)[row].shape(2);//rows of op(A)
            int n = R[row - 1].shape(2);//col of op(B)
            int k = R[row - 1].shape(0) * R[row - 1].shape(1);

            tmp5.resize( shape( peps(row,Lx-1).shape(3),peps(row,Lx-1).shape(3),env.gl(0)[row].shape(1),env.gl(0)[row].shape(2),R[row-1].shape(2) ) );

            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,m,n,k,1.0, env.gl(0)[row].data(),k, R[row - 1].data(), n,0.0,tmp5.data(), n);

            R[row - 2].clear();
            Contract(1.0,tmp5,shape(2,3,4),env.gr(0)[row],shape(1,2,3),0.0,R[row - 2]);

         }


      }
      else{//right

         //first the rightmost operator
         DArray<4> tmp4;
         DArray<5> tmp5;

         //tmp comes out index (l,r)
         Contract(1.0,env.gl(Lx - 2)[Ly - 1],shape(1,2),env.gr(Lx - 2)[Ly - 1],shape(1,2),0.0,tmp4);

         //reshape tmp to a 2-index array
         R[Ly - 3] = tmp4.reshape_clear(shape(env.gl(Lx - 2)[Ly - 1].shape(0),peps(Ly-1,Lx-1).shape(3),peps(Ly-1,Lx-1).shape(3)));

         //now construct the rest
         for(int row = Ly - 2;row > 1;--row){

            tmp5.clear();
            Contract(1.0,env.gl(Lx - 2)[row],shape(3),R[row-1],shape(0),0.0,tmp5);

            R[row - 2].resize(env.gl(Lx - 2)[row].shape(0),peps(row,Lx-1).shape(3),peps(row,Lx-1).shape(3));

            int m = tmp5.shape(0);//rows of op(A)
            int n = env.gr(Lx - 2)[row].shape(0);//col of op(B)
            int k = tmp5.shape(1) * tmp5.shape(2) * tmp5.shape(3) * tmp5.shape(4);

            blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans,m,n,k,1.0, tmp5.data(),k, env.gr(Lx - 2)[row].data(), k,0.0, R[row - 2].data(), n);

         }

      }

   }


   /**
    * update left renormalized operator on site (row,col )
    * @param option 'H'orizonal or 'V'ertical
    * @param li large index, if 'H' then row, if 'V' then col
    * @param si small index, if 'H' then col, if 'V' then row
    * @param peps the input PEPS object
    * @param LO input old left renormalized operator, output new left renormalized operator
    */
   void update_L(char option,int li,int si,const PEPS<double> &peps,DArray<3> &LO){
      /*
         if(option == 'H'){

         if(si == 0){

         DArray<4> tmp4;
         env.construct_double_layer('H',peps(li,0),tmp4);

      //paste top environment on
      DArray<5> tmp5;
      Contract(1.0,env.gt(li)[0],shape(1),tmp4,shape(1),0.0,tmp5);

      //then bottom enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),env.gb(li-1)[0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LO = tmp6.reshape_clear(shape(env.gt(li)[0].shape(2),tmp4.shape(3),env.gb(li-1)[0].shape(2)));

      }
      else{

      enum {i,j,k,m,n,o};

      //first attach top to left unity
      DArray<4> I4;
      Contract(1.0,env.gt(li)[si],shape(0),LO,shape(0),0.0,I4);

      DArray<4> tmp4;
      env.construct_double_layer('H',peps(li,si),tmp4);

      DArray<4> I4bis;
      Contract(1.0,I4,shape(i,j,k,o),tmp4,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LO.clear();
      Contract(1.0,I4bis,shape(2,3),env.gb(li-1)[si],shape(0,1),0.0,LO);

      }

      }
      else{//Vertical

      if(si == 0){

      DArray<4> tmp4;
      env.construct_double_layer('V',peps(0,li),tmp4);

      //paste top environment on
      DArray<5> tmp5;
      Contract(1.0,env.gr(li)[0],shape(1),tmp4,shape(1),0.0,tmp5);

      //then bottom enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),env.gl(li-1)[0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LO = tmp6.reshape_clear(shape(env.gr(li)[0].shape(2),tmp4.shape(3),env.gl(li-1)[0].shape(2)));

      }
      else{

      enum {i,j,k,m,n,o};

      //first attach top to left unity
      DArray<4> I4;
      Contract(1.0,env.gr(li)[si],shape(0),LO,shape(0),0.0,I4);

      DArray<4> tmp4;
      env.construct_double_layer('V',peps(si,li),tmp4);

      DArray<4> I4bis;
      Contract(1.0,I4,shape(i,j,k,o),tmp4,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LO.clear();
      Contract(1.0,I4bis,shape(2,3),env.gl(li-1)[si],shape(0,1),0.0,LO);

   }

   }
   */
   }

}

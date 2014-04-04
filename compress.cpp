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

namespace compress {

   /**
    * initialize the renormalized operator
    * @param dir left or right renormalized operator
    * @param ro input empty, output will contain left or right ro
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   template<typename T>
      void init_ro(const BTAS_SIDE &dir,std::vector< TArray<T,2> > &ro,const MPS<T> &bra,const MPS<T> &ket){

         int L = bra.size();

         ro.resize(L - 1);

         if(dir == Left){

            Contract((T)1.0,ket[0],shape(0,1),bra[0],shape(0,1),(T)0.0,ro[0]);

            TArray<T,3> I;

            for(int i = 1;i < L - 1;++i){

               I.clear();

               Contract((T)1.0,ket[i],shape(0),ro[i-1],shape(0),(T)0.0,I);

               Contract((T)1.0,I,shape(2,0),bra[i],shape(0,1),(T)0.0,ro[i]);

            }

         }
         else{

            Contract((T)1.0,ket[L - 1],shape(1,2),bra[L - 1],shape(1,2),(T)0.0,ro[L - 2]);

            TArray<T,3> I;

            for(int i = L - 2;i > 0;--i){

               I.clear();

               Contract((T)1.0,ket[i],shape(2),ro[i],shape(0),(T)0.0,I);

               Contract((T)1.0,I,shape(1,2),bra[i],shape(1,2),(T)0.0,ro[i-1]);

            }

         }

      }

   /**
    * update the left renormalized operator
    * @param site index of the site
    * @param LO left renormalized operator
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   template<typename T>
      void update_L(int site,std::vector< TArray<T,2> > &LO,const MPS<T> &bra,const MPS<T> &ket){

         LO[site].clear();

         if(site == 0)
            Contract((T)1.0,ket[0],shape(0,1),bra[0],shape(0,1),(T)0.0,LO[0]);
         else{

            TArray<T,3> I;

            Contract((T)1.0,ket[site],shape(0),LO[site - 1],shape(0),(T)0.0,I);

            Contract((T)1.0,I,shape(2,0),bra[site],shape(0,1),(T)0.0,LO[site]);

         }

      }

   /**
    * update the right renormalized operator
    * @param site index of the site
    * @param RO Right renormalized operator
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   template<typename T>
      void update_R(int site,std::vector< TArray<T,2> > &RO,const MPS<T> &bra,const MPS<T> &ket){

         RO[site - 1].clear();

         int L = bra.size();

         if(site == L - 1)
            Contract((T)1.0,ket[L - 1],shape(1,2),bra[L - 1],shape(1,2),(T)0.0,RO[L - 2]);
         else{

            TArray<T,3> I;

            Contract((T)1.0,ket[site],shape(2),RO[site],shape(0),(T)0.0,I);

            Contract((T)1.0,I,shape(1,2),bra[site],shape(1,2),(T)0.0,RO[site-1]);

         }

      }



   template void init_ro<double>(const BTAS_SIDE &dir,std::vector< TArray<double,2> > &ro,const MPS<double> &bra,const MPS<double> &ket);
   template void init_ro< complex<double> >(const BTAS_SIDE &dir,std::vector< TArray< complex<double> ,2> > &ro,const MPS< complex<double> > &bra,const MPS< complex<double> > &ket);

   template void update_L<double>(int site,std::vector< TArray<double,2> > &ro,const MPS<double> &bra,const MPS<double> &ket);
   template void update_L< complex<double> >(int site,std::vector< TArray< complex<double> ,2> > &ro,const MPS< complex<double> > &bra,const MPS< complex<double> > &ket);

   template void update_R<double>(int site,std::vector< TArray<double,2> > &ro,const MPS<double> &bra,const MPS<double> &ket);
   template void update_R< complex<double> >(int site,std::vector< TArray< complex<double> ,2> > &ro,const MPS< complex<double> > &bra,const MPS< complex<double> > &ket);
}

#ifndef COMPRESS_H
#define COMPRESS_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace compress {

   template<typename T>
      void init_ro(const BTAS_SIDE &,std::vector< TArray<T,2> > &,const MPS<T> &,const MPS<T> &);

   template<typename T>
      void update_L(int,std::vector< TArray<T,2> > &,const MPS<T> &,const MPS<T> &);

   template<typename T>
      void update_R(int,std::vector< TArray<T,2> > &,const MPS<T> &,const MPS<T> &);

}

#endif

#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

#include "MPO.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class MPS is a class written for the construction of matrix product states without symmetry.
 * More specifically it will be used for the contraction of PEPS networks. Where the reduction to a MPS-like for is done.
 */
template<typename T>
class MPS : public vector< TArray<T,3> > {

   public:

      MPS(int L,int D);

      //constructor
      MPS(const PEPS<T> &,const PEPS<T> &);

      //copy constructor
      MPS(const MPS &);

      //destructor
      virtual ~MPS();

      int gD() const;

      void gemv(char , const MPO<T> &);

   private:

      //!dimension of the bonds
      int D;


};

/**
 * output stream operator overloaded for MPS<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const MPS<T> &mps_p){

   for(int s = 0;s < mps_p.size();++s){

         output << std::endl;
         output << "Tensor on site (" << s << ")\t" << std::endl;
         output << std::endl;
         output << mps_p[s] << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/

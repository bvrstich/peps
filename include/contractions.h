#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

//some much repeated contractions are put in separate functions here
namespace contractions {

   void init_ro(bool,char option,const PEPS<double> &,vector< DArray<3> > &R);

   void init_ro(bool,char option,int rc,const PEPS<double> &,vector< DArray<4> > &RO);

   void update_L(char option,int col,const PEPS<double> &,DArray<3> &L);

   void update_L(char option,int row,int col,const PEPS<double> &,DArray<4> &LO);

}

#endif

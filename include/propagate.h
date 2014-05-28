#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &,int);

   void step_no_env(PEPS<double> &,int);

   void construct_reduced_tensor(char,char,const DArray<5> &,DArray<4> &,DArray<3> &);

   void construct_double_layer(char,const DArray<4> &,DArray<5> &);

   void get_X(DArray<4> &,DArray<3> &);

   void invert(DArray<2> &);

   void calc_N_eff(char,int,const DArray<2> &,const DArray<4> &,const DArray<2> &,const DArray<4> &,DArray<4> &);

   void calc_N_eff(char,int,int,const DArray<3> &,const DArray<4> &,const DArray<3> &,const DArray<4> &,DArray<4> &);

   void canonicalize(DArray<3> &X,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR);

   void update(int,DArray<3> &,DArray<3> &);

   void init_ro(char option,vector< DArray<2> > &R);

   void init_ro(char option,int row,const PEPS<double> &,vector< DArray<3> > &R);

   void update_L(char option,int col,DArray<2> &L);

   void update_L(char option,int row,int col,const PEPS<double> &,DArray<3> &LO);

   void apply_stag_field(PEPS<double> &,double B);

}

#endif

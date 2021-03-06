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
 * empty constructor
 */
Hamiltonian::Hamiltonian(){ }

/**
 * copy constructor
 * @param ham_c object to copy
 */
Hamiltonian::Hamiltonian(const Hamiltonian &ham_c){

   this->delta = ham_c.gdelta();

   this->L = ham_c.gL();
   this->R = ham_c.gR();

   this->B = ham_c.gB();

   this->coef = ham_c.gcoef();

   this->is_local = ham_c.gis_local();

}

/**
 * empty destructor
 */
Hamiltonian::~Hamiltonian(){ }

/**
 * @return the number of terms in the nn-interaction
 */
int Hamiltonian::gdelta() const {

   return delta;

}

/**
 * @return the flag, if true there is local interaction, if false not
 */
bool Hamiltonian::gis_local() const {

   return is_local;

}

/**
 * @param i index of the nn term
 * @return coefficient 'i' to the nn-interaction
 */
const double &Hamiltonian::gcoef(int i) const {

   return coef[i];

}

/**
 * @return array of 'delta' coefficients to the nn-interaction
 */
const std::vector<double> &Hamiltonian::gcoef() const {

   return coef;

}

/**
 * @param i index of operator to return
 * @return the left operator with index 'i'
 */
const DArray<2> &Hamiltonian::gL(int i) const {

   return L[i];

}

/**
 * @param i index of operator to return
 * @return the right operator with index 'i'
 */
const DArray<2> &Hamiltonian::gR(int i) const {

   return R[i];

}

/**
 * @return the local operator
 */
const DArray<2> &Hamiltonian::gB() const {

   return B;

}

/**
 * @return the array of 'delta' left operators
 */
const std::vector< DArray<2> > &Hamiltonian::gL() const {

   return L;

}

/**
 * @return the array of 'delta' left operators
 */
const std::vector< DArray<2> > &Hamiltonian::gR() const {

   return R;

}

/**
 * initialize the operators on the nn-Heisenberg model
 * @param ladder if true, use ladder operators (+,-,z), if false, use x,y,z
 */
void Hamiltonian::set_heisenberg(bool ladder) {

   delta = 3;

   is_local = false;

   L.resize(delta);
   R.resize(delta);

   coef.resize(delta);

   if(ladder){

      L[0].resize(global::d,global::d); //S+
      L[1].resize(global::d,global::d); //S-
      L[2].resize(global::d,global::d); //Sz

      R[0].resize(global::d,global::d); //S-
      R[1].resize(global::d,global::d); //S+
      R[2].resize(global::d,global::d); //Sz

      L[0] = 0.0; L[1] = 0.0; L[2] = 0.0;
      R[0] = 0.0; R[1] = 0.0; R[2] = 0.0;

      //S+
      L[0](1,0) = 1.0;
      R[1](1,0) = 1.0;

      //S-
      L[1](0,1) = 1.0;
      R[0](0,1) = 1.0;

      //Sz
      L[2](0,0) = -0.5;
      L[2](1,1) = 0.5;

      R[2](0,0) = -0.5;
      R[2](1,1) = 0.5;

      //coefficients:
      coef[0] = -0.5;//minus sign because of the Marshall sign rule
      coef[1] = -0.5;//minus sign because of the Marshall sign rule
      coef[2] = 1.0;//minus sign because of the Marshall sign rule

   }
   else{

      L[0].resize(global::d,global::d); //Sx
      L[1].resize(global::d,global::d); //iSy
      L[2].resize(global::d,global::d); //Sz

      R[0].resize(global::d,global::d); //Sx
      R[1].resize(global::d,global::d); //iSy
      R[2].resize(global::d,global::d); //Sz

      L[0] = 0.0; L[1] = 0.0; L[2] = 0.0;
      R[0] = 0.0; R[1] = 0.0; R[2] = 0.0;

      //Sx
      L[0](0,1) = 0.5;
      L[0](1,0) = 0.5;

      R[0](0,1) = 0.5;
      R[0](1,0) = 0.5;

      //iSy
      L[1](0,1) = -0.5;
      L[1](1,0) = 0.5;

      R[1](0,1) = -0.5;
      R[1](1,0) = 0.5;

      //Sz
      L[2](0,0) = -0.5;
      L[2](1,1) = 0.5;

      R[2](0,0) = -0.5;
      R[2](1,1) = 0.5;

      //coefficients:
      coef[0] = -1.0;//minus sign because of the Marshall sign rule
      coef[1] = 1.0;//minus sign of the Marshall sign rule cancelled because of squared 'i'
      coef[2] = 1.0;//minus sign because of the Marshall sign rule

   }

}

/**
 * initialize the operators on the nn-Heisenberg model
 * @param B_in strength of the magnetic field in the x direction
 */
void Hamiltonian::set_transverse_field_ising(double B_in) {

   delta = 1;

   is_local = true;

   L.resize(delta);
   R.resize(delta);

   coef.resize(delta);

   L[0].resize(global::d,global::d); //S+
   R[0].resize(global::d,global::d); //S-

   L[0] = 0.0;
   R[0] = 0.0;

   //Sz
   L[0](0,0) = -1.0;
   L[0](1,1) = 1.0;

   R[0](0,0) = -1.0;
   R[0](1,1) = 1.0;

   //coefficients:
   coef[0] = -1.0;

   //local term
   B.resize(global::d,global::d); //Sx

   B = 0.0;

   //Sx
   B(0,1) = B_in;
   B(1,0) = B_in;

}

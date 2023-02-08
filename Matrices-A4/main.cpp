//============================================================
// file: main.cpp
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>

//MATRIX and VECTOR classes assignment
#include "vector3dT_soln.h"
// #include "matrix3dT.h"

// void test_matrices_and_vectors() {
//   print("\n====================  TESTING MATRICES and VECTORS  ========================");
//   vector3D p("p", 2, {1, 2});
//   matrix3dD m("m", 2, {1, 2,   3, 4});
//   show_vect(p);
//   show_mat(m);
//   assert(p * m == m * p);

//   vector3D q("q", 3, {1, 2, 3});
//   matrix3dD n("n", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
//   show_vect(q);
//   show_mat(n);
//   assert(q * n == n * q);
//   print("...test_matrices_and_vectors assertions passed");
//   print("====================  FINISHED testing matrices and vectors  ========================");
// }


int main(int argc, const char * argv[]) {
  vector3D::run_tests();
  // test_matrices();
  // test_matrices_and_vectors();
    
  return 0;
}



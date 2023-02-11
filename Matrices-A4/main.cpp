//============================================================
// file: main.cpp
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>

//MATRIX and VECTOR classes assignment
#include "vector3dT.h"
#include "matrix3dT.h"

int main(int argc, const char * argv[]) {
  vector3D::run_tests();
  matrix3D::run_tests(); 
    
  return 0;
}



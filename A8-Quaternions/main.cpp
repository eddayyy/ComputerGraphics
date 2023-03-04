
#include "quaternion_T.h"

void test_quaternions() {
    print("\n====================  TESTING QUATERNIONS  ========================");
    quaternion<double>::run_tests();
    print("...test_matrices_and_vectors assertions passed");
    print("====================  FINISHED testing quaternions  ========================");
}

int main(int argc, const char * argv[]) {
//  test_vectors();
//  test_matrices();
//  test_matrices_and_vectors();
  test_quaternions();
  print("... program completed...\n");
  return 0;
}
// Author: Eduardo Nunez
// Author email: eduardonunez@csu.fullerton.edu

#include "quaternion_T.h"
#include <stdio.h>

void test_quaternions() {
    printf("\n====================  TESTING QUATERNIONS  ========================");
    quaternion<double>::run_tests();
    printf("...test_matrices_and_vectors assertions passed");
    printf("====================  FINISHED testing quaternions  ========================");
}

int main(int argc, const char * argv[]) {
    FILE *fp;
    fp = freopen("output.txt", "w", stdout); // redirect stdout to a file
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    test_quaternions();
    printf("... program completed...\n");
    fclose(fp);
    return 0;
}
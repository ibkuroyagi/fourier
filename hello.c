#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double arr_test(double arr1[],double arr2[],int N){
  for(int i=0;i<N;i++){
    arr2[i] = arr1[i]*2;
  }
}
int main(void) {
  double arr1[16],arr2[16];
  double a[16][16];
  // for(int i=0;i<16;i++){
  //   arr1[i] = i;
  // }

  // arr_test(arr1,arr2,16);
  // for(int i=0;i<16;i++){
  //   printf("arr1[%d] = %lf, arr2[%d] = %lf\n",i,arr1[i],i,arr2[i]);
  // }
  
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      // a[i][j] = 10*i + j;
      a[i][j] = 1;
    }    
  }
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      printf("%d, ",a[i][j]);
    }
    printf("\n");
  }
  return 0;
}


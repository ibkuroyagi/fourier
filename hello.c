#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double arr_test(double arr1[],double arr2[],int N){
  for(int i=0;i<N;i++){
    arr2[i] = arr1[i]*2+0.01;
  }
}

double make_arr(double arr1[],double arr2[]){
  double a[16][16];
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      a[i][j] = 20.0*(double)i + (double)j + 0.101;
    }
  }  
  for(int k=0;k<10;k++){
    arr2[k] = a[k][1];
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
      a[i][j] = 10.0*(double)i + (double)j;
      // a[i][j] = 1;
    }    
  }
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      printf("%lf, ",a[i][j]);
    }
    printf("\n");
  }
  make_arr(arr1,arr2);
  printf("%lf\n",arr2[1]);
  return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void open_1(void);
void line(void){
  int i;
  for(i=0;i<10;i++){
    printf("no\n");
  }
}

int main(void){
  printf("a1\n");
  open_1();
  line();
  printf("a2\n");
  return 0;
}

void open_1(void) {
	FILE *fp; // FILE型構造体
  int i = 0;
	char fname[] = "test.csv";
	float f1, Re[1024];
  printf("aaa\n");
	fp = fopen(fname, "r"); // ファイルを開く。失敗するとNULLを返す。
	if(fp == NULL) {
		printf("%s file not open!\n", fname);
		exit(1);
	}
 
	while(fscanf(fp, "%f", &f1) != EOF) {
		printf("%f\n", f1);
    Re[i] = f1;
    i++;
	}
  i = i - 1;
  printf("Re[%d] = %f\n",i,Re[i]);
	fclose(fp); // ファイルを閉じる 
}
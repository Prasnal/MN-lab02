#include <iostream>
#include <math.h>
#include "nr.h"
#include "nrutil.c"
#include "nrutil.h"
#include "ludcmp.c"
#include "lubksb.c"
#define SIZE 3
using namespace std;

void print_matrix(int high,int len, float **tab){

  for(int i=1; i<high+1; i++){
    for(int j=1; j<len+1; j++){
      std::cout<<tab[i][j]<<" ";
    }
    std::cout<<""<<std::endl;
  }
}

float max_el(int high, int len, float**tab){
  float max=0.;
  for (int i=1; i<high+1; i++){
    for (int j=1; j<len+1; j++){
      if(fabs(tab[i][j])>max){
        max=tab[i][j];
      }
    }
  }
  return max;
}

float** multi_matrix(float **mat_one, float** mat_two, int row){
  float** result;
  result=matrix(1,3,1,3);
  for(int i=1;i<row+1; i++){
    for(int j=1;j<row+1; j++){
      for(int k=1;k<row+1; k++){
        result[i][j]=result[i][j]+mat_one[i][k]*mat_two[k][j];
      }
    }
  }
  return result;
}

float** reverse(float **matri, int size, int* indexA){ // DZIALA
  float** result=matrix(1,size,1,size);
  float* wektor=vector(1,size);
  //int* indexA=ivector(1,size);
  for(int j=1; j<size+1;j++){
      for(int i=1;i<size+1;i++){
      wektor[i]=0.0;
    }
    wektor[j]=1.0;  
    lubksb(matri,size,indexA,wektor);
    for(int i=1;i<size+1;i++){
      result[i][j]=wektor[i];
      }
  } 
  return result;
}

int main(){
  FILE *fp=fopen("dane.dat","w");
  float **A;
  float **B;
  float **Acp;
  float **Bcp;
  float **Aodw;
  float **Bodw;

  float **multimatrixA;
  float **multimatrixB;
  
  int n=SIZE;
  int* indexA=ivector(1,n); 
  int* indexB=ivector(1,n); 
  float da; 
  float db;
  float* ba1=vector(1,n);
  float* bb1=vector(1,n);
  int tmp=1;
  
  //Tworzenie macierzy A i B oraz ich kopii
  A=matrix(1,n,1,n);
  B=matrix(1,n,1,n);
  Acp=matrix(1,n,1,n);
  Bcp=matrix(1,n,1,n);
  Aodw=matrix(1,n,1,n);
  Bodw=matrix(1,n,1,n);
  multimatrixA=matrix(1,n,1,n);
  multimatrixB=matrix(1,n,1,n);
  for(int i=1; i<n+1; i++){
    for(int j=1; j<n+1; j++){
      A[i][j]=tmp;
      B[i][j]=tmp;
      Acp[i][j]=tmp;
      Bcp[i][j]=tmp;
      tmp++;
    }
  }
  B[1][1]=1.1;
  Bcp[1][1]=1.1;
  //============================
  
  //MACIERZ A:
  printf("Macierz A: \n");
  print_matrix(n,n,A);
 
  //MACIERZ B:
  printf("Macierz B: \n");
  print_matrix(n,n,B);  
  
  //rozklad LU macierzy A
  ludcmp(A,n,indexA,&da);
  printf("Rozklad LU macierzy A: \n");
  print_matrix(n,n,A);
  //printf("WYZNACNZIK A:");
  //cout<<wyznacznik(A,n)<<endl;
  
  //ZAPISYWANIE ROZKLADU LU DO PLIKU
  fprintf(fp,"ROZKLAD LU MACIERZY A: \n");
  
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
       fprintf(fp," %f ",A[i][j]);
    }
    fprintf(fp,"\n");
    }

  //rozklad LU macierzy B
  ludcmp(B,n,indexB,&db);
  printf("Rozklad LU macierzy B: \n");
  print_matrix(n,n,B);
  //printf("WYZNACNZIK B:");
  //cout<<wyznacznik(B,n)<<endl;

  //ZAPISYWANIE ROZKLADU LU DO PLIKU
  fprintf(fp,"ROZKLAD LU MACIERZY B: \n");
  
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
      fprintf(fp," %f ",B[i][j]);
    }
    fprintf(fp,"\n");
  }
  
  //Odwracanie macierzy A
  printf("A^-1: \n");
  Aodw=reverse(A,n,indexA);
  print_matrix(n,n,Aodw);


  //ZAPISYWANIE ODWROCONEJ MACIERZY A
  fprintf(fp,"A^(-1):: \n");
  
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
       fprintf(fp," %f ",Aodw[i][j]);
    }
    fprintf(fp,"\n");
    }
  
  //Odwracanie macierzy B
  printf("B^ -1: \n");
  Bodw=reverse(B,n,indexB);
  print_matrix(n,n,Bodw);

  //ZAPISYWANIE ODWROCONEJ MACIERZY A
  fprintf(fp,"B^(-1):: \n");
  
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
       fprintf(fp," %f ",Bodw[i][j]);
    }
    fprintf(fp,"\n");
  }
  
  //Wskazniki uwarunkowania macierzy A

  printf("Macierz A: \n");
  print_matrix(n,n,Acp);
  cout<<"max el macierzy A: "<<max_el(n,n,Acp)<<endl;
  cout<<"max el macierzy A^-1: "<<max_el(n,n,Aodw)<<endl;
  cout<<"wskaznik uwarunkowania A: "<<max_el(n,n,Acp)*max_el(n,n,Aodw)<<endl;
  //zapis wskaznika uwarunkowania do pliku:
  fprintf(fp,"Wskaznik uwarunkowania macierzy A: %f \n",max_el(n,n,Acp)*max_el(n,n,Aodw));
  
  //Wskaznik uwarunkowania macierzy B 
  printf("Macierz B: \n");
  print_matrix(n,n,Bcp);
  cout<<"max element macierzy B: "<<max_el(n,n,Bcp)<<endl;
  cout<<"max element macierzy B^-1: "<<max_el(n,n,Bodw)<<endl;
  cout<<"wskaznik uwarunkowania B: "<<max_el(n,n,Bcp)*max_el(n,n,Bodw)<<endl;
  fprintf(fp,"Wskaznik uwarunkowania macierzy B: %f \n",max_el(n,n,Bcp)*max_el(n,n,Bodw));


  //Iloczyn macierzy A*A-1
  print_matrix(n,n,Acp);
  printf("* \n");
  print_matrix(n,n,Aodw);
  printf("A*A^(-1)= \n");
  multimatrixA=multi_matrix(Acp,Aodw,n);
  fprintf(fp,"MACIERZ A*A^(-1): \n");
  
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
       fprintf(fp," %f ",multimatrixA[i][j]);
    }
    fprintf(fp,"\n");
    }

  print_matrix(n,n,multimatrixA);
  //Iloczyn macierzy B*B-1
  print_matrix(n,n,Bcp);
  printf("* \n");
  print_matrix(n,n,Bodw);
  printf("\n B*B^(-1)= \n");
  multimatrixB=multi_matrix(Bcp,Bodw,n);

  fprintf(fp,"MACIERZ B*B^(-1): \n");
  for(int i=1; i<n+1; i++){ //zapisywanie macierzy do pliku
    for(int j=1; j<n+1; j++){
       fprintf(fp," %f ",multimatrixB[i][j]);
    }
    fprintf(fp,"\n");
    }

  fclose(fp);
  return 0;
}

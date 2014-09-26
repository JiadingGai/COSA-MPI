#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "data_io.h"
#include "matrix_util.h"

typedef struct{
    int objectID;
    float distance;
} NN_T; // nearest neighbor struct

float Calc_Weight_Error(MATRIX_T* w, MATRIX_T* w_old);
float max(float a, float b);

int N;
int n;
static int Noit;
static int Niter;
static int K;
static float lambda, alpha, threshold;
int main(int argc, char* argv[])
{
  char* data_file_name;
  if(argc == 10)
  {
    data_file_name = argv[1];
    K = atoi( argv[2]  );
    lambda = atof( argv[3] );
    alpha = atof( argv[4] );
    Niter = atoi( argv[5] );
    Noit  = atoi( argv[6] );    
    N = atoi( argv[7] );
    n = atoi( argv[8] );
    threshold = atof( argv[9] );
  }
  else
  {
	printf("Usage: ./cosa2 <data_file_name> <K> <lambda> <alpha> <Niter> <Noit> <N> <n> <threshold>\n");
	printf("  E.g: ./cosa2 c22_ceu.dat 1 0.2 0.1 6 3 1000 13000 0.1\n");
	printf("  E.g: ./cosa2 e6.dat 1 0.2 0.1 6 3 7 10 0.1\n");
	exit(1);
  }

  //Generate_random_data_text("e8.txt",N,n);
  //ConvertData_text2binary("e8.txt", "e8.dat",N,n,1);
 
  clock_t start = clock();

  // IO
  FILE *input_file;
  input_file = fopen(data_file_name,"rb");

  MATRIX_T* x = Matrix_Allocate(n,N,1);// x is an N-by-n data matrix
  int i, j, k;
  
  for(i=0;i<N;i++)
  {
      for(j=0;j<n;j++)
      {
          fread(&Entry2D(x,i,j), sizeof(float), 1, input_file );
          //printf("%f ",Entry2D(x,i,j));
      }
      //printf("\n");
  }
  fclose(input_file);
         
  // Initialization
  //int K = (int) sqrt((float) N); // # of KNN; Ref: COSA paper, section 7, last paragraph
  MATRIX_T* w = Matrix_Allocate(N, n, 1);// w is an n-by-N matrix
  MATRIX_T* w_old = Matrix_Allocate(N, n, 1);// w_old is an n-by-N matrix
  MATRIX_T* s = Matrix_Allocate(1, n, 1);// s is an n-by-1 vector
  MATRIX_T* D = Matrix_Allocate(N, N, 1);// D is an N-by-N distance matrix between all pairs
  //MATRIX_T* d = Matrix_Allocate(N, N, n);// d is an N-by-N-by-n distance matrix per attribute
  NN_T* knn = (NN_T*) malloc(N * K * sizeof(NN_T));// knn is an N-by-K matrix (i.e., height=N; width=K).
  MATRIX_T* S = Matrix_Allocate(N, n, 1);// S is an n-by-N matrix
  MATRIX_T* Sum = Matrix_Allocate(1, N, 1);// Sum is an N-by-1 vector
  //float lambda = 0.2, alpha = 0.10;// COSA paper, section 6
  float eta = lambda;// COSA paper
  int noit = 0, niter = 0;// outer iterations
  float weight_error = FLT_MAX;// Initialize to a huge number; Average error.
  //float threshold = 1e-5;// w stablizes when (weight_error < threshold)
  char wtcomb[] = "each";// weight combining (30) or (33)

  for(k=0;k<n;k++)
  for(i=0;i<N;i++)
    Entry2D(w,k,i) = 1.0 / ((float) n); // 1.0;// weights add up to the number of attribute

  while( noit<Noit && weight_error>threshold )
  {
      printf("Outer iteration# %d (weight_error=%f)\n", noit, weight_error);
      noit++;
      Copy_Matrix(w, w_old);

      niter = 0;
      while(niter<Niter)
      {
          printf("Inner iteration# %d\n ",niter);
          niter++;

          for(i=0;i<N;i++)
          for(j=0;j<N;j++)
              Entry2D(D,i,j) = 0.0;
          
          for(k=0;k<n;k++)
          {
              Entry1D(s,k) = 0.0;
              for(i=0;i<N;i++)
              for(j=0;j<N;j++)
                  Entry1D(s,k) += fabs(Entry2D(x,i,k) - Entry2D(x,j,k));
              Entry1D(s,k) /= ((float) N*N);

              for(i=0;i<N;i++)
              for(j=0;j<N;j++)
              {
                  //Entry3D(d,i,j,k) = fabs(Entry2D(x,i,k) - Entry2D(x,j,k)) / Entry1D(s,k);
                  float dijk = fabs(Entry2D(x,i,k) - Entry2D(x,j,k)) / Entry1D(s,k);
                  Entry2D(D,i,j) += Entry2D(w,k,i) * exp( -dijk / eta );
              }
          }
          
          for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
              //if(i==j) 
              //    Entry2D(D,i,j) = 0.0;
              //else
              Entry2D(D,i,j) = -eta * log( Entry2D(D,i,j) );                
          }
          
          for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
              if(!strcmp(wtcomb,"each"))
                      Entry2D(D,i,j) = max( Entry2D(D,i,j),Entry2D(D,j,i) );
              else
              {
                  Entry2D(D,i,j) = 0.0;
                  for(k=0;k<n;k++)
                  {
                      float dijk = fabs(Entry2D(x,i,k) - Entry2D(x,j,k)) / Entry1D(s,k);
                      Entry2D(D,i,j) += max( Entry2D(w,k,i),Entry2D(w,k,j) ) * dijk;
                  }
              }
          }
      
          for(i=0;i<N;i++)
          {
              for(j=0;j<N;j++)
              {
                  if(j<K)
                  {
                      (*(knn+i*K+j)).objectID = j;
                      (*(knn+i*K+j)).distance = Entry2D(D,i,j);                    
                  }
                  else
                  {
                      // Compare Entry2D(D,i,j) with the max distance
                      // among the current K nearest neighbors.
                      int max_pos = -1;
                      float max_distance = -1.0*FLT_MAX;
                      for(k=0;k<K;k++)
                      {
                          if( max_distance<(*(knn+i*K+k)).distance )
                          {
                              max_pos = k;
                              max_distance = (*(knn+i*K+k)).distance;
                          }
                      }
                      if( Entry2D(D,i,j)<max_distance )
                      {
                          (*(knn+i*K+max_pos)).objectID = j;
                          (*(knn+i*K+max_pos)).distance = Entry2D(D,i,j);
                      }        
                  }
              }
          }
          
          for(k=0;k<n;k++)
          for(i=0;i<N;i++)
          {
              Entry2D(S,k,i) = 0.0;
              for(j=0;j<K;j++)
              {
                  float dijk = fabs(Entry2D(x,i,k) - Entry2D(x,(*(knn+i*K+j)).objectID,k)) / Entry1D(s,k);
                  Entry2D(S,k,i) += dijk / ((float) K);
              }
              Entry2D(w,k,i) = exp( -Entry2D(S,k,i) / lambda );
          }

          for(i=0;i<N;i++)// Normalize w[k][i]
          {
              Entry1D(Sum,i) = 0.0;
              for(k=0;k<n;k++)
                  Entry1D(Sum,i) += Entry2D(w,k,i);

              for(k=0;k<n;k++)
                  Entry2D(w,k,i) /= Entry1D(Sum,i);

          }
      }
      eta = eta + alpha * lambda;
      weight_error = Calc_Weight_Error(w,w_old);
  }

  // IO
  FILE *dist_file;
  dist_file = fopen("mydist.txt","w");
  
  for(j=0;j<N;j++)
  {
    for(i=0;i<N;i++)
    {
      //if(i>j) // lower half
        fprintf(dist_file, "%f,",Entry2D(D,i,j));
    }
    fprintf(dist_file,"\n");
  }
  fclose(dist_file);

  FILE *weight_file;
  weight_file = fopen("myweight.txt","w");

  for(k=0;k<n;k++)
  {
    for(i=0;i<N;i++)
    {
      float temp = Entry2D(w,k,i) * ((float) n);
      fprintf(weight_file, "%f\t",temp);
    }
      fprintf(weight_file,"\n");
  }
  fclose(weight_file);

  // take timing
  clock_t end = clock();
  FILE* time_file = fopen("time_serial_cosa2.txt","w");
  fprintf(time_file,"Elapsed time = %f seconds\n",(end-start)/(float) CLOCKS_PER_SEC);
  fclose(time_file);

  Free_Matrix(x);
  Free_Matrix(w);
  Free_Matrix(w_old);
  Free_Matrix(s);
  Free_Matrix(D);
  //Free_Matrix(d);
  free(knn);
  Free_Matrix(S);
  Free_Matrix(Sum);
  return 1;
}

float Calc_Weight_Error(MATRIX_T* w, MATRIX_T* w_old)
{
    float weight_error = 0.0f;
    int i, j;
    
    for (i = 0; i < (w->height); i++)
    for (j = 0; j < (w->width); j++)
      weight_error += fabs(Entry2D(w,i,j) - Entry2D(w_old,i,j));

    weight_error = weight_error / (float)(w->height * w->width);

    return weight_error;
}
float max(float a, float b)
{
    float bigger = a > b ? a : b;
    return bigger;
}

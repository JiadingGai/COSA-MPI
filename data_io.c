#include "data_io.h"
#include "matrix_util.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// Convert a data file (text format) "in" to "out" (binary format)
void ConvertData_text2binary(const char* in, const char* out, int N, int n, int AutoScale)
{
    FILE* in_file  = fopen(in,"r");
    FILE* out_file = fopen(out,"wb");

    int i,j;
    float temp;

    if(AutoScale){
        int robust = 0;
        MATRIX_T *x = Matrix_Allocate(n,N,1);
        for(i=0;i<N;i++)
        {
            for(j=0;j<n;j++)
            {
                fscanf( in_file,"%f\t",&temp );
                Entry2D(x,i,j) = temp;
            }
        }

        Auto_Scale(x,robust);

        for(i=0;i<N;i++)
        {
            for(j=0;j<n;j++)
            {
                //if(j<1055)
                    fwrite( &Entry2D(x,i,j),sizeof(float),1,out_file );
            }
        }
        Free_Matrix(x);
    }
    else{

        for(i=0;i<N;i++)
        {
            for(j=0;j<n;j++)
            {
                fscanf( in_file,"%f\t",&temp );
                //printf("%f ",temp);
                fwrite( &temp,sizeof(float),1,out_file );
            }
            //printf("\n");
        }

    }

    fclose(in_file);
    fclose(out_file);
}

// Generate random test data "file_name" of Nxn (text format)
void Generate_random_data_text( const char* file_name, int N, int n  )
{
  FILE* file = fopen(file_name, "w");
  int i;
  float temp;

  srand( time(NULL) );
  for(i=0;i<N*n;i++)
  {
    temp = (float) ( 99.0*rand() / RAND_MAX );
    fprintf(file,"%f\t",temp);
  }
  fclose(file);
}

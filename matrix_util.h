#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#define Entry1D(b,i) (*((b->entries)+i))
#define Entry2D(x,i,j)   (*((x->entries)+i*(x->width)+j))
#define Entry3D(d,i,j,k) (*((d->entries)+k*(d->width)*(d->height)+i*(d->width)+j))

typedef struct
{
  int width;
  int height;
  int depth;
  float* entries;
} MATRIX_T;

MATRIX_T* Matrix_Allocate(int width, int height, int depth);
void Free_Matrix(MATRIX_T* a_matrix);
void Copy_Matrix(MATRIX_T* src, MATRIX_T* dest);
void write_matrix_to_file(MATRIX_T* a_matrix, const char* filename);
void print_matrix(MATRIX_T* a_matrix);
float GetMedian(float *a, int p, int q);
float IQR(float *a, int p, int q);
float GetMean(float *a, int p, int q);
float GetVar(float *a, int p, int q);
void Auto_Scale(MATRIX_T* x, int robust);
#endif

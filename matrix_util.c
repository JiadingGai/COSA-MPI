#include "matrix_util.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*
 * Local/Static function declarations go inside the source
 * and not the header file!
 *
 */
static int partition(float *a, int p, int q);
static int RANDOMIZED_partition(float *a, int p, int q);
static float RANDOMIZED_select(float *a, int p, int q, int i);
static void swap(float *a, float *b);
static void RANDOMIZED_quick_sort(float *a, int p, int q);


MATRIX_T* Matrix_Allocate(int width, int height, int depth)
{
    MATRIX_T* a_matrix = (MATRIX_T*) malloc(sizeof(MATRIX_T));
    a_matrix->width  = width;
    a_matrix->height = height;
    a_matrix->depth  = depth;

    a_matrix->entries = (float*) malloc(width*height*depth*sizeof(float));    
    
    return a_matrix;
}
void Free_Matrix(MATRIX_T* a_matrix)
{
    free(a_matrix->entries);
    free(a_matrix);
}
void Copy_Matrix(MATRIX_T* src, MATRIX_T* dest)
{//ASSUMPTION: src and dest have been properly
 //            initialized to the same size.
 
    int i,j,k;

    if(1==src->depth && 1==src->width)
    {// 1D column vector duplication
        for(i=0;i<(src->height);i++)
            Entry1D(dest,i) = Entry1D(src,i);
    }
    else if(1==src->depth && 1<src->width)
    {// 2D matrix duplication
        for(i=0;i<(src->height);i++)
        for(j=0;j<(src->width);j++)
            Entry2D(dest,i,j) = Entry2D(src,i,j);    
    }
    else
    {// 3D matrix duplication
        for(i=0;i<(src->height);i++)
        for(j=0;j<(src->width);j++)
        for(k=0;k<(src->depth);k++)
            Entry3D(dest,i,j,k) = Entry3D(src,i,j,k);    
    }
}
void write_matrix_to_file(MATRIX_T* a_matrix, const char* filename)
{
    FILE* a_file;
    int i,j,k;
    a_file = fopen(filename,"w");
        
    int h = a_matrix->height;
    int w = a_matrix->width;
    int d = a_matrix->depth;

    if(1==d && 1==w)
    {
        for(i=0;i<h;i++)
            fprintf(a_file, "%f\t", Entry1D(a_matrix,i));
    }
    else if(1==d)
    { // 2D matrix
        for(i=0;i<h;i++){
            for(j=0;j<w;j++){
                fprintf(a_file, "%f\t", Entry2D(a_matrix,i,j));
            }
            fprintf(a_file,"\n");
        }
    }
    else
    {
        for(k=0;k<d;k++)
        for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(a_file, "%f\t", Entry3D(a_matrix,i,j,k));
    }
            
    fclose(a_file);
}
void print_matrix(MATRIX_T* a_matrix)
{
    int i,j,k;
        
    int h = a_matrix->height;
    int w = a_matrix->width;
    int d = a_matrix->depth;

    if(1==d && 1==w)
    {   // 1D matrix
        for(i=0;i<h;i++)
            printf("%f ", Entry1D(a_matrix,i));
    }
    else if(1==d)
    {   // 2D matrix
        for(i=0;i<h;i++)
        {
            for(j=0;j<w;j++)
                printf("%f ", Entry2D(a_matrix,i,j));
            printf("\n");
        }
    }
    else
   {    // 3D matrix 
        for(k=0;k<d;k++)
        {
            for(i=0;i<h;i++)
            {
                for(j=0;j<w;j++)
                    printf("%f ", Entry3D(a_matrix,i,j,k));
                printf("\n");
            }
            printf("\n");
            printf("\n");
        }
    }
    fflush(stdout);
}


/* Median Filter */
static int partition(float *a, int p, int q)
{
    float x = a[p];
    int i = p, j;

    for(j=p+1;j<=q;j++)
    {
        if(a[j]<=x){
            i = i + 1;
            swap(&a[i],&a[j]);
        }
    }
    
    swap(&a[p],&a[i]);
    return i;
}

static int RANDOMIZED_partition(float *a, int p, int q)
{
    float temp = ( ((float) rand()) / RAND_MAX ); // i is a uniform r.v. from [0,1]
    temp = ( (float) (q-p))*temp + (float) p; // now i is a uniform r.v. from [p,q]
    int i = (int) temp;

    swap(&a[p],&a[i]);
    return partition(a, p, q);
}

static float RANDOMIZED_select(float *a, int p, int q, int i)
{
    if(p==q) return a[p];
    int r = RANDOMIZED_partition(a,p,q);
    int k = r-p+1;
    if(i==k) 
        return a[r];
    if(i<k)  
        return RANDOMIZED_select(a,p,r-1,i);
    if(i>k)  
        return RANDOMIZED_select(a,r+1,q,i-k);
    
    return -1; // avoid warning: control reaches end of non-void function;
               // At the end of the function, add a return statement that 
               // returns a suitable return value, even if control never reaches there.
}

// Get the median out of T* a, starting at p, ending at q
float GetMedian(float *a, int p, int q)
{
    assert( p<=q );
    if(p==q)
        return a[p];

    int N = q-p+1;

    if( !(N%2) )
    {
        int m2 = N / 2;
        int m1 = m2 - 1;

        float v1 = RANDOMIZED_select(a,p,q,m1+1);// m1 denotes the index when array is 0-index, 
        float v2 = RANDOMIZED_select(a,p,q,m2+1);// RANDOMIZE_select uses the index when array is 1-index
                                                 // Hence, m1+1 -> RANDOMIZED_select

        return (v1+v2) / 2.0;
    }
    else
    {
        int m = round( (q-p)/2.0 + 1 ); // Since, in c++, array is 0 index, hence the position of the
                                        // median in a CPP array should be (q-p)/2.0+1 instead of (q-p)/2.0
        return RANDOMIZED_select(a,p,q,m);
    }
}


// Interquantile Range algorithm:
// Percentiles are specified using percentages, from 0 to 100. 
// For an n-element vector X, prctile computes percentiles as follows:
// 1. The sorted values in X are taken to be the 100(0.5/N), 100(1.5/N), ..., 100([N-0.5]/N) percentiles.
// 2. Linear interpolation is used to compute percentiles for percent values between 100(0.5/N) and 100([N-0.5]/N).
// 3. The minimum or maximum values in X are assigned to percentiles for percent values outside that range.
float IQR(float *a, int p, int q)
{
    assert( p<=q );// p points to the start
    int N = q-p+1;
    float *b = (float *) malloc( N*sizeof(float) );
    memcpy(b,a,N*sizeof(float));
    RANDOMIZED_quick_sort(b,p,q);

    int i;
    float *prctiles = (float *) malloc( N*sizeof(float) );

    for(i=1;i<=N;i++)
        prctiles[i-1] = 100.0*( (float) i-0.5) / ((float) N);

    // y1, y2 are the indexes corresponding to
    // 25.0% and 75.0%
    float y1 = (25.0 * N) / 100.0 + 0.5 - 1;
    float a1;// (y1,a1) pair.
    float y1_l = ( floor( y1 ) );
    float y1_r = (  ceil( y1 ) );

    if( y1_r==y1_l )
        a1 = b[ (int) y1_l ];
    else{
        float a1_l = b[ (int) y1_l ];
        float a1_r = b[ (int) y1_r ];

        a1 = y1 * (a1_r-a1_l) / (y1_r-y1_l) + (a1_l*y1_r-a1_r*y1_l) / (y1_r-y1_l);
    }
        

    float y2 = (75.0 * N) / 100.0 + 0.5 - 1;
    float a2;
    float y2_l = ( floor( y2 ) );
    float y2_r = (  ceil( y2 ) );

    if( y2_r==y2_l )
        a2 = b[ (int) y2_l ];
    else{
        float a2_l = b[ (int) y2_l ];
        float a2_r = b[ (int) y2_r ];

        a2 = y2 * (a2_r-a2_l) / (y2_r-y2_l) + (a2_l*y2_r-a2_r*y2_l) / (y2_r-y2_l);
    }
     
    free(b);
    return (a2-a1);
}

static void RANDOMIZED_quick_sort(float *a, int p, int q)
{
    int r;

    if(p<q){
        r = RANDOMIZED_partition(a,p,q);
        RANDOMIZED_quick_sort(a,p,r-1);
        RANDOMIZED_quick_sort(a,r+1,q);
    }
}

static void swap(float *a, float *b)
{
    float temp = *a;
            *a = *b;
            *b = temp;
}

float GetMean(float *a, int p, int q)
{
    assert( p<=q );
    int N = q-p+1;
    int i;
    float sum = 0.0;

    for(i=0;i<N;i++)
        sum += a[i];
    
    float mean_value = sum / ( (float)N );

    return mean_value;
}
float GetVar(float *a, int p, int q)
{
    assert( p<=q );
    if(p==q)
        return 0.0;
    
    float mean_value = GetMean(a, p, q);
    int N = q-p+1;
    int i;
    float sum = 0.0;

    for(i=0;i<N;i++)
        sum += (a[i]-mean_value)*(a[i]-mean_value);

    // unbiased estimator
    float variance = sum / ((float)(N-1));

    return variance;
}

// Auto-scale the data matrix x,
// which is a Nxn matrix (N: the 
// number of samples; n: the number
// of attributes).
// Adapted from autoscale func in normfuncs.q
void Auto_Scale(MATRIX_T* x, int robust)
{
    int N = x->height;
    int n = x->width;

    assert( (1<=N) && (1<=n) );
    int i,j;
    for(j=0;j<n;j++)
    {
        float *a = (float *) malloc( N*sizeof(float) );
        for(i=0;i<N;i++)
            a[i] = Entry2D(x,i,j);

        float scl, mn, md, std;
        if( robust ){
            scl = IQR(a,0,N-1) / 1.349;
            md  = GetMedian(a,0,N-1);

            if(scl>0){
                for(i=0;i<N;i++)
                    Entry2D(x,i,j) = ( Entry2D(x,i,j) - md ) / scl;
            }
            
            if(0==scl){
                scl = GetVar(a,0,N-1);// variance
                mn  = GetMean(a,0,N-1);// mean

                if(scl>0){
                    std = sqrt( scl );
                    for(i=0;i<N;i++)
                        Entry2D(x,i,j) = ( Entry2D(x,i,j) - mn ) / std;
                }
            }
        }
        else{
            scl = GetVar(a,0,N-1);// variance
            mn  = GetMean(a,0,N-1);// mean

            if(scl>0){
                std = sqrt( scl );
                for(i=0;i<N;i++)
                    Entry2D(x,i,j) = ( Entry2D(x,i,j) - mn ) / std;
            }
        }

        free(a);
    }
}

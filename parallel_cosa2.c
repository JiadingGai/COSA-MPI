#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
#include <unistd.h>
#include "matrix_util.h"
#define GAI_DEBUG 0

typedef struct{
    int objectID;
    float distance;
} NN_T; // nearest neighbor struct

typedef struct
{
    int p;
    MPI_Comm comm;//entire grid
    MPI_Comm slice_comm_along_z;// free_coords=[1,1,0]
    MPI_Comm slice_comm_along_x;// free_coords=[0,1,1]
    MPI_Comm depth_comm;// free_coords=[0,0,1]
    MPI_Comm row_comm;// free_coords=[1,0,0]
    MPI_Comm col_comm;// free_coords=[0,1,0]
    int pDims[3];//processes layout in 3D grid; y-x-z order (i.e., row-column-depth order)
    int my_coords[3];//x-y-z order, NOT row-column-depth order
    int my_rank;//my rank in comm
    int my_row_rank;
    int my_col_rank;
    int my_slice_comm_along_z_rank;
    int my_slice_comm_along_x_rank;
    int my_depth_comm_rank;
} GRID_INFO_T;


float Calc_Weight_Error(MATRIX_T* w, MATRIX_T* w_old);
void Setup_grid(GRID_INFO_T* grid);

// Distribute_input is a naive scheme for reading and distributing input data;
// In it, root in grid.com reads and dissipates data to all other processes;
// Distribute_input_scheme_2 is a clever scheme for reading and distributing
// input data. In it, each root in grid.row_comm is in charge of reading and
// distributing x_2 to other processes in the same row_comm and each root in
// grid.col_comm is in charge of reading and distributing x_1 to other processes
// in the same col_comm.
void Distribute_input(const char* filename,GRID_INFO_T grid,MATRIX_T* x_1,MATRIX_T* x_2);
void Distribute_input_scheme_2(const char* filename,GRID_INFO_T grid,MATRIX_T* x_1,MATRIX_T* x_2);

// Given N and p, calculate how to divide up n among p processes UNEVENLY,
// and yet in a most possible load-balanced fashion
// INPUT: N and p
// OUTPUT: N_res
int calculate_resource_allocation(int num,int p);

void Build_derived_KNN_type(NN_T* nn_t_ptr, MPI_Datatype* mpi_knn_t_ptr);
void Reshape(NN_T* global_knn, int N_bar, int K, int pDims0);

void rearrange_weight_matrix( MATRIX_T* w );
void rearrange_distance_matrix( MATRIX_T* D );
float max(float a, float b);

/* 
global variables: read process grid dimensions  
from command line.
*/
int px,py,pz;
int N;// # of samples
int n;// # of attributes
int N_res;
int n_res;
int N_bar;
int n_bar;
static int Noit;
static int Niter;
static int K;
static float lambda, alpha, threshold;

int main(int argc, char* argv[])
{

    char* data_file_name;
    if(argc>1)
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
      px = atoi( argv[10] );
      py = atoi( argv[11] );
      pz = atoi( argv[12] );
    }
    else
    {
      data_file_name = "e7.dat";
      N = 8;
      n = 10;
      // by default, process grid
      // consists of only 1 process
      px = 1;
      py = 1;
      pz = 1;
    }

    N_res = calculate_resource_allocation(N,px);//px=py
    n_res = calculate_resource_allocation(n,pz);

    if(0!=N_res){
        N_bar = (N-N_res) / (px-1);
        assert( N == (N_bar*(px-1)+N_res) );
    }
    else{
        N_bar = N / px;
        assert( (N == (N_bar*px)) );
    }

    if(0!=n_res){
        n_bar = (n-n_res) / (pz-1);
        assert( n == (n_bar*(pz-1)+n_res) );
    }
    else{
        n_bar = n / pz;
        assert( (n == (n_bar*pz)) );        
    }

    /* Initialization */
    MPI_Init(&argc,&argv);

    GRID_INFO_T grid;
    Setup_grid(&grid);

    /* if(0==grid.my_rank) */
    /* { */
    /*     //ConvertData_text2binary("e4.txt","e4.dat",4,6); */
    /*     ConvertData_text2binary("c22_ceu.txt","c22_ceu.dat",1000,13000); */
    /* } */

    if(0==grid.my_rank)
    {
      printf("(%d,%d,%d,%d)\n",N_bar,n_bar,N_res,n_res);
    }

    //take timing
    double start, finish;
    MPI_Barrier(grid.comm);
    start = MPI_Wtime();

    int i, j, k, ss, t, v;// s already used.

    MATRIX_T* x_1; 
    MATRIX_T* x_2;

    int N_y_bar, N_x_bar, n_z_bar;

    if( grid.my_coords[0]==(px-1) && 0!=N_res ) {
        N_x_bar    = N_res;
    }
    else{
        N_x_bar    = N_bar;
    }

    if( grid.my_coords[1]==(py-1) && 0!=N_res ) {
        N_y_bar    = N_res;
    }
    else{
        N_y_bar    = N_bar;
    }

    if( grid.my_coords[2]==(pz-1) && 0!=n_res ){
        n_z_bar   = n_res;
    }
    else{
        n_z_bar   = n_bar;
    }
       
    x_1 = Matrix_Allocate(n_z_bar,N_x_bar,1);// x1 is an N_x_bar-by-n_z_bar data matrix
    x_2 = Matrix_Allocate(n_z_bar,N_y_bar,1);// x2 is an N_y_bar-by-n_z_bar data matrix


    /* input IO */
    Distribute_input_scheme_2(data_file_name, grid, x_1, x_2);

#if GAI_DEBUG
    if(grid.my_rank==5)
    {
        printf("%d\n",x_2->height);
        printf("%d\n",x_2->width);
        
        for(i=0;i<x_2->height;i++){
        for(j=0;j<x_2->width;j++){
            printf("%f ",Entry2D(x_2,i,j));
            }
        printf("\n");
        }
    }
#endif    
    //Distribute_input(data_file_name, grid, x_1, x_2);
    
    //int K = (int) sqrt((float) N); // # of KNN; Ref: COSA paper, section 7, last paragraph
    MATRIX_T* w = Matrix_Allocate(N,n,1);// w is an n-by-N matrix
    MATRIX_T* local_w = Matrix_Allocate(N_y_bar,n_z_bar,1);// w is an n-by-N matrix
    MATRIX_T* local_w_old = Matrix_Allocate(N_y_bar,n_z_bar,1);// w_old is an n-by-N matrix 
    MATRIX_T* s = Matrix_Allocate(1,n_z_bar,1);// s is an n-by-1 vector
    MATRIX_T* local_s = Matrix_Allocate(1,n_z_bar,1);
    MATRIX_T* D = Matrix_Allocate(N_x_bar,N_y_bar,1);// D is an N-by-N distance matrix between all pairs
    MATRIX_T* local_D = Matrix_Allocate(N_x_bar,N_y_bar,1);
    MATRIX_T* other_D = Matrix_Allocate(N_y_bar,N_x_bar,1);
    MATRIX_T* collect_D = Matrix_Allocate(N,N,1);
    /* //MATRIX_T* d = Matrix_Allocate(N,N,n);// d is an N-by-N-by-n distance matrix per attribute */
    NN_T* local_knn  = (NN_T*) malloc(N_y_bar*K*sizeof(NN_T));// knn is an N_bar-by-K matrix.
    for(i=0;i<N_y_bar;i++)
    for(j=0;j<K;j++)
        (*((local_knn)+i*K+j)).distance = FLT_MAX;

    NN_T* global_knn = (NN_T*) malloc(N_y_bar*(grid.pDims[0])*K*sizeof(NN_T));// global_knn is an N-by-pK matrix.
    NN_T* knn  = (NN_T*) malloc(N_y_bar*K*sizeof(NN_T));// knn is an N_bar-by-K matrix.
    MATRIX_T* local_S = Matrix_Allocate(N_y_bar,n_z_bar,1);// S is an n-by-N matrix 
    MATRIX_T* S = Matrix_Allocate(N_y_bar,n_z_bar,1);// S is an n-by-N matrix 
    MATRIX_T* local_Sum = Matrix_Allocate(1,N_y_bar,1);// Sum is an N-by-1 vector 
    MATRIX_T* Sum = Matrix_Allocate(1,N_y_bar,1);// Sum is an N-by-1 vector 
    //float lambda = 0.2, alpha = 0.10;// COSA paper, section 6
    float eta = lambda;// COSA paper
    int noit = 0, niter = 0;// outer iterations 
    float local_weight_error = FLT_MAX;// Initialize to a huge number; Average error.
    float weight_error = FLT_MAX;// Initialize to a huge number; Average error.
    //float threshold = 1e-5;// w stablizes when (weight_error < threshold)
    /* char wtcomb[] = "each";// weight combining (30) or (33) */


    // Initialize the weights
    for(k=0;k<n_z_bar;k++)
    for(i=0;i<N_y_bar;i++)
        Entry2D(local_w,k,i) = 1.0 / ((float) n);
                                 //1.0;// weights add up to the number of attributes
    while( noit<Noit /*&& weight_error>threshold*/ )
    {

        if(0==grid.my_rank)
          printf("Outer iteration# %d (weight_error = %f)\n",noit,weight_error);
        
        noit++;
        Copy_Matrix(local_w,local_w_old);
        
        niter = 0;
        while(niter<Niter)
        {
            if(0==grid.my_rank)
                printf("Inner iteration# %d\n", niter);
            niter++;
            // Entering the while loop
            for(k=0;k<n_z_bar;k++)
            {
                Entry1D(local_s,k) = 0.0;
                for(i=0;i<N_y_bar;i++)
                for(j=0;j<N_x_bar;j++)
                {
                    Entry1D(local_s,k) += fabs(Entry2D(x_2,i,k) - Entry2D(x_1,j,k)) / ((float) N*N);
                }
            }

            MPI_Allreduce( &Entry1D(local_s,0),&Entry1D(s,0),n_z_bar,MPI_FLOAT,MPI_SUM,grid.slice_comm_along_z );

            for(i=0;i<N_y_bar;i++)
            for(j=0;j<N_x_bar;j++)
                Entry2D(local_D,i,j) = 0.0;

            for(k=0;k<n_z_bar;k++)
            for(i=0;i<N_y_bar;i++)
            for(j=0;j<N_x_bar;j++)
            {
                //Entry3D(d,i,j,k) = fabs(Entry2D(x,i,k) - Entry2D(x,j,k)) / Entry1D(s,k);
                float dijk = fabs(Entry2D(x_2,i,k) - Entry2D(x_1,j,k)) / Entry1D(s,k);
                Entry2D(local_D,i,j) += Entry2D(local_w,k,i) * exp( -dijk / eta );
            }    

            MPI_Allreduce( &Entry2D(local_D,0,0),&Entry2D(D,0,0),N_y_bar*N_x_bar,MPI_FLOAT,MPI_SUM,grid.depth_comm );

            for(i=0;i<N_y_bar;i++)
            for(j=0;j<N_x_bar;j++)
            {
                Entry2D(D,i,j) = -eta * log( Entry2D(D,i,j) );
            }

            if(grid.my_coords[0]==grid.my_coords[1])
            { // me on the diagonal of the process grid
                for(i=0;i<N_y_bar;i++)
                for(j=0;j<N_x_bar;j++)
                {
                    if(i>=j)
                    {
                        float temp_max = max( Entry2D(D,i,j),Entry2D(D,j,i) );
                        Entry2D(D,i,j) = temp_max;
                        Entry2D(D,j,i) = temp_max;                
                    }
                }
            }
            else
            {
                MPI_Status status;
                if( grid.my_coords[0]>grid.my_coords[1] )
                {
                    int other_side_coords[3];// my symmetric process across the diagonal
                    other_side_coords[0] = grid.my_coords[1];
                    other_side_coords[1] = grid.my_coords[0];
                    other_side_coords[2] = grid.my_coords[2];

                    int other_side_rank;
                    MPI_Cart_rank(grid.comm,other_side_coords,&other_side_rank);
                    MPI_Send(&Entry2D(D,0,0),N_y_bar*N_x_bar,MPI_FLOAT,other_side_rank,1,grid.comm);
                    MPI_Recv(&Entry2D(other_D,0,0),N_y_bar*N_x_bar,MPI_FLOAT,other_side_rank,1,grid.comm,&status);

#if GAI_DEBUG                    
                    printf("%d (%d,%d,%d)-> %d\n",grid.my_rank,grid.my_coords[0],grid.my_coords[1],
                          grid.my_coords[2],other_side_rank);
#endif
                    // flip other_D
                    for(i=0;i<N_y_bar;i++)
                    for(j=0;j<N_x_bar;j++)
                        Entry2D(D,i,j) = Entry2D(other_D,j,i);
                }
                else
                {
                    int other_side_coords[3];
                    other_side_coords[0] = grid.my_coords[1];
                    other_side_coords[1] = grid.my_coords[0];
                    other_side_coords[2] = grid.my_coords[2];

                    int other_side_rank;
                    MPI_Cart_rank(grid.comm,other_side_coords,&other_side_rank);
                    MPI_Recv(&Entry2D(other_D,0,0),N_y_bar*N_x_bar,MPI_FLOAT,other_side_rank,1,grid.comm,&status);

                    for(i=0;i<N_y_bar;i++)
                    for(j=0;j<N_x_bar;j++)
                        Entry2D(D,i,j) = max( Entry2D(D,i,j),Entry2D(other_D,j,i) );
            
                    MPI_Send(&Entry2D(D,0,0),N_y_bar*N_x_bar,MPI_FLOAT,other_side_rank,1,grid.comm);
                }        
            }

            // KNN step -- assume N_bar > K
            for(i=0;i<N_y_bar;i++)
            {
                for(j=0;j<N_x_bar;j++)
                {
                    if(j<K)
                    {
                        (*(local_knn+i*K+j)).objectID = j + grid.my_coords[0]*N_bar;//global data point ID (in row_comm)
                        (*(local_knn+i*K+j)).distance = Entry2D(D,i,j);
                    }
                    else
                    {
                        // Find max distance among current K nearest neighbors
                        int max_pos = -1;
                        float max_distance = -1.0*FLT_MAX;
                
                        for(k=0;k<K;k++)
                        {
                            if( max_distance<(*(local_knn+i*K+k)).distance )
                            {
                                max_pos = k;
                                max_distance = (*(local_knn+i*K+k)).distance;
                            }
                        }
                        // See if the current D[i][j] can replace max_distance
                        // in the current K nearest neighbors.
                        if( Entry2D(D,i,j)<max_distance )
                        {
                            (*(local_knn+i*K+max_pos)).objectID = j + grid.my_coords[0]*N_bar;//global data point ID
                            (*(local_knn+i*K+max_pos)).distance = Entry2D(D,i,j);
                        }    
                    }
                }
            }

            NN_T an_NN_T_Object;// an NN_T object, merely used to feed Build_derived_KNN_type(.)
            MPI_Datatype mpi_knn_t;
            Build_derived_KNN_type( &(an_NN_T_Object),&mpi_knn_t );
            
            MPI_Allgather(local_knn,N_y_bar*K,mpi_knn_t,global_knn,N_y_bar*K,mpi_knn_t,grid.row_comm);

            // Regroup global_knn so that it lays out as
            // if "2D" local_knn stacked up together.
            Reshape(global_knn,N_y_bar,K,grid.pDims[0]);

            for(i=0;i<N_y_bar;i++)
            {
                for(j=0;j<((grid.pDims[0])*K);j++)
                {
                    if(j<K)
                    {
                        (*(knn+i*K+j)).objectID = (*(global_knn+i*(grid.pDims[0])*K+j)).objectID;
                        (*(knn+i*K+j)).distance = (*(global_knn+i*(grid.pDims[0])*K+j)).distance;
                    }
                    else
                    {
                        int max_pos  = -1;
                        float max_distance = -1.0*FLT_MAX;
                        for(k=0;k<K;k++)
                        {
                            if( max_distance<(*(knn+i*K+k)).distance )
                            {
                                max_pos = k;
                                max_distance = (*(knn+i*K+k)).distance;
                            }
                        }
                        if( (*(global_knn+i*(grid.pDims[0])*K+j)).distance<max_distance )
                        {
                            (*(knn+i*K+max_pos)).objectID = (*(global_knn+i*(grid.pDims[0])*K+j)).objectID;
                            (*(knn+i*K+max_pos)).distance = (*(global_knn+i*(grid.pDims[0])*K+j)).distance;
                        }                
                    }
                }
            }

            /* // not necessary since every process computes knn simultaneously */
            /* MPI_Bcast(knn,N_bar*K,mpi_knn_t,0,grid.row_comm); */
#if GAI_DEBUG
            if(0==grid.my_rank)
            {
                printf("****************************************\n");
                printf("sizeof(int) = %d\n",sizeof(int));
                printf("sizeof(float) = %d\n",sizeof(float));
                printf("sizeof(NN_T) = %d\n",sizeof(NN_T));
                printf("sizeof(MPI_Aint) = %d\n",sizeof(MPI_Aint));
                NN_T paul_NN;
                MPI_Aint start_address;
                MPI_Aint address1, address2;
                MPI_Address( &(paul_NN), &start_address);
                MPI_Address( &(paul_NN.objectID), &address1 );
                MPI_Address( &(paul_NN.distance), &address2 );
                printf("address of paul_NN = %d\n",start_address);
                printf("address of paul_NN.objectID = %d\n", address1);
                printf("address of paul_NN.distance = %d\n", address2);
                printf("****************************************\n");

                printf("*********D(i,j)**********\n");
                for(i=0;i<N_bar;i++)
                {
                    for(j=0;j<N_bar;j++)
                        printf("%f   ",Entry2D(D,i,j));                         
                    printf("\n");
                }

                printf("**********local_knn**********\n");
                for(i=0;i<N_bar;i++)
                {
                    for(j=0;j<K;j++)
                        printf("%f(%Ld)    ",(*(local_knn+i*K+j)).distance,(*(local_knn+i*K+j)).objectID);                     
                    printf("\n");
                }

                printf("**********global_knn**********\n");
                for(i=0;i<N_bar;i++)
                {
                    for(j=0;j<((grid.pDims[0])*K);j++)
                        printf("%f(%Ld)    ",(*(global_knn+i*(grid.pDims[0])*K+j)).distance,(*(global_knn+i*(grid.pDims[0])*K+j)).objectID);                printf("\n");
                }

                printf("**********knn**********\n");
                for(i=0;i<N_bar;i++)
                {
                    for(j=0;j<K;j++)
                        printf("%f(%Ld)    ",(*(knn+i*K+j)).distance,(*(knn+i*K+j)).objectID);                     
                    printf("\n");
                }
            }
#endif
            
            for(k=0;k<n_z_bar;k++)
            for(i=0;i<N_y_bar;i++)
            {
                Entry2D(local_S,k,i) = 0.0;
                for(j=0;j<K;j++)
                {
                    // which rank in grid.row_comm has this data point j:
                    int which_rank_in_row_comm = (*(knn+i*K+j)).objectID / N_bar;
                    int local_objectID = ( (*(knn+i*K+j)).objectID-which_rank_in_row_comm*N_bar );
                    if( which_rank_in_row_comm==grid.my_row_rank )
                    {
                        float dijk = fabs(Entry2D(x_2,i,k) - Entry2D(x_1,local_objectID,k)) / Entry1D(s,k);
                        Entry2D(local_S,k,i) += dijk / ((float) K);
                    }
                }
            }

            MPI_Allreduce( &Entry2D(local_S,0,0), &Entry2D(S,0,0), n_z_bar*N_y_bar, MPI_FLOAT, MPI_SUM, grid.row_comm );

            for(k=0;k<n_z_bar;k++)
            for(i=0;i<N_y_bar;i++)
            {
                Entry2D(local_w,k,i) = exp( -Entry2D(S,k,i) / lambda );
            }

            for(i=0;i<N_y_bar;i++)
            {
                Entry1D(local_Sum,i) = 0.0;
                Entry1D(Sum,i) = 0.0;// Set Sum[i] to zero for the MPI_Allreduce 5 lines below.
                for(k=0;k<n_z_bar;k++)
                    Entry1D(local_Sum,i) += Entry2D(local_w,k,i);
            }

            MPI_Allreduce( &Entry1D(local_Sum,0), &Entry1D(Sum,0), N_y_bar, MPI_FLOAT, MPI_SUM, grid.depth_comm );
    
            for(i=0;i<N_y_bar;i++)
            {
                for(k=0;k<n_z_bar;k++)
                {
                    Entry2D(local_w,k,i) /= Entry1D(Sum,i);
                }
            }
            
        }
        eta = eta + alpha*lambda;
        local_weight_error = Calc_Weight_Error(local_w,local_w_old);
        MPI_Allreduce( &local_weight_error, &weight_error, 1, MPI_FLOAT, MPI_MAX, grid.comm );
    }

    /* output IO */
    int x, y, z;
    int local_width, local_height;
    int* recvcounts_w = (int *) malloc( pz*py*sizeof(int) );
    int* displace_w = (int *) malloc( pz*py*sizeof(int) );
    int offset_w = 0;
    int process_index = 0;
    for(y=0;y<py;y++)
    for(z=0;z<pz;z++)
    {
        if( z==(pz-1) && (0!=n_res) )
            local_width = n_res;
        else
            local_width = n_bar;
        
        if( y==(py-1) && (0!=N_res) )
            local_height = N_res;
        else
            local_height = N_bar;

        displace_w[process_index] = offset_w;
        recvcounts_w[process_index] = local_width*local_height;
        offset_w += recvcounts_w[process_index];
        process_index += 1;
    }

    MPI_Gatherv( &Entry2D(local_w,0,0), n_z_bar*N_y_bar, MPI_FLOAT, 
                 &Entry2D(w,0,0), recvcounts_w, displace_w, MPI_FLOAT, 
                 0, grid.slice_comm_along_x);

    int* recvcounts_D = (int *) malloc( px*py*sizeof(int) );
    int* displace_D   = (int *) malloc( px*py*sizeof(int) );
    int offset_D = 0;
    process_index = 0;
    for(x=0;x<px;x++)
    for(y=0;y<py;y++)
    {
        if( x==(px-1) && (0!=N_res) )
            local_width = N_res;
        else
            local_width = N_bar;

        if( y==(py-1) && (0!=N_res) )
            local_height = N_res;
        else
            local_height = N_bar;

        displace_D[process_index] = offset_D;
        recvcounts_D[process_index] = local_width*local_height;
        offset_D += recvcounts_D[process_index];
        process_index += 1;
    }

    MPI_Gatherv( &Entry2D(D,0,0), N_y_bar*N_x_bar, MPI_FLOAT, 
                 &Entry2D(collect_D,0,0), recvcounts_D, displace_D, MPI_FLOAT, 
                 0, grid.slice_comm_along_z);
#if GAI_DEBUG
    if(1==grid.my_rank)
    {
        printf("========= local_w ============\n");
        printf("(%d: %d,%d,%d)\n",grid.my_rank,grid.my_row_rank,grid.my_slice_comm_along_z_rank,grid.my_slice_comm_along_x_rank);
        print_matrix(local_w);
        printf("========= D (local version) ============\n");
        printf("(%d: %d,%d,%d)\n",grid.my_rank,grid.my_row_rank,grid.my_slice_comm_along_z_rank,grid.my_slice_comm_along_x_rank);
        print_matrix(D);
        fflush(stdout);
    }
#endif

    if(0==grid.my_rank)
    {
        /* printf("===========  collect_D  =============\n"); */
        /* print_matrix(collect_D); */
        /* printf("===========  w before re-arrangment =============\n"); */
        /* print_matrix(w); */


        rearrange_weight_matrix( w );
        /* printf("========== w after rearrangement ===========\n"); */
        /* print_matrix(w); */
        char weight_file_name[50];
        //sprintf(weight_file_name,"myweight_par_px=%d_py=%d_pz=%d_noit=%d_niter=%d.txt",px,py,pz,noit,niter);
        sprintf(weight_file_name,"myweight_par.txt");
        write_matrix_to_file(w,weight_file_name);

        rearrange_distance_matrix( collect_D );
        /* printf("========== collect_D after rearrangement ===========\n"); */
        /* print_matrix(collect_D); */
        //write_matrix_to_file(collect_D,"mydist_par.txt");
        char distance_file_name[50];
        sprintf(distance_file_name,"mydist_par_px=%d_py=%d_pz=%d_noit=%d_niter=%d.txt",px,py,pz,noit,niter);


        FILE* dist_file_lower_diagonal = fopen(distance_file_name,"w");
        for(j=0;j<N;j++)
        for(i=0;i<N;i++)
        {
          if(i>j)
            fprintf(dist_file_lower_diagonal, "%f,", Entry2D(collect_D,i,j));
        }
        fclose(dist_file_lower_diagonal);

        FILE* dist_file_full = fopen("mydist_par.txt","w");        
        for(i=0;i<N;i++)
            {
          for(j=0;j<N;j++)
          {
            fprintf(dist_file_full, "%f,",Entry2D(collect_D,i,j));
          }
          fprintf(dist_file_full,"\n");
        }
        fclose(dist_file_full);
    }


    // Finally, take timings
    MPI_Barrier(grid.comm);
    finish = MPI_Wtime();

    if(0==grid.my_rank)
    {
        
        char time_file_name[50]; 
        sprintf(time_file_name,"time_%s_px=%d_py=%d_pz=%d_noit=%d.txt",data_file_name,px,py,pz,noit);
        FILE* time_file = fopen(time_file_name,"w");
        fprintf(time_file,"Elapsed time = %f seconds\n",finish-start);
        fprintf(time_file,"noit = %d \n",noit);
        fprintf(time_file,"weight error = %f \n",weight_error);

        fclose(time_file);
    }
    

    MPI_Finalize();

    Free_Matrix(local_w);
    Free_Matrix(local_w_old);
    Free_Matrix(w);
    Free_Matrix(s);
    Free_Matrix(local_s);
    Free_Matrix(D);
    Free_Matrix(other_D);
    Free_Matrix(local_D);
    free(knn);
    free(local_knn);
    free(global_knn);
    Free_Matrix(local_S);
    Free_Matrix(S);
    Free_Matrix(local_Sum);
    Free_Matrix(Sum);
    Free_Matrix(x_1);
    Free_Matrix(x_2);
    return 1;
}
float Calc_Weight_Error(MATRIX_T* w, MATRIX_T* w_old)
{
    float weight_error = 0.0;
    int i,j;
    
    for(i=0;i<(w->height);i++)
    for(j=0;j<(w->width );j++)
        weight_error += fabs( Entry2D(w,i,j)-Entry2D(w_old,i,j) );
    
    weight_error = weight_error / (float)(w->height*w->width);
    return weight_error;
}
void Setup_grid(GRID_INFO_T* grid)
{
    int wrap_around[3] = {0,0,0};
    int free_coords[3];// x-y-z major order

    MPI_Comm_size(MPI_COMM_WORLD,&(grid->p));
    grid->pDims[0] = py;// # of rows
    grid->pDims[1] = px;// # of columns
    grid->pDims[2] = pz;// # of depth

    /*Create 3D grid topology on all processes; Get my rank and grid coords */
    MPI_Cart_create(MPI_COMM_WORLD,3,grid->pDims,wrap_around,1,&(grid->comm));
    MPI_Comm_rank(grid->comm,&(grid->my_rank));
    MPI_Cart_coords(grid->comm,grid->my_rank,3,grid->my_coords);
    
    /*Set up per slice communicator along z axis: free_coords=(1,1,0)*/
    free_coords[0] = 1;
    free_coords[1] = 1;
    free_coords[2] = 0;    
    MPI_Cart_sub(grid->comm,free_coords,&(grid->slice_comm_along_z));
    MPI_Comm_rank(grid->slice_comm_along_z, &(grid->my_slice_comm_along_z_rank));

    /*Set up per slice communicator along x axis: free_coords=(0,1,1)*/
    free_coords[0] = 0;
    free_coords[1] = 1;
    free_coords[2] = 1;    
    MPI_Cart_sub(grid->comm,free_coords,&(grid->slice_comm_along_x));
    MPI_Comm_rank(grid->slice_comm_along_x,&(grid->my_slice_comm_along_x_rank));

    /*Set up depth communicator: free_coords=(0,0,1)*/
    free_coords[0] = 0;
    free_coords[1] = 0;
    free_coords[2] = 1;    
    MPI_Cart_sub(grid->comm,free_coords,&(grid->depth_comm));
    MPI_Comm_rank(grid->depth_comm,&(grid->my_depth_comm_rank));

    /*Set up row communicator: free_coords=(1,0,0)*/
    free_coords[0] = 1;
    free_coords[1] = 0;
    free_coords[2] = 0;    
    MPI_Cart_sub(grid->comm,free_coords,&(grid->row_comm));
    MPI_Comm_rank(grid->row_comm,&(grid->my_row_rank));

    /*Set up column communicator: free_coords=(0,1,0)*/
    free_coords[0] = 0;
    free_coords[1] = 1;
    free_coords[2] = 0;    
    MPI_Cart_sub(grid->comm,free_coords,&(grid->col_comm));
    MPI_Comm_rank(grid->col_comm,&(grid->my_col_rank));
}
float max(float a, float b)
{
    float bigger = a>b?a:b;
    return bigger;
}
void Build_derived_KNN_type(NN_T* nn_t_ptr, MPI_Datatype* mpi_knn_t_ptr)
{
    /*The number of elements in each "block" of the new type*/
    int block_lengths[2];
    block_lengths[0] = block_lengths[1] = 1;

    /*Displacement of each element from start of new type*/
    MPI_Aint displacements[2];

    /*MPI types of the elements.*/
    MPI_Datatype typelist[2];

    /*Use for calculating displacements*/
    MPI_Aint start_address;
    MPI_Aint address;

    /* Build a derived datatype consisting of */
    /* one int and one float */
    typelist[0] = MPI_LONG_LONG;
    typelist[1] = MPI_FLOAT;

    /* First element, objectID, is a displacement 0 */
    MPI_Address(nn_t_ptr,&start_address);
    MPI_Address(&(nn_t_ptr->objectID), &address);
    displacements[0] = address - start_address;

    /* Find address of distance and displacement from objectID*/
    MPI_Address(&(nn_t_ptr->distance), &address);
    displacements[1] = address - start_address;

    /* Build A_NN_T_Object-the derived MPI KNN type */
    MPI_Type_struct(2,block_lengths,displacements,typelist, mpi_knn_t_ptr);

    /* Commit it -- tell system we will be using it for communication */
    MPI_Type_commit(mpi_knn_t_ptr);
}

void Reshape(NN_T* global_knn, int N_y_bar, int K, int pDims0)
{
    NN_T* temp_knn = (NN_T*) malloc(N_y_bar*K*pDims0*sizeof(NN_T));
    memcpy(temp_knn,global_knn,N_y_bar*K*pDims0*sizeof(NN_T));

    int i,j;
    for(i=0;i<N_y_bar;i++)
    for(j=0;j<(pDims0*K);j++)
    {
        int x,y,z;
        z = j / K;// which local block
        y = i;
        x = j % K;
        *(global_knn+i*(pDims0*K)+j) = *(temp_knn+z*(N_y_bar*K)+y*K+x);
    }

    free(temp_knn);
}
void Distribute_input(const char* filename, GRID_INFO_T grid, MATRIX_T* x_1, MATRIX_T* x_2)
{
    // repeat the computation of N,n,N_bar,n_bar
    int N_bar = x_1->height;
    int n_bar = x_1->width;
    int N = N_bar * (grid.pDims[0]);
    int n = n_bar * (grid.pDims[2]);        
    int coords[3];// y-x-z major order (i.e., row-column-depth order)
    MATRIX_T* temp_x_1 = Matrix_Allocate(n_bar,N_bar,1);// x is an N-by-n data matrix
    MATRIX_T* temp_x_2 = Matrix_Allocate(n_bar,N_bar,1);// x is an N-by-n data matrix
    int i,j,k,ss,t,v;
    
    if(0==(grid.my_rank))
    {
        /* FILE* e4_file = fopen("e0.dat","wb"); */
        
        /* for(i=1;i<=N*n;i++) */
        /* { */
        /*     float temp = (float) ( 99.0*rand() / RAND_MAX ); */
        /*     //printf("---write--to-file: %f \n",temp); */
        /*     fwrite(&temp,sizeof(float),1,e4_file); */
        /* } */
        /* fclose(e4_file); */

        FILE* input_file = fopen(filename,"rb");            

#if GAI_DEBUG
        printf("--------- Print input data ---------\n");
        for(i=0;i<N;i++)
        {
            for(j=0;j<n;j++)
            {
                float temp;
                fread(&temp,sizeof(float),1,input_file);
                printf("%f ",temp);
            }
            printf("\n");
        }
#endif
        int dest_process;
        
        for(i=0;i<(grid.pDims[0]);i++)
        for(j=0;j<(grid.pDims[1]);j++)
        for(k=0;k<(grid.pDims[2]);k++)
        {
            coords[0] = i;//row
            coords[1] = j;//column
            coords[2] = k;//depth
            
            MPI_Cart_rank(grid.comm,coords,&dest_process);
#if GAI_DEBUG
            printf("~~dest_process~~%d~~~~~~~~\n",dest_process);
#endif
            for(t=0;t<N_bar;t++)
            {
                for(v=0;v<n_bar;v++)
                {
                    fseek(input_file,((j*N_bar+t)*n+k*n_bar+v)*sizeof(float),SEEK_SET);
                    fread(&Entry2D(temp_x_1,t,v),sizeof(float),1,input_file);
#if GAI_DEBUG
                    printf("%f ",Entry2D(temp_x_1,t,v));
#endif
                }
#if GAI_DEBUG
                printf("\n");
#endif
            }
#if GAI_DEBUG
            printf("----------------\n");
#endif

            for(ss=0;ss<N_bar;ss++)
            {
                for(v=0;v<n_bar;v++)
                {
                        fseek(input_file,((i*N_bar+ss)*n+k*n_bar+v)*sizeof(float),SEEK_SET);
                        fread(&Entry2D(temp_x_2,ss,v),sizeof(float),1,input_file);
#if GAI_DEBUG
                        printf("%f ",Entry2D(temp_x_2,ss,v));
#endif
                }
#if GAI_DEBUG
                printf("\n");
#endif
            }
#if GAI_DEBUG
            printf("=============\n");
#endif
            if(0==dest_process)
            {
                Copy_Matrix(temp_x_1,x_1);// I am reading my data => x_1 = temp_x_1;
                Copy_Matrix(temp_x_2,x_2);// I am reading my data => x_2 = temp_x_2;
            }
            else
            {
                MPI_Send(&Entry2D(temp_x_1,0,0),n_bar*N_bar,MPI_FLOAT,dest_process,1, grid.comm);
                MPI_Send(&Entry2D(temp_x_2,0,0),n_bar*N_bar,MPI_FLOAT,dest_process,1, grid.comm);
            }
        }
        fclose(input_file);
    }
    else
    {
        MPI_Status status;
        MPI_Recv(&Entry2D(x_1,0,0),n_bar*N_bar,MPI_FLOAT,0,1,grid.comm,&status);
        MPI_Recv(&Entry2D(x_2,0,0),n_bar*N_bar,MPI_FLOAT,0,1,grid.comm,&status);
    }

    Free_Matrix(temp_x_1);
    Free_Matrix(temp_x_2);
}
void Distribute_input_scheme_2(const char* filename,GRID_INFO_T grid,MATRIX_T* x_1,MATRIX_T* x_2)
{
    
    int i,j,k,ss,t,v;

    if(0==grid.my_col_rank)
    {
        FILE* input_file = fopen(filename,"rb");

        // Note: grid.my_coords is in x-y-z order
        i = grid.my_coords[1];
        j = grid.my_coords[0];
        k = grid.my_coords[2];

        for(t=0;t<(x_1->height);t++)
        {
            for(v=0;v<(x_1->width);v++)
            {
                fseek(input_file,((j*N_bar+t)*n+k*n_bar+v)*sizeof(float),SEEK_SET);
                fread(&Entry2D(x_1,t,v),sizeof(float),1,input_file);
                //printf("%f ",Entry2D(temp_x_1,t,v));
            }
            //printf("\n");
        }

        fclose(input_file);
    }
    MPI_Bcast( &Entry2D(x_1,0,0),(x_1->height)*(x_1->width),MPI_FLOAT,0,grid.col_comm );

    //------------------------------------------------------------------//

    if(0==grid.my_row_rank)
    {
        FILE* input_file = fopen(filename,"rb");

        // Note: grid.my_coords is in x-y-z order
        i = grid.my_coords[1];
        j = grid.my_coords[0];
        k = grid.my_coords[2];

        for(ss=0;ss<(x_2->height);ss++)
        {
            for(v=0;v<(x_2->width);v++)
            {
                fseek(input_file,((i*N_bar+ss)*n+k*n_bar+v)*sizeof(float),SEEK_SET);
                fread(&Entry2D(x_2,ss,v),sizeof(float),1,input_file);
                //printf("%f ",Entry2D(x_2,ss,v));
            }
            //printf("\n");
        }

        fclose(input_file);
    }
    MPI_Bcast( &Entry2D(x_2,0,0),(x_2->height)*(x_2->width),MPI_FLOAT,0,grid.row_comm );
}
int calculate_resource_allocation(int num,int p)
{
    if(1==p) return 0;
    if( !(num%p) ) return 0;

    float x = ( (float) num ) / p;
    int num_res = (int) x;

    while( (num-num_res)%(p-1) ) {
        num_res = num_res-1;
    }
    return num_res;
}
void rearrange_weight_matrix( MATRIX_T* w )
{
    MATRIX_T * temp_w = Matrix_Allocate(N,n,1);
    int i, j;
    for(i=0;i<n;i++)
    for(j=0;j<N;j++)
            Entry2D(temp_w,i,j) = Entry2D(w,i,j);
    
    for(i=0;i<n;i++)
    {
        for(j=0;j<N;j++)
        {
            int grid_i = i / n_bar;// (i,j) position in the process 
            int grid_j = j / N_bar;// grid of position (grid_i,grid_j).
            
            int my_local_width, my_local_height;

            if( grid_j==(py-1) && (0!=N_res) )
                my_local_width = N_res;
            else
                my_local_width = N_bar;
            
            if( grid_i==(pz-1) && (0!=n_res) )
                my_local_height = n_res;
            else
                my_local_height = n_bar;

            int local_z = (i-grid_i*n_bar) % my_local_height;
            int local_y = (j-grid_j*N_bar) % my_local_width;

            int NumofElement_b4_me = grid_j*N_bar*n + grid_i*n_bar*my_local_width;

            Entry2D(w,i,j) = *(temp_w->entries+NumofElement_b4_me+local_z*my_local_width+local_y);
            Entry2D(w,i,j) *= (float) n;
        }
    }
    Free_Matrix(temp_w);
}
void rearrange_distance_matrix( MATRIX_T* D )
{
    MATRIX_T * temp_D = Matrix_Allocate(N,N,1);
    int i, j;
    for(i=0;i<N;i++)
    for(j=0;j<N;j++)
            Entry2D(temp_D,i,j) = Entry2D(D,i,j);
    
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            int grid_i = i / N_bar;// (i,j) position in the process 
            int grid_j = j / N_bar;// grid of position (grid_i,grid_j).
            
            int my_local_width, my_local_height;

            if( grid_j==(px-1) && (0!=N_res) )
                my_local_width = N_res;
            else
                my_local_width = N_bar;
            
            if( grid_i==(py-1) && (0!=N_res) )
                my_local_height = N_res;
            else
                my_local_height = N_bar;

            int local_y = (i-grid_i*N_bar) % my_local_height;
            int local_x = (j-grid_j*N_bar) % my_local_width;

            int NumofElement_b4_me = grid_j*N_bar*N + grid_i*N_bar*my_local_width;

            Entry2D(D,i,j) = *(temp_D->entries+NumofElement_b4_me+local_y*my_local_width+local_x);
        }
    }
    Free_Matrix(temp_D);
}

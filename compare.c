#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

float compare_2_matrix_file(  char* f1_name, char* f2_name);
int main(int argc, char* argv[])
{

    //compare_2_matrix_file("myweight.txt","myweight_par.txt");
    compare_2_matrix_file("mydist_par_px=2_py=2_pz=3_noit=20_niter=5.txt","mydist_par_px=1_py=1_pz=5_noit=20_niter=5.txt");
    _Exit(0);

    int NumOfTests = 0;
    int noit = 20;
    int niter = 5;

    srand( time(0) );

    while( 1 ){

        sleep(1);

      int N = (int) (7.0+80.0 * (rand() / (RAND_MAX + 1.0)));
      int n = (int) (17.0+136.0 * (rand() / (RAND_MAX + 1.0)));

        int px = (int) (1.0+4.0 * (rand() / (RAND_MAX + 1.0)));
        while( ((int)((float) N / px)) < ((int) sqrt( (float) N)) )
            px = (int) (1.0+4.0 * (rand() / (RAND_MAX + 1.0)));
        int py = px;


        int pz = (int) (1.0+4.0 * (rand() / (RAND_MAX + 1.0)));
        while( ((int)((float) n / pz)) < 2 )
            pz = (int) (1.0+4.0 * (rand() / (RAND_MAX + 1.0)));


        int p  = px*py*pz;

        printf("(N=%d,n=%d,px=%d,py=%d,pz=%d)\n",N,n,px,py,pz);
        fflush(stdout);
        char command1[50];
        sprintf(command1,"make clean;make;./cosa2 e8.dat %d %d %d %d",N,n,noit,niter);
        printf("%s\n",command1);
        system(command1);

        char command2[50];
        sprintf(command2,"mpirun -np %d parallel_cosa2 e8.dat %d %d %d %d %d %d %d",p,N,n,px,py,pz,noit,niter);
        printf("%s\n",command2);
        system(command2);

        //char weight_file_name[50];
        //sprintf(weight_file_name,"myweight_par_px=%d_py=%d_pz=%d_noit=%d_niter=%d.txt",px,py,pz,noit,niter);

        //char distance_file_name[50];
        //sprintf(distance_file_name,"mydist_par_px=%d_py=%d_pz=%d_noit=%d_niter=%d.txt",px,py,pz,noit,niter);

        float error_weight   = compare_2_matrix_file( "myweight.txt","myweight_par.txt" );
        float error_distance = compare_2_matrix_file( "mydist.txt"  ,"mydist_par.txt" );

        if( (error_weight>=0.01) || (error_distance>=0.01) || (error_weight!=error_weight) || (error_distance!=error_distance) )
        {
            printf("ALERT: Testing case failed: (N=%d,n=%d,px=%d,py=%d,pz=%d)\n",N,n,px,py,pz);
            printf("%d tests run ... Before failure\n",NumOfTests);            
            fflush(stdout);
            _Exit(0);
        }
    }
    printf("%d tests run ... All passed\n",NumOfTests);
    return 1;
}

float compare_2_matrix_file(  char* f1_name, char* f2_name)
{

  FILE* f1 = fopen(f1_name,"r");
  FILE* f2 = fopen(f2_name,"r");
  float n1, n2;
  float error = 0.0;
  int num_of_element_read = 0;
  while ( (EOF!=fscanf(f1,"%f,",&n1)) && (EOF!=fscanf(f2,"%f,",&n2)) )
  {
    //printf("(n1=%f,n2=%f)\n",n1,n2);
    error += (n1-n2)*(n1-n2);
    num_of_element_read += 1;
  }
    
  error /= num_of_element_read;
  printf("The number of elements read is %d\n",num_of_element_read);
  printf("The difference between [%s,%s] is %f\n",f1_name,f2_name,error);

  fclose(f1);
  fclose(f2);

  return error;
}

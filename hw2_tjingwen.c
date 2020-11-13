
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#define PI 3.1415926535
#define max_rows = 10000
#define max_cols = 10000
double A[10000][10000];
double A1[10000][10000];

double f(double x){
    double y;
    int i;
    y = x;
    for (i = 1; i < 11; i ++){
        y = y + sin(x*i)/pow(2.0,i);
    }
    return y;
}

double ztoa(double x){
    double temp = 100;
    double temp2 = -100;
    if (x <= 100){
        temp = x;
    }
    if (temp > -100){
        temp2  = temp;
    }
    return temp2;
}

int main(int argc, char* argv[])
{
    int my_id, root_process, num_procs, ierr;
    int m, n, i, j, start_row, end_row, avg_rows, num_rows;
    double z, soln[2], answer[2];
    double partial_sum=0;
    double partial_sum2=0;

    /* Let process 0 be the root process. */

    root_process = 0;
    
    /* Now replicate this process to create parallel processes. */

    ierr = MPI_Init(&argc, &argv);
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    /* Find out MY process ID, and how many processes were started. */

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = 0;
    if(my_id == root_process) {
        time_start = MPI_Wtime();
        
     /* I must be the root process, so I will query the user
        to determine how many interpolation intervals to use. */
    }
    avg_rows = floor(m / num_procs);
    start_row = my_id*avg_rows;
    end_row  = (my_id+1)*avg_rows - 1;
    if (my_id == num_procs - 1){
        end_row = m-1;
    }
    num_rows = end_row - start_row + 1;
    for (i = start_row; i < end_row+1;i++){
        for (j = 0;j<n;j++){
            A[i][j] = i*sin(i) + j*cos(j) + sqrt(i+j);
        }
    }
    
    
    
    for (int k=0;k<10;k++){
        
        if (num_procs > 1){
            if (my_id == root_process){
                MPI_Send(&A[end_row][0], n, MPI_DOUBLE, 1, k, MPI_COMM_WORLD);
                MPI_Recv(&A[end_row+1][0], n, MPI_DOUBLE, 1, k, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            
            else if (my_id == num_procs-1){
                MPI_Send(&A[start_row][0], n, MPI_DOUBLE, num_procs-2, k, MPI_COMM_WORLD);
                MPI_Recv(&A[start_row-1][0], n, MPI_DOUBLE, num_procs-2, k, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else{
                MPI_Send(&A[start_row][0], n, MPI_DOUBLE, my_id-1, k, MPI_COMM_WORLD);
                MPI_Send(&A[end_row][0], n, MPI_DOUBLE, my_id+1, k, MPI_COMM_WORLD);
                MPI_Recv(&A[start_row-1][0], n, MPI_DOUBLE, my_id-1, k, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&A[end_row+1][0], n, MPI_DOUBLE, my_id+1, k, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        
         
        /*Then we calculate the matrix*/
        for (i = start_row; i < end_row+1; i++){
            for (j = 0; j < n; j++){
                if(i==0 | i==(m-1) | j==0 | j==(n-1)){
                    A1[i][j] = A[i][j];
                }
                else{
                    z = (f(A[i-1][j]) + f(A[i+1][j]) + f(A[i][j-1]) + f(A[i][j+1]) + f(A[i][j]))/5;
                    A1[i][j] = ztoa(z);
                }
            }
        }
        
        for (i = start_row; i < end_row+1;i++){
            for (j = 0;j<n;j++){
                A[i][j] = A1[i][j];
            }
        }
    }
    
    
    
    for (i = start_row; i < end_row+1; i++){
        for (j = 0; j < n; j++){
            partial_sum += A1[i][j];
            partial_sum2 += pow(A1[i][j], 2);
        }
    }
    soln[0] = partial_sum;
    soln[1] = partial_sum2;
    ierr = MPI_Reduce(&soln, &answer, 2, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);

    /* and, if I am the root process, print the result. */

    if(my_id == root_process) {
        double duration = MPI_Wtime() - time_start;
        printf("time:%fs %lf %lf\n", (float)duration, answer[0], answer[1]);
    }

    /* Close down this processes. */

    ierr = MPI_Finalize();
}

#include <iostream>
#include <deque>
#include <omp.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f.h"
using namespace std;
struct Node{
    double l,r,fl,fr,M;
};

int main(int argc, char* argv[])
{
    double M,M1;
    double middle,Mc,Md,fm;
    
    deque<Node> in; //beginning and ending position of interval
    Node po,pu;
    po.l = a;
    po.r = b;
    po.fl = f(a);
    po.fr = f(b);
    po.M = (po.fl + po.fr + s*(b-a))/2;
    M = fmax(po.fl,po.fr);
    in.push_back(po);
    


    
    int numth = omp_get_max_threads();//number of existing threads
    vector <int> work(numth,0);//vector work or not 1: work; 0:idle
    
    
    omp_lock_t lck; //lock for intervals deque
    omp_init_lock (&lck);
    
    // time start
    double start_time = omp_get_wtime();
    int i;
    for (i = 0;i<numth;i++){
        work[i] = 0;
    }
    
    //then use DFS to search for max value
    // Beginning of parallel region
    #pragma omp parallel default(shared) private(M1,Mc,Md,middle,fm,po,pu)
    {
        int tid = omp_get_thread_num();
        int i;
        bool check = true;

        while(true){
            check = true;
            omp_set_lock (&lck);
            if (!in.empty()){
                check = false;
                work[tid] = 1;

                po = in.front();
                in.pop_front();
            }
            else {
                work[tid] = 0;
                for (i=0;i<numth;i++){
                    if (work[i] == 1){
                        check = false;
                    }
                }
            }
            omp_unset_lock (&lck);
            if (check == true){
                break;
            }
            if (work[tid] == 0){
                continue;
            }
            
            #pragma omp atomc read
            M1 = M;
            if (po.M < M1 + epsilon){
                continue;
            }
            
            middle = (po.l + po.r)/2;
            fm = f(middle);
            Mc = 0.5*(po.fl + fm + s*((middle-po.l)/2));
            Md = 0.5*(fm + po.fr + s*((po.r-middle)/2));
            
            #pragma omp atomic write
            M = fmax(fm,M);
            
            omp_set_lock (&lck);
            if (Mc >= M + epsilon){
                pu.l = po.l;
                pu.r = middle;
                pu.fl = po.fl;
                pu.fr = fm;
                pu.M = Mc;
                in.push_front(pu);
            }
            if (Md >= M + epsilon){
                pu.l = middle;
                pu.r = po.r;
                pu.fl = fm;
                pu.fr = po.fr;
                pu.M = Md;
                in.push_front(pu);
            }
            omp_unset_lock (&lck);

            }
        }// Ending of parallel region

    omp_destroy_lock(&lck);
    cout<<"Maximum value: "<<M<<endl;
    double elapse_time = omp_get_wtime()-start_time;
    cout<<"Time:"<<elapse_time<<endl;
    return 0;
}



//double findmax(double (*f) (),double a,double b, double epsi, double s){
    
//}

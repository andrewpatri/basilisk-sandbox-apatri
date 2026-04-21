// this file is a test for the declaration in a point of q_source

// #include "embed.h"


scalar q[];


int maxlevel = 7; int minlevel = 2;
double D0 = 3.e-2;
double tend = 10.; 

int main(){

    // q changes iteratively with time
    q[top] = 10. * pow(t, 1.5);

    // square domain 
    L0 = D0 *3;

    origin(-L0/2, 0);

    run();

}

//define rectangle
event init(i=0){
    q[top] = 10.;

    //solid    
}

event output( t =+ 1) {

    fprintf (stderr, "%g/n, t");
}

event stop ( t = tend)
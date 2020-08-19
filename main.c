#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>
#include    <string.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_blas.h> //not used will remove later 
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h> 
#include    <suitesparse/umfpack.h> //not used
// HEADER FILES I CREATED

#include "read_header.h"
#include "data_types.h"

// constants

int main()
{ 
    //this is to read data from the file 
    COMB_STRUCT All_data  =  Reader_fn();
    
    //creating a variable of suitable dimention blocks to get the initial values as initial values has many things nside it
    INITIAL * Initial_state = Gen_data(All_data);
    
    //getting Y_BUS
    Y_STRUCT Y = Y_BUS(All_data);
    
    //SIMULTANEOUS SOLVER FUNCTION
    double ** jack = jacobian(All_data,Initial_state,Y,10,0.1);
    
    printf("\n done \n");    
    
    return 0;    
}




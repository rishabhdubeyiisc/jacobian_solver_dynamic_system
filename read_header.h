#ifndef READ_HEADER
#define READ_HEADER

#include "data_types.h"


;
COMB_STRUCT Reader_fn();

INITIAL * Gen_data(COMB_STRUCT);

void write(COMB_STRUCT);

V_VECTOR Voltage_vector_maker(COMB_STRUCT);

Y_STRUCT Y_BUS(COMB_STRUCT);

I_VECTOR Current_vector_maker(V_VECTOR ,Y_STRUCT ,COMB_STRUCT );

gsl_complex * rotator(gsl_complex ,gsl_complex ,gsl_complex ,double );

gsl_complex*rotate_all_vectors_initial_state(COMB_STRUCT,INITIAL*);

gsl_complex ** complex_matrix_invertor(Y_STRUCT ,COMB_STRUCT);
//MYCONTAINER Manager(char * filename); // jus the file is required and pass the file pointer to the other fn
double partial_f_wrt_x(double ,double );

double Power_at_ith(double ,double ,double ,double ,double ,double);

double ** jacobian(COMB_STRUCT ,INITIAL*,Y_STRUCT ,int ,double);


#endif

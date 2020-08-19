#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h>

// HEADER FILES I CREATED

#include    "read_header.h"
#include    "data_types.h"
//torque solver
double Torque(double Power,double omega_operating)

{
    double Torque;
    Torque=(Power/omega_operating);
    return Torque;
};
//partial solver
double partial_f_wrt_x(double y,double x)

{   
    double par_diff=1; 
    return par_diff; 
};
//power eqn solver
double Power_at_ith(double ED_das_das_ith,double EQ_das_das_ith,double VD_ith,double VQ_ith,double G_ii,double B_ii)

{   
    double term1=((ED_das_das_ith*ED_das_das_ith)-(ED_das_das_ith*VD_ith))*G_ii;
    double term2=((EQ_das_das_ith*EQ_das_das_ith)-(EQ_das_das_ith*VQ_ith))*G_ii;
    double term3=((EQ_das_das_ith*VD_ith)-(ED_das_das_ith*VQ_ith))*B_ii;;
    return term1+term2+term3;
};
//MATRIX INVERSE
double * matrix_inverse(double*A)

{   double*mat_inverse;
    return mat_inverse;
}; 
 //SIMULTANEOUS SOLVER FUNCTION

STATES * Simultaneous_solver(COMB_STRUCT All_data,INITIAL*Initial_state,Y_STRUCT Y,int T,double time_step )
 
{   
    int x=0;//change this later
    double D;
    double omega_initial;
    double omega_base;
    double del_t_by_2=(time_step*(0.5));
    int number_of_generator=All_data.constants[0].NumberofGens;
    int number_of_buses=All_data.constants[0].LOAD_FLOW_number;
        
    double Efd0[number_of_generator];
    for (int i = 0; i < number_of_generator; i++)
    {
       Efd0[i]=gsl_complex_abs(Initial_state[i].Efd_0);
    }

    
     
    //Y_BUS_SPLITTER
    double G[number_of_buses][number_of_buses];
    double B[number_of_buses][number_of_buses];

    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; i < number_of_buses; j++)
        {
            G[i][j]=GSL_REAL(Y.MAT[i][j]);
            B[i][j]=GSL_IMAG(Y.MAT[i][j]);
        }
        
    }
    
    //MAKE A 6 variable array so can take all state variable vector of all gens
    //Ed_das,Eq_das,Ed_das_das,Eq_das_das,delta,slip
    STATES state_array[number_of_generator];

    //frame changers for ED and EQ to calculate 
    double ED_das_das_ith_gen[number_of_generator];
    double EQ_das_das_ith_gen[number_of_generator];
    for (int i = 0; i < number_of_generator; i++)
    {
        ED_das_das_ith_gen[i]= ((state_array[i].Ed_das_das)*cos(state_array[i].delta))+((state_array[i].Eq_das_das)*sin(state_array[i].delta));
        EQ_das_das_ith_gen[i]= ((state_array[i].Ed_das_das)*(-1)*sin(state_array[i].delta))+((state_array[i].Eq_das_das)*cos(state_array[i].delta));
    }
           
    
   //double VD,VQ  //DEFINE A NETWORK VECTOR OF BUSES 
    NW_STATES network_var_array[number_of_buses];
   

    //defining VD,VQ

    //POWER EQN AND MODIFY THE EQUATIONS OF STATES
    //done
        
    double mech_power_array[number_of_generator];
    for (int i = 0; i < number_of_generator; i++)
    {   double t_m[number_of_generator];
        t_m[i]=GSL_REAL(Initial_state[i].Tm_0);
        mech_power_array[i]=((t_m[i])/omega_initial);
    }
    //done
    
    //MAKE A TORQUE SOLVER
    //done
        
    //GET A INVERSE FN
    //modify this

    //Define differential of states_vector
    DIFF_STATES_VECTOR f_of_state_vector[number_of_generator];
    //defining some variable
    double tdo1=(All_data.Generator_ps[x].Tdo1);
    double tqo1=(All_data.Generator_ps[x].Tqo1);
    double tdo2=(All_data.Generator_ps[x].Tdo2);
    double tqo2=(All_data.Generator_ps[x].Tqo2);
    double xd2_mul_xq2_of_gen=(All_data.Generator_ps[x].Xd2)*(All_data.Generator_ps[x].Xq2);

    double xd0_of_gen=All_data.Generator_ps[x].Xd0;
    double xq0_of_gen=All_data.Generator_ps[x].Xq0;
    double xd1_of_gen=All_data.Generator_ps[x].Xd1;
    double xq1_of_gen=All_data.Generator_ps[x].Xq1;
    double xd2_of_gen=All_data.Generator_ps[x].Xd2;
    double xq2_of_gen=All_data.Generator_ps[x].Xq2;
    double H=All_data.Generator_ps[x].H;
        
    //frame transform 
        
    double Vd[number_of_generator];
    double Vq[number_of_generator];
    for (int i = 0; i < number_of_generator; i++)
    {
        Vd[x]= network_var_array[x].VD * cos(state_array[x].delta)-(network_var_array[x].VQ * sin(state_array[x].delta));
        Vd[x]= network_var_array[x].VD * sin(state_array[x].delta)+(network_var_array[x].VQ * cos(state_array[x].delta));
    }
       
    //define all f_of_state_vectors these will be 6 for 1 generator 
    f_of_state_vector[x].f_of_Ed_das=(-1/tqo1)*(state_array[x].Ed_das)+((xq1_of_gen-xq0_of_gen)/tqo1)*((state_array[x].Ed_das_das-Vq[x])/xq2_of_gen);
    f_of_state_vector[x].f_of_Eq_das=(-1/tdo1)*(state_array[x].Eq_das)+(1/tdo1)*(gsl_complex_abs(Initial_state[x].Efd_0))+((xd0_of_gen-xd1_of_gen)/tdo1)*((Vq[x]-state_array[x].Eq_das_das)/xd2_of_gen);
    f_of_state_vector[x].f_of_Ed_das_das=(-1/tqo2)*(state_array[x].Ed_das)+(-1/tqo2)*(state_array[x].Ed_das_das)+((xq2_of_gen-xq0_of_gen)/tqo2)*((state_array[x].Ed_das_das-Vd[x])/xq2_of_gen);
    f_of_state_vector[x].f_of_Eq_das_das=(-1/tdo2)*(state_array[x].Eq_das)+(-1/tdo2)*(state_array[x].Eq_das_das)+((xd0_of_gen-xd2_of_gen)/tdo2)*((Vq[x]-state_array[x].Eq_das_das)/xd2_of_gen);
    f_of_state_vector[x].f_of_delta=omega_base*(state_array[x].slip);
    f_of_state_vector[x].f_of_slip=((mech_power_array[x]-Power_at_ith(ED_das_das_ith_gen[x],EQ_das_das_ith_gen[x],network_var_array[x].VD,network_var_array[x].VQ,G[x][x],B[x][x]))/(2*All_data.Generator_ps[x].H));

    //define the 2 state vectors IgQi and IgDi depends on VQ_kth , VD_kth and transfer impedance in all simple I added up to all

    //DEFINE A MISMATCH VECTOR OF NETWORK-----g
    //these are just required later
    double IgQ_ith[number_of_buses];
    double IgD_ith[number_of_buses];
    
    NW_STATES G_mismatch[number_of_buses];
    
    for (int i = 0; i < number_of_buses ; i++)
    {   double delta=state_array[i].delta;
        double E_ith=sqrt(pow(state_array[i].Ed_das_das,2)+pow(state_array[i].Eq_das_das,2));
        IgQ_ith[i]=E_ith*(cos(delta)*G[i][i]-B[i][i]*sin(delta)); 
        IgD_ith[i]=E_ith*(cos(delta)*B[i][i]+G[i][i]*sin(delta));  
        double summer_Q=0;
        double summer_D=0;
            for (int j = 0; i < number_of_buses; i++)
            {
                summer_Q = (G[i][j])*(network_var_array[i].VD)-(B[i][j])*(network_var_array[i].VQ);
                summer_D = (G[i][j])*(network_var_array[i].VQ)+(B[i][j])*(network_var_array[i].VD);
                summer_Q+=summer_Q;
                summer_Q+=summer_Q;
            }
        
        G_mismatch[i].VQ =summer_Q-IgQ_ith[i];   
        G_mismatch[i].VD =summer_D-IgD_ith[i];
    }
    
    //define next states of generatos   
    //define mismatch of states 
  
    //making Jacobian 
    //GEN BLOCK
   
    //vref changer ther should be
    
        
    //NOW MAKE A LOOP

    //DONE

    //PRINT FUNCTION
};

double ** jacobian(STATES  * state_array ,COMB_STRUCT All_data,INITIAL*Initial_state,Y_STRUCT Y,int T,double time_step )

{   
    int x=0;//change this later
    
    double omega_base;
    double del_t_by_2=(time_step*(0.5));
    int number_of_generator=All_data.constants[0].NumberofGens;
    int number_of_buses=All_data.constants[0].LOAD_FLOW_number;
    
    
    double ED_das_das_ith_gen[number_of_generator];
    double EQ_das_das_ith_gen[number_of_generator];
    for (int i = 0; i < number_of_generator; i++)
    {
        ED_das_das_ith_gen[i]= ((state_array[i].Ed_das_das)*cos(state_array[i].delta))+((state_array[i].Eq_das_das)*sin(state_array[i].delta));
        EQ_das_das_ith_gen[i]= ((state_array[i].Ed_das_das)*(-1)*sin(state_array[i].delta))+((state_array[i].Eq_das_das)*cos(state_array[i].delta));
    }

    // J_GEN_BLOCK is a multidimentional array 
    //set at zeo first
    //ten put values generator wise
    //these are 3 diagonal blocks
    double J_GEN_BLOCK[number_of_generator][6][6];
    //set this to zero
    for (int i = 0; i < number_of_generator; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                J_GEN_BLOCK[i][j][k]=0;
            }
            
        }
        
    }
    
    double D[number_of_generator]; //define this damping ratio

    for (int i = 0; i < number_of_generator; i++)
    {
        double tdo1=(All_data.Generator_ps[i].Tdo1);
        double tqo1=(All_data.Generator_ps[i].Tqo1);
        double tdo2=(All_data.Generator_ps[i].Tdo2);
        double tqo2=(All_data.Generator_ps[i].Tqo2);
        double H=All_data.Generator_ps[i].H;

        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {   
                J_GEN_BLOCK[i][0][0]=(-1)-((1/tqo1)*(del_t_by_2));
                J_GEN_BLOCK[i][1][1]=(-1)-((1/tdo1)*(del_t_by_2));
                J_GEN_BLOCK[i][2][0]=(-1)-((1/tqo2)*(del_t_by_2));
                J_GEN_BLOCK[i][2][2]=(-1)-((1/tqo2)*(del_t_by_2));
                J_GEN_BLOCK[i][3][1]=(-1)-((1/tdo2)*(del_t_by_2));
                J_GEN_BLOCK[i][3][3]=(-1)-((1/tdo2)*(del_t_by_2));
                J_GEN_BLOCK[i][4][4]=(-1);
                J_GEN_BLOCK[i][4][5]=(del_t_by_2)*(omega_base);
                J_GEN_BLOCK[i][5][4]=0.5;
                J_GEN_BLOCK[i][5][5]=(-1)-((-D[x]/(2*H))*(del_t_by_2));   
            }
            
        }
        
    }
    
    
    //splitiing Y_BUS
    double G_n[number_of_buses][number_of_buses];
    double B_n[number_of_buses][number_of_buses];

    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            G_n[i][j]=GSL_REAL(Y.MAT[i][j]);
            B_n[i][j]=GSL_IMAG(Y.MAT[i][j]);
        }
        
    }
   

    Y_2X2 Y_MAT[number_of_buses][number_of_buses];
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_MAT[i][j].Y_2X2[0][0]=G_n[i][j];
            Y_MAT[i][j].Y_2X2[0][1]=(-B_n[i][j]);   
            Y_MAT[i][j].Y_2X2[1][0]=B_n[i][j];   
            Y_MAT[i][j].Y_2X2[1][1]=G_n[i][j];      
        }
        
    }
    
    double GEN_BUS_BLOCK[6*number_of_generator][2*number_of_buses];
    for (int i = 0; i < 6*number_of_generator; i++)
    {
        for (int j = 0; j < 2*number_of_buses; j++)
        {
            GEN_BUS_BLOCK[i][j]=0;
        }
        
    }
    
    
    for (int i = 5; i < 6*number_of_generator; i=i+6)
    {
        for (int j = 0; j < number_of_generator; j++)
        {
            double H=All_data.Generator_ps[j].H;
            //for VQ then VD 
            GEN_BUS_BLOCK[i][j]=(del_t_by_2)*(0.5/H)*(-1)*(ED_das_das_ith_gen[j]*B_n[j][j]+EQ_das_das_ith_gen[j]*G_n[j][j]);
            GEN_BUS_BLOCK[i][j+1]=(del_t_by_2)*(0.5/H)*(EQ_das_das_ith_gen[j]*B_n[j][j]-ED_das_das_ith_gen[j]*G_n[j][j]);
            //see the equations correctly 
        }
        
    }
    
    
   printf("done");

};


INITIAL * Gen_data(COMB_STRUCT NetData)
    
    
    {   
        
        CONSTANT   * constant    = &NetData.constants[0];
        
        int Number_of_gens       = constant[0].NumberofGens;

        INITIAL *  Gen_initial;
        Gen_initial=malloc(Number_of_gens*sizeof(INITIAL));
        
        gsl_complex I_0      [Number_of_gens];
        gsl_complex Eq_0     [Number_of_gens];
        gsl_complex id_0     [Number_of_gens];
        gsl_complex iq_0     [Number_of_gens];
        gsl_complex vd_0     [Number_of_gens];
        gsl_complex vq_0     [Number_of_gens];
        gsl_complex Efd_0    [Number_of_gens];
        gsl_complex Eq_dash_0[Number_of_gens];
        gsl_complex Ed_dash_0[Number_of_gens];
        gsl_complex Te_0     [Number_of_gens];
        gsl_complex Tm_0     [Number_of_gens];
        
        
        double delta_0_arr   [Number_of_gens];
        double Ed_das_das    [Number_of_gens]; 
        double Eq_das_das    [Number_of_gens];   
  


        MCD        *  Gen_array      = &NetData.Generator_ps[0];
        LOAD_FLOW  *  Load_flow_arr  = &NetData.Load_flow_ps[0];
        
        for (int x = 0; x < Number_of_gens  ; x++, Gen_array++, Load_flow_arr++)
        {
            double theta             = Load_flow_arr->theta; 
            double angle             = ((Load_flow_arr->theta)*PI)/180;
            gsl_complex z            = gsl_complex_rect(0, angle); 
            
            gsl_complex e_j_theta    = gsl_complex_exp(z);
            
            gsl_complex Vt_0         = gsl_complex_rect(Load_flow_arr->V,0.0) ;
            double      mag_Vt_0     = gsl_complex_abs(Vt_0);

            gsl_complex Pg           = gsl_complex_rect(Load_flow_arr->Pg,0.0);
            gsl_complex Qg           = gsl_complex_rect(0.0,Load_flow_arr->Qg);
            gsl_complex Xq0          = gsl_complex_rect(0.0,Gen_array->Xq0);
            
            
            I_0[x] =gsl_complex_div(gsl_complex_sub(Pg,Qg),gsl_complex_mul(Vt_0,e_j_theta));
            Gen_initial[x].I_0=I_0[x];
            
            double mag_I_0    = gsl_complex_abs(I_0[x]);
            double phi_0      = gsl_complex_arg(I_0[x]);
            
            Eq_0[x]= gsl_complex_add((gsl_complex_mul(I_0[x],Xq0)),(gsl_complex_mul(Vt_0,e_j_theta)));
            Gen_initial[x].Eq_0=Eq_0[x];
            
            double delta_0    = gsl_complex_arg(Eq_0[x]);
            delta_0_arr[x]=delta_0;
            Gen_initial[x].delta_0=delta_0_arr[x];

            id_0[x]     = gsl_complex_rect(( -mag_I_0   ) * ( sin (delta_0-phi_0) ),0.0) ;
            iq_0[x]     = gsl_complex_rect((  mag_I_0   ) * ( cos (delta_0-phi_0) ),0.0) ;
            
            vd_0[x]     = gsl_complex_rect(( -mag_Vt_0  ) * ( sin (delta_0-theta) ),0.0) ;
            vq_0[x]     = gsl_complex_rect((  mag_Vt_0  ) * ( cos (delta_0-theta) ),0.0) ;

            Gen_initial[x].id_0=id_0[x];
            Gen_initial[x].iq_0=iq_0[x];
            Gen_initial[x].vd_0=vd_0[x];
            Gen_initial[x].vq_0=vq_0[x];


            Efd_0[x]    = gsl_complex_rect
                            (gsl_complex_abs( Eq_0[x]) - ( ( (Gen_array->Xd0  - Gen_array->Xq0 ) ) * gsl_complex_abs(id_0[x] ) )
                            ,0.0);

            Eq_dash_0[x]= gsl_complex_rect
                            (gsl_complex_abs(Efd_0[x]) - ( ( (Gen_array->Xd0) - (Gen_array->Xd1) ) * gsl_complex_abs(id_0[x] ) )
                            ,0.0);
            
            Ed_dash_0[x]= gsl_complex_rect
                            (
                                (0 - ( (  (Gen_array->Xq0) - (Gen_array->Xq1) )                    * gsl_complex_abs(iq_0[x] ) ) )
                            ,0.0);
            
            gsl_complex temp=gsl_complex_rect((Gen_array->Xd1)-(Gen_array->Xq1),0.0);

            Te_0[x] =gsl_complex_rect
                        ((( (( gsl_complex_abs(Eq_dash_0[x]) ) * ( gsl_complex_abs(iq_0[x]) ) )
                        + (( gsl_complex_abs(Ed_dash_0[x]) ) * (gsl_complex_abs(id_0[x]) ) ) ) 
                        + ((gsl_complex_abs(temp))*(gsl_complex_abs(id_0[x]))*(gsl_complex_abs(iq_0[x])))),
                        0.0);
            
            
            Tm_0[x]     = Te_0[x];

            Ed_das_das[x]=GSL_REAL(Ed_dash_0[x])-GSL_REAL(iq_0[x])*(Gen_array->Xq1-Gen_array->Xq2);
            
            Eq_das_das[x]=GSL_REAL(Eq_dash_0[x])+GSL_REAL(id_0[x])*(Gen_array->Xd1-Gen_array->Xd2);
            
            
            Gen_initial[x].Efd_0    =Efd_0[x];
            Gen_initial[x].Eq_dash_0=Eq_dash_0[x];
            Gen_initial[x].Ed_dash_0=Ed_dash_0[x];
            Gen_initial[x].Te_0     =Te_0[x];
            Gen_initial[x].Tm_0     =Tm_0[x];
            
            Gen_initial[x].Ed_das_das=Ed_das_das[x];
            Gen_initial[x].Eq_das_das=Eq_das_das[x];
            Gen_initial[x].delta_0   =delta_0_arr[x];


        }
        
        return Gen_initial ;
    };
  





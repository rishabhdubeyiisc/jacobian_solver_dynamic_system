#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES I CREATED

#include    "read_header.h"
#include    "data_types.h"


COMB_STRUCT Reader_fn()
{   
    COMB_STRUCT     Array_of_pointers;
    COMB_STRUCT   * pointer_to_arr_blocks = & Array_of_pointers;
    
    CONSTANT     * constant_arr;
    MCD          * generator_arr;
    EXCT         * exciter_arr;
    LOAD_FLOW    * Load_flow_arr;
    TX_PARA      * trans_line_para_arr;
   
    
    FILE         * get_data;
    get_data     = fopen("RISHABH_DATA.txt","r");
    
    if(get_data == NULL) 
    {  
        printf("DATA NOT EXECUTED\n");  
        exit(1); 
    };
    

    constant_arr        = malloc(sizeof(CONSTANT));
    
    fscanf(get_data,"%d",&constant_arr[0].NumberofGens);
    fscanf(get_data,"%d",&constant_arr[0].Number_of_lines);
    fscanf(get_data,"%d",&constant_arr[0].LOAD_FLOW_number);
    
    int NumberofGens        = constant_arr[0].NumberofGens;
    int Number_of_lines     = constant_arr[0].Number_of_lines;
    int LOAD_FLOW_number    = constant_arr[0].LOAD_FLOW_number;
    
    generator_arr       = malloc(NumberofGens     * sizeof(MCD));
    exciter_arr         = malloc(NumberofGens     * sizeof(EXCT));
    Load_flow_arr       = malloc(LOAD_FLOW_number * sizeof(LOAD_FLOW));
    trans_line_para_arr = malloc(Number_of_lines  * sizeof(TX_PARA));
    
  
    
    for ( int i = 0; i < NumberofGens ; i++)
    {   
        fscanf(get_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" ,
                &generator_arr[i].H,
                &generator_arr[i].Xd0,
                &generator_arr[i].Xd1,
                &generator_arr[i].Xd2,
                &generator_arr[i].Xq0,
                &generator_arr[i].Xq1,
                &generator_arr[i].Xq2,
                &generator_arr[i].Tdo1,
                &generator_arr[i].Tdo2,
                &generator_arr[i].Tqo1,
                &generator_arr[i].Tqo2);
    
    }
    
    for (int i = 0; i < NumberofGens; i++)
    {
        fscanf( get_data,"%lf %lf %lf %lf %lf %lf",
                &exciter_arr[i].ka,
                &exciter_arr[i].ta,
                &exciter_arr[i].ke,
                &exciter_arr[i].te,
                &exciter_arr[i].kf,
                &exciter_arr[i].tf 
                                    );
    }
    
    for (int i = 0; i < Number_of_lines; i++)
    {
         fscanf(get_data,"%d %d %lf %lf %lf",
                &trans_line_para_arr[i].bus1,
                &trans_line_para_arr[i].bus2,
                &trans_line_para_arr[i].R,
                &trans_line_para_arr[i].X,
                &trans_line_para_arr[i].B
                                            );
    }
    
    for (int i = 0; i < LOAD_FLOW_number; i++)
    {
        fscanf( get_data,"%d %d %lf %lf %lf %lf %lf %lf",
                &Load_flow_arr[i].bus_number,
                &Load_flow_arr[i].bus_type,
                &Load_flow_arr[i].V,
                &Load_flow_arr[i].theta,
                &Load_flow_arr[i].Pg,
                &Load_flow_arr[i].Qg,
                &Load_flow_arr[i].Pl,
                &Load_flow_arr[i].Ql
                                    );
    }
    
   
    pointer_to_arr_blocks->constants        = constant_arr         ;
    pointer_to_arr_blocks->Generator_ps     = generator_arr        ;
    pointer_to_arr_blocks->exciter_ps       = exciter_arr          ;
    pointer_to_arr_blocks->Load_flow_ps     = Load_flow_arr        ;
    pointer_to_arr_blocks->trans_line_para  = trans_line_para_arr  ;

    //new here

   
    fclose(get_data)        ;
    
    return Array_of_pointers;
};

Y_STRUCT Y_BUS(COMB_STRUCT data)
{   FILE       * y_data;
    y_data     = fopen("Y_BUS.txt","r");
    
    if(y_data == NULL) 
    {  
        printf("DATA in Y BUS text NOT EXECUTED\n");  
        exit(1); 
    };
    int Bus_number=data.constants[0].LOAD_FLOW_number;
    int rows=Bus_number;
    int column=Bus_number;
    /*
    gsl_complex ** MAT;
    MAT=malloc(rows*sizeof(gsl_complex));

    for (int i = 0; i < column; i++)
    {
        MAT[i]=malloc(column*sizeof(gsl_complex*));
    }
    for (int i = 0; i < rows ; i++)
    {
      for (int j = 0; j < column; j++)
      {
          fscanf(y_data,"%lf %lf",&MAT[i][j].dat[0],&MAT[i][j].dat[1]);
      }
      
    }*/
    Y_STRUCT MATRIX;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < column; j++)
        {
            fscanf(y_data,"%lf %lf",&MATRIX.MAT[i][j].dat[0],&MATRIX.MAT[i][j].dat[1]);
        }
        
    }
    
    
    return MATRIX;
    
   /* gsl_complex  MAT[rows][column]; 
    
    //MAT = malloc(rows * column * sizeof(gsl_complex));

    //fscanf(get_data,"%lf %lf",&z.dat[0],&z.dat[1]);
    
    for (int x = 0; x <rows; x++)
    {
        for (int y = 0; y < 2*column; y++)
        {
            fscanf(y_data,"%lf %lf",&MAT[x][y].dat[0],&MAT[x][y].dat[1]);
        }
        
    }*/

};


void write(COMB_STRUCT NetData)
{

   // CONSTANT     * constant_arr          = &NetData.constants[0];
    MCD          * generator_arr         = &NetData.Generator_ps[0];
    EXCT         * exciter_arr           = &NetData.exciter_ps[0];
   // LOAD_FLOW    * Load_flow_arr         = &NetData.Load_flow_ps[0];
   // TX_PARA      * trans_line_para_arr   = &NetData.trans_line_para[0];

        int a = NetData.constants[0].NumberofGens;
        int b = NetData.constants[0].Number_of_lines;
        int c = NetData.constants[0].LOAD_FLOW_number;
        printf("/n NUMBER OF GENERATOR : %d \n ",a);
        printf("/n NUMBER OF LINES     : %d \n ",b);
        printf("/n NUMBER OF BUSES     : %d \n ",c);

    for (int i = 0; i < a ; i++ )
    {
        printf("Generator %d: H%d = %lf\n",i+1,i+1,generator_arr[i].H);
        printf("Generator %d: Xd0%d = %lf\n",i+1,i+1,generator_arr[i].Xd0);
        printf("Generator %d: Xd1%d = %lf\n",i+1,i+1,generator_arr[i].Xd1);
        printf("Generator %d: Xd2%d = %lf\n",i+1,i+1,generator_arr[i].Xd2);
        printf("Generator %d: Xq0%d = %lf\n",i+1,i+1,generator_arr[i].Xq0);
        printf("Generator %d: Xq1%d = %lf\n",i+1,i+1,generator_arr[i].Xq1);
        printf("Generator %d: Xq2%d = %lf\n",i+1,i+1,generator_arr[i].Xq2);
        printf("Generator %d: Tdo1%d = %lf\n",i+1,i+1,generator_arr[i].Tdo1);
        printf("Generator %d: Tdo2%d = %lf\n",i+1,i+1,generator_arr[i].Tdo2);
        printf("Generator %d: Tqo1%d = %lf\n",i+1,i+1,generator_arr[i].Tqo1);
        printf("Generator %d: Tqo2%d = %lf\n",i+1,i+1,generator_arr[i].Tqo2);
    }

    for (int i = 0; i < a; i++)
    {
       printf("Exciter %d: Ka%d = %lf\n",i+1,i+1,exciter_arr[i].ka);
       printf("Exciter %d: Ke%d = %lf\n",i+1,i+1,exciter_arr[i].ke);
       printf("Exciter %d: Kf%d = %lf\n",i+1,i+1,exciter_arr[i].kf);
       printf("Exciter %d: Ta%d = %lf\n",i+1,i+1,exciter_arr[i].ta);
       printf("Exciter %d: Te%d = %lf\n",i+1,i+1,exciter_arr[i].te);
       printf("Exciter %d: Tf%d = %lf\n",i+1,i+1,exciter_arr[i].tf);
    }
   
    
   /* 
    for (int i = 0; i < b; i++)
    {
        printf("Transmission line Between BUS (%d) and BUS(%d) Resistance  %lf Reactance  %lf Suceptance %lf \n ",
                trans_line_para_arr[i].bus1,
                trans_line_para_arr[i].bus2,
                trans_line_para_arr[i].R,
                trans_line_para_arr[i].X,
                trans_line_para_arr[i].B
                );

    }
    for (int i = 0; i < c; i++)
    {
        printf("Load flow data is as follow for BUS NUMBER %d BUS TYPE %d\n
                Pg=%lf  Qg=%lf  Pl=%lf  Ql=%lf ",
                Load_flow_arr[i].bus_number,
                Load_flow_arr[i].bus_type,
                Load_flow_arr[i].Pg,
                Load_flow_arr[i].Qg,
                Load_flow_arr[i].Pl,
                Load_flow_arr[i].Ql
                );
    }
    */
};


/*

EXCT* read_exciter_data(int NumberofGens)
{ 
    static EXCT array_var[NumberofGens];
    FILE * reader_variable;
    reader_variable=fopen("data_folder/Exciter_Data.txt","r");
    if (reader_variable == NULL)
        {
            printf("EXITER DATA NOT READABLE");
            exit(1);
        };
    for ( int i = 0; i < NumberofGens ; i++)
    {   
        fscanf(reader_variable,"%lf %lf %lf %lf %lf %lf " ,
                &array_var[i].ka,
                 &array_var[i].ta,
                  &array_var[i].ke,
                   &array_var[i].te,
                    &array_var[i].kf,
                     &array_var[i].tf);

    }
    fclose(reader_variable);
    return array_var;
};
TX_PARA * read_transmission_line(int Number_of_lines)
{   //defining array_var named TX_PARA of variables as structure of TX_PARA_all
    static TX_PARA TX_PARA_array[Number_of_lines];
    FILE *get_data;
    get_data = fopen("data_folder/Transmission_Lline_Data.txt","r");
    if(get_data == NULL) 
    {
        printf("FILE  GENERATOR DATA NOT EXECUTED\n");  
        exit(1);
    };
    for ( int i = 0; i < Number_of_lines ; i++)
    {   
        fscanf( get_data,
                "%d %d %lf %lf %lf",
                &TX_PARA_array[i].bus1,
                &TX_PARA_array[i].bus2,
                &TX_PARA_array[i].R,
                &TX_PARA_array[i].X,
                &TX_PARA_array[i].B);
    
    }
    
    fclose(get_data);
    return TX_PARA_array;

};
LOAD_FLOW* read_LOAD_FLOW_data(int LOAD_FLOW_number)
{   
    FILE *variable_to_file;
    variable_to_file = fopen("data_folder/LOAD_FLOW_Data.txt","r");
    if (variable_to_file == NULL){
        printf("LOAD FLOW NON EXECUTABLE\n");
        exit(1);
    };
    //int temp;
    //fscanf(variable_to_file,"%d",&temp);
    static LOAD_FLOW array_var[LOAD_FLOW_number];
    for (int i = 0; i < LOAD_FLOW_number; i++)
    {
        fscanf( variable_to_file,"%lf %lf %lf %lf %lf %lf %lf",
                &array_var[i].bus_type,
                &array_var[i].V,
                &array_var[i].theta,
                &array_var[i].Pg,
                &array_var[i].Qg,
                &array_var[i].Pl,
                &array_var[i].Ql);
    }
    fclose(variable_to_file);
    return array_var;
};

*/


/*CONSTANT * read_constants()
{
    FILE * filer;
    filer = fopen("RISHABH_DATA.txt","r");
    if(filer == NULL) 
    {  
        printf("CANNOT READ FILE\n"); 
        exit(1); 
    };
    
    CONSTANT * constant_array;
    constant_array=malloc(sizeof(CONSTANT));
    
    fscanf(filer,"%d %d %d ",
                            &constant_array[0].NumberofGens     ,
                            &constant_array[0].Number_of_lines  ,
                            &constant_array[0].LOAD_FLOW_number )
                                                                    ;
    
    return constant_array   ;
};*/
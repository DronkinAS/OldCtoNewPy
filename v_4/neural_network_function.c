#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "neural_network_function.h"

double* read_matrix(char file_name[], int N, int M) 
{
    int i,j;
    double *A;
    FILE *f;
    
    f = fopen(file_name,"r");
    
    A = (double*)malloc(N * M * sizeof(double));   
    
    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            fscanf(f,"%lf",(A + i*M + j));
        }
    }
    
    fclose(f);
    
    return A;
}

void write_matrix(char file_name[], double* A, int N, int M) 
{
    int i,j;
    FILE *f;
    
    f = fopen(file_name,"w");   
    
    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            fprintf(f,"%lf ",*(A + i*M + j));
        }
        fprintf(f,"\n");
    }
    
    fclose(f);
} 

double* straight(double* A, double* B, double *b, int N, int M)
{
    int i, j;
    double *C, s = 0.0;
    
    C = (double*)malloc(N * sizeof(double));
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
        *(C + i) = 0.0;
    }

    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            *(C + i) += (*(A + j)) * (*(B + i*M + j));
        }
        s = (*(C + i)) + (*(b + i));
        *(C + i) = 2/(1 + exp(-1 * s)) - 1;
    }
    return C;
}

double* nnout(double* A, double* B, double *b, int N, int M)
{
    int i, j;
    double *C, s;
    
    C = (double*)malloc(N * sizeof(double));
    s = 0.0;
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
        *(C + i) = 0.0;
    }

    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            *(C + i) += (*(A + j)) * (*(B + i*M + j));
        }
        s = (*(C + i)) + (*(b + i));
        *(C + i) = 2/(1 + exp(-1 * s)) - 1;
    }
    return C;
}

double* delta_out(double out, double error)
{
    int i;
    double *x;
    
    x = (double*)malloc(sizeof(double));
    
    if (( fabs(out) > 0.99 ) && (fabs(error) >= 0.0001))
    
    {
        out = out * 0.99;
    }

    *x = 0.5*error*(1 - pow(out,2));

    return x;
}

double* back(double* W, double* delta,
             double *out,
             int N, int M)
{
    int i,j;
    double *x, s;
    
    x = (double*)malloc(M * sizeof(double));
    
    for ( i = 0 ; i <= M-1 ; i++)
    {
        *(x + i) = 0.0;
    }
    
    for ( i = 0 ; i <= M-1 ; i++)
    {
        s = 0.0;
        
        for ( j = 0 ; j <= N-1 ; j++)
        {
            s += (*(delta + j)) * (*(W + i + j*M));
        }
        
        if (( fabs(*(out + i)) > 0.95 ) && ( fabs(s) >= 0.0001))
        {
            *(out + i) = *(out + i) * 0.95;
        }
        
        *(x + i) = 0.5*s*(1 - pow(*(out + i),2));
    }
        
    return x;
}

void correction_w(double E, double A, double* w, double* dw, double* d, double* out, int N, int M)
{
    int i, j;
    double grad_1, grad_2, grad;

    for ( i = 0 ; i <= N-1 ; i++)
    {
        for ( j = 0 ; j <= M-1 ; j++)
        {
            grad_1 = *(out + j);
            grad_2 = *(d + i);
            grad = grad_1 * grad_2;
            
            *(dw + i*M + j) = (E * grad) + (A * (*(dw + i*M + j)));
            *(w + i*M + j) += *(dw + i*M + j);
        }
    }
}

void correction_b(double E, double A, double* b, double* db, double* d, int N)
{
    int i;
    double grad;    
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
            grad = (*(d + i));
            *(db + i) = E*grad;
            *(db + i) += A * (*(db + i));
            *(b + i) += *(db + i);
    }
}
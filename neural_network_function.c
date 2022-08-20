#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "neural_network_function.h"

//Считывание матрицы из файла
double* read_matrix(char file_name[], int N, int M) 
{
    int i,j;
    double *A;
    FILE *f;
    
    //Открытие файла
    f = fopen(file_name,"r");
    
    //Выделение памяти под матрицу
    A = (double*)malloc(N * M * sizeof(double));   
    
    //Заполнение мытрицы данными из файла
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

//Занесение матрицы в файл
void write_matrix(char file_name[], double* A, int N, int M) 
{
    int i,j;
    FILE *f;
    
    //Открытие файла
    f = fopen(file_name,"w");   
    
    //Запись в файл данных из матрицы
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

//Функция прямого распространения нейронной сети
double* straight(double* A, double* B, double *b, int N, int M)
{
    int i, j;
    double *C, s = 0.0;
    
    //Выделение памяти под выходные значения слоя
    C = (double*)malloc(N * sizeof(double));
    
    //Обнуление выходных значений
    for ( i = 0 ; i <= N-1 ; i++)
    {
        *(C + i) = 0.0;
    }
    
    //Расчет выходов нейронов
    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            *(C + i) += (*(A + j)) * (*(B + i*M + j));
        }
        //Суммарное возмущение нейрона
        s = (*(C + i)) + (*(b + i));
        //Выход после прохождения функции активации
        *(C + i) = 2/(1 + exp(-0.1 * s)) - 1;
    }
    return C;
}

//Расчет выходного нейрона
double* nnout(double* A, double* B, double *b, int N, int M)
{
    int i, j;
    double *C, s;
    
    //Выделение памяти под выход нейрона
    C = (double*)malloc(N * sizeof(double));
    s = 0.0;
    
    //Обнуление выходных значений
    for ( i = 0 ; i <= N-1 ; i++)
    {
        *(C + i) = 0.0;
    }
    
    //Расчет выхода нейрона
    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            *(C + i) += (*(A + j)) * (*(B + i*M + j));
        }
        //Суммарное возмущение нейрона
        s = (*(C + i)) + (*(b + i));
        //Выход после прохождения активационной функции
        *(C + i) = 2/(1 + exp(-0.005 * s)) - 1;
    }
    return C;
}

//Расчет дельты выходного нейрона
double* delta_out(double out, double error)
{
    int i;
    double *x;
    
    //Выделение памяти
    x = (double*)malloc(sizeof(double));
    
    //Искусственный запрет на достижение предельного значения
    if (( fabs(out) > 0.99 ) && (fabs(error) >= 0.0001))
    
    {
        out = out * 0.99;
    }
    
    //Расчет дельты
    *x = 0.0025*error*(1 - pow(out,2));

    return x;
}

//Функция обратного расспространения нейронной сети
double* back(double* W, double* delta,
             double *out,
             int N, int M)
{
    int i,j;
    double *x, s;
    
    //Выделение памяти под дельты нейронов
    x = (double*)malloc(M * sizeof(double));
    
    //Зануление дельт
    for ( i = 0 ; i <= M-1 ; i++)
    {
        *(x + i) = 0.0;
    }
    
    //Расчет дельты для каждого нейрона слоя
    for ( i = 0 ; i <= M-1 ; i++)
    {
        s = 0.0;
        
        for ( j = 0 ; j <= N-1 ; j++)
        {
            s += (*(delta + j)) * (*(W + i + j*M));
        }
        
        //Искусственный запрет на достижение предельного значения
        if (( fabs(*(out + i)) > 0.99 ) && ( fabs(s) >= 0.0001))
        {
            *(out + i) = *(out + i) * 0.99;
        }
        
        *(x + i) = 0.05*s*(1 - pow(*(out + i),2));
    }
        
    return x;
}

//Коррекция весов нейронов
void correction_w(double E, double A, double* w, double* dw, double* d, double* out, int N, int M)
{
    int i, j;
    double grad_1, grad_2, grad;

    for ( i = 0 ; i <= N-1 ; i++)
    {
        for ( j = 0 ; j <= M-1 ; j++)
        {
            //Расчет градиента на основе дельты
            grad_1 = *(out + j);
            grad_2 = *(d + i);
            grad = grad_1 * grad_2;
            
            //Расчет текущей величины коррекции
            *(dw + i*M + j) = (E * grad) + (A * (*(dw + i*M + j)));
            
            //Изменение веса
            *(w + i*M + j) += *(dw + i*M + j);
        }
    }
}

//Коррекция значений смещений нейронов
void correction_b(double E, double A, double* b, double* db, double* d, int N)
{
    int i;
    double grad;    
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
        //Расчет градиента на основе дельты
        grad = (*(d + i));
        
        //Расчет величины коррекции
        *(db + i) = E*grad;
        *(db + i) += A * (*(db + i));
        
        //Изменение значения смещения
        *(b + i) += *(db + i);
    }
}
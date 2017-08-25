//
//  alocacao.h
//  Util
//
//  Created by Felipe Ferreira on 30/11/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _ALOCACAO_H_
#define _ALOCACAO_H_

#include <stdlib.h>
#include <stdio.h>

int *aloca_vetor_int(int N);
int **aloca_matriz_int(int N, int K);
double *aloca_vetor_double(int N);
double **aloca_matrizd(int N, int K);
double **aloca_matriz_double(int N, int K);
double ***aloca_matriz3d(int n, int N, int K);

void desaloca_matriz_int(int **pt, int N);
void desaloca_matrizd(double **pt, int N);
void desaloca_matriz3d(double ***pt, int n, int N);

int *aloca_vetor_int(int N){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro inteiro de dimens„o N.           *
     *******************************************************************/
    int *x,i;
    
    x = (int *)malloc(N*sizeof(int));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<N;i++)
        x[i] = 0;
    
    return (x);
}

double *aloca_vetor_double(int N){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro double de dimens„o N.           *
     *******************************************************************/
    double *x;
    int i;
    
    x = (double *)malloc(N*sizeof(double));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<N;i++)
        x[i] = 0.0;
    
    return (x);
}

double **aloca_matrizd(int N, int K){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro double de dimens„o N x K         *
     *******************************************************************/
    int i,j;
    double **x;
    
    x = (double **)malloc(N*sizeof(double*));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<N;i++){
        x[i] =(double *)malloc(K * sizeof(double));
        if(!x[i]) {
            printf("\nFalha na alocaÁ„o de memÛria\n");
            exit(1);
        }
    }
    
    for(i=0;i<N;i++)
        for(j=0;j<K;j++)
            x[i][j] = 0.0;
    
    return (x);
}

double ***aloca_matriz3d(int n, int N, int K){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro double de dimens„o n x N x K     *
     *******************************************************************/
    int i,j,c;
    double ***x;
    
    x = (double ***)malloc(n*sizeof(double*));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<n;i++)
    {
        x[i] =(double **)malloc(N * sizeof(double));
        if(!x[i])
        {
            printf("\nFalha na alocaÁ„o de memÛria\n");
            exit(1);
        }
        
        for(j=0;j<N;j++)
        {
            x[i][j] =(double *)malloc(K * sizeof(double));
            if(!x[i][j])
            {
                printf("\nFalha na alocaÁ„o de memÛria\n");
                exit(1);
            }
        }
    }
    
    for(i=0;i<n;i++)
        for(j=0;j<N;j++)
            for(c=0;c<K;c++)
                x[i][j][c] = 0.0;
    
    return (x);
}

double **aloca_matriz_double(int N, int K){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro double de dimens„o N x K         *
     *******************************************************************/
    int i,j;
    double **x;
    
    x = (double **)malloc(N*sizeof(double*));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<N;i++){
        x[i] =(double *)malloc(K * sizeof(double));
        if(!x[i]) {
            printf("\nFalha na alocaÁ„o de memÛria\n");
            exit(1);
        }
    }
    
    for(i=0;i<N;i++)
        for(j=0;j<K;j++)
            x[i][j] = 1.0;
    
    return (x);
}

int **aloca_matriz_int(int N, int K){
    /*******************************************************************
     * FunÁ„o que  retorna um ponteiro int de dimens„o N x K            *
     *******************************************************************/
    int i,j;
    int **x;
    
    x = (int **)malloc(N*sizeof(int*));
    if(!x){
        printf("\nFalha na alocaÁ„o de memÛria\n");
        exit(1);
    }
    
    for(i=0;i<N;i++){
        x[i] =(int *)malloc(K * sizeof(int));
        if(!x[i]) {
            printf("\nFalha na alocaÁ„o de memÛria\n");
            exit(1);
        }
    }
    
    for(i=0;i<N;i++)
        for(j=0;j<K;j++)
            x[i][j] = 0;
    
    return (x);
}


void desaloca_matrizd(double **pt, int N){
    /******************************************************************
     * Desaloca o ponteiro da matriz double (dim N x K).               *
     ******************************************************************/
    
    int i;

    for(i=0;i<N;i++)
        free(pt[i]);
    
    free(pt);
}

void desaloca_matriz3d(double ***pt, int n, int N){
    /******************************************************************
     * Desaloca o ponteiro da matriz double (dim n x N x K).               *
     ******************************************************************/
    
    int i, j;
    
    for(i=0;i<n;i++)
    {
        for(j=0;j<N;j++)
            free(pt[i][j]);
        
        free(pt[i]);
    }
    
    free(pt);
    
    return;
}

void desaloca_matriz_int(int **pt, int N){
    /******************************************************************
     * Desaloca o ponteiro da matriz int (dim N x K).                  *
     ******************************************************************/
    
    int i;
    
    for(i=0;i<N;i++)
        free(pt[i]);
    free(pt);
    
    return;
}

#endif /* _ALOCACAO_H_ */
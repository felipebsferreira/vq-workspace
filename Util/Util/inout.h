#ifndef	_INOUT_H_
#define	_INOUT_H_	1

#include "alocacao.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

double **ler_arquivo(const char *nome, int N, int K);
double **ler_arquivo(const char *nome, int N, int qtd, int K);
double **ler_arquivo_treino(const char *nome, int *nvet, int K);
double *ReadFileByPattern(const char *nome, int *qtyPixels);

int **LoadImagePixels(const char *name, int *sizeX, int *sizeY);

void escreve_arquivo(const char *nome, double **x, int N, int K, const char *tipo);
void escreve_arquivo_float(char nome[], float **x, int N, int K);
void escreve_arquivo_lsf(char nome[], double **x, int N, int K);
void escreve_arquivo_float_niter(char nome[], float **x, int N, int K, int niter);

double **ler_arquivo(const char *nome, int N, int K){
/*******************************************************************
* Função que ler retorna um ponteiro de dimenão N x K com os veto-*
* res lidos a partir do dicionário armazenado no arquivo "nome"    *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;
    double **x;

    /* Leitura dos vetores do dicionário a partir do arquivo */

    arquivo = fopen(nome, "r");
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    /* Alocação de memória para o vetor */

    x = aloca_matrizd(N,K);
     
    /* Leitura propriamente dita */

    for(i = 0; i < N; i++)
        for(j = 0; j < K; j++)
            fscanf(arquivo,"%lf",(x[i]+j));
    
    fclose(arquivo);

    return (x);
}

double **ler_arquivo(const char *nome, int N, int qtd, int K){
    /*******************************************************************
     * Função que ler retorna um ponteiro de dimens„o N x K com os veto-*
     * res lidos a partir do dicion·rio armazenado no arquivo "nome"    *
     *******************************************************************/
    
    FILE *arquivo;
    int h, i, j;
    double **x;
    
    /* Leitura dos vetores do dicion·rio a partir do arquivo */
    
    arquivo = fopen(nome, "r");
    if (!arquivo) {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }
    
    /* Alocação„o de memória para o vetor */
    
    x = aloca_matrizd(N * qtd, K);
    
    /* Leitura propriamente dita */
    
    for (h = 0; h < qtd; h++)
    {
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < K; j++)
            {
                fscanf(arquivo,"%lf",(&x[(h * N) + i][j]));
            }
        }
    }
    
    fclose(arquivo);
    
    return (x);
}


double **ler_arquivo_treino(const char *nome, int *nvet, int K){
/*******************************************************************
* Função que ler retorna um ponteiro de dimens„o nvet,K com os veto*
* res lidos a partir do dicion·rio armazenado no arquivo "nome".   *
* A Função tambÈm retorna o n˙mero de vetores k-dimensionais do ar *
* quivo de treino atravÈs da vari·vel N.                           *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;
    double **y,aux;

    /* Leitura dos vetores do dicion·rio a partir do arquivo */

    arquivo = fopen(nome, "r");
    
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    /* Leitura dos vetores de treino:  
       1.o Passo) Faz-se a leitura do arquivo, sÛ para contar o n˙mero de
                vetores necess·rios;
       2.o Passo) Aloca-se ent„o memÛria para os vetores;
        3.o Passo) Faz-se a leitura dos vetores
    */

    /* 1.o Passo */

    *nvet = 0;
    
    for(i = 1; !feof(arquivo); i++)
    {
        fscanf(arquivo,"%lf",&aux);
        if(feof(arquivo)) break;
      
        for(j=1;j<K;j++)
            fscanf(arquivo,"%lf",&aux);
        
        *nvet = i;
    }

    rewind(arquivo); /* Retorna o ponteiro para o inÌcio do arquivo */

    /* 2.o Passo */

    /* AlocaÁ„o de memÛria para os vetores de treino */
    /* nvet vetores de dimens„o K */

    y = aloca_matrizd(*nvet,K); 
        
    /* 3.o Passo */

    for(i=0;i<*nvet;i++)
        for(j=0;j<K;j++)
            fscanf(arquivo,"%lf",(y[i]+j));
    
    fclose(arquivo);

    return (y);
}

double *ReadFileByPattern(const char *nome, int *qtyPixels)
{
    /*******************************************************************
     * Função que ler retorna um ponteiro de dimens„o nvet,K com os veto*
     * res lidos a partir do dicion·rio armazenado no arquivo "nome".   *
     * A Função tambÈm retorna o n˙mero de vetores k-dimensionais do ar *
     * quivo de treino atravÈs da vari·vel N.                           *
     *******************************************************************/
    
    FILE *arquivo;
    int i;
    double *y,aux;
    char tempChar[50];
    int tempInt;
    
    /* Leitura dos vetores do dicionário a partir do arquivo */
    
    arquivo = fopen(nome, "r");
    
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }
    
    fscanf(arquivo, "%s %i %i %i ", tempChar, &tempInt, &tempInt, &tempInt);
    
    /* Leitura dos vetores de treino:
     1.o Passo) Faz-se a leitura do arquivo, sÛ para contar o n˙mero de
     vetores necess·rios;
     2.o Passo) Aloca-se ent„o memÛria para os vetores;
     3.o Passo) Faz-se a leitura dos vetores
     */
    
    /* 1.o Passo */
    
    *qtyPixels = 0;
    
    for(i = 0; !feof(arquivo); i++)
    {
        fscanf(arquivo,"%lf",&aux);
    }
    
    *qtyPixels = i;
    
    rewind(arquivo); /* Retorna o ponteiro para o inÌcio do arquivo */
    
    fscanf(arquivo, "%s %i %i %i ", tempChar, &tempInt, &tempInt, &tempInt);
    
    /* 2.o Passo */
    
    /* AlocaÁ„o de memÛria para os vetores de treino */
    /* nvet vetores de dimens„o K */
    
    y = aloca_vetor_double(*qtyPixels);
    
    /* 3.o Passo */
    
    for(i = 0; i < *qtyPixels; i++)
        fscanf(arquivo, "%lf", &y[i]);
    
    fclose(arquivo);
    
    return (y);
}

int **LoadImagePixels(const char *name, int *sizeX, int *sizeY)
{
    int **pixels;
    int tempInt;
    char tempChar[50];
    FILE *file;
    
    file = fopen(name, "r");
    fscanf(file, "%s %i %i %i ", tempChar, sizeX, sizeY, &tempInt);
    
    pixels = aloca_matriz_int(*sizeX, *sizeY);
    
    for (int i = 0; i < *sizeX; i++)
    {
        for (int j = 0; j < *sizeY; j++)
        {
            fscanf(file, "%i ", &pixels[i][j]);
        }
    }
    
    fclose(file);
    
    return pixels;
}

//Entrada:
//    nome - Nome do arquivo
//    x - Matriz a ser escrita
//    N - Quantidade de vetores
//    K - Dimensão
//    tipo - Tipo de abertura do arquivo ("w" ou "a")
void escreve_arquivo(const char *nome, double **x, int N, int K, const char *tipo)
{
    FILE *arquivo;
    int i, j;

    arquivo = fopen(nome, tipo);
        
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < K; j++)
            fprintf(arquivo,"%f ",x[i][j]);
        
        fprintf(arquivo,"\n");
    }

    fclose(arquivo);
}

void escreve_arquivo_float(char nome[], float **x, int N, int K){
/*******************************************************************
* Função que escreve os N vetores K-dimensionais armazenados no    *
* ponteiro x no arquivo "nome".                                    *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;

    arquivo = fopen(nome, "w");
    
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    for(i = 0; i < N; i++)
    {
        for( j = 0; j < K; j++)
            fprintf(arquivo,"%f ",x[i][j]);

        fprintf(arquivo,"\n");
    }

    fclose(arquivo);

    return;
}

void escreve_arquivo_float_niter(char nome[], float **x, int N, int K, int niter){
/*******************************************************************
* Função que escreve os N vetores K-dimensionais armazenados no    *
* ponteiro x no arquivo "nome".                                    *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;
    char teste[40], iteracao[10];

    teste[0] = '\0';

    sprintf(iteracao,"%d",niter);

    strcat(teste,nome);
    strcat(teste,iteracao);

    arquivo = fopen(teste, "w");
    
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    for(i=0;i<N;i++){
      for(j=0;j<K;j++)
         fprintf(arquivo,"%f ",x[i][j]);
      fprintf(arquivo,"\n");
    }

    fclose(arquivo);

    return;
}

void escreve_arquivo_double_niter(char nome[], double **x, int N, int K, int niter){
/*******************************************************************
* Função que escreve os N vetores K-dimensionais armazenados no    *
* ponteiro x no arquivo "nome".                                    *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;
    char teste[40], iteracao[10];

    teste[0] = '\0';

    sprintf(iteracao, "%d", niter);

    strcat(teste,nome);
    strcat(teste,iteracao);

    arquivo = fopen(teste, "w");
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < K; j++)
            fprintf(arquivo,"%f ", x[i][j]);
        
        fprintf(arquivo, "\n");
    }

    fclose(arquivo);

    return;
}

void escreve_arquivo_lsf(char nome[], double **x, int N, int K){
/*******************************************************************
* Função que escreve os N vetores K-dimensionais armazenados no    *
* ponteiro x no arquivo "nome".                                    *
*******************************************************************/ 

    FILE *arquivo;
    int i, j;

    arquivo = fopen(nome, "w");
    
    if (!arquivo)
    {
        printf("\nErro na abertura do arquivo %s!!!\a\n",nome);
        exit(1);
    }

    for(i = 0; i < N; i++)
        for(j = 0; j < K; j++)
            fprintf(arquivo,"%.14f\n",x[i][j]);

    fclose(arquivo);

    return;
}

#endif	/* _INOUT_H_  */

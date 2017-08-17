//
//  codec.h
//  Util
//
//  Created by Felipe Ferreira on 15/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _CODEC_H_
#define _CODEC_H_

#include <math.h>
#include "alocacao.h"
#include "inout.h"
#include "general_vq.h"

void codificar(char entrada[], char saida[], double **x, int N, int K);
void decodificar(char entrada[], char saida[], double **x, int N, int K);

int *codificar(double **dic, double **treino, int N, int K, int nTreino);
double **decodificar(double **dic, int *p, int K, int N, int n);

int buscar_vmp(double **pta, double *y, int N, int K);

void dat_to_pgm(double **matriz, int x, int y, const char *nome, int Ximg, int Yimg);

double * gerar_imagem_erro(const char * img_original, const char * img_quantizada, const char * img_erro);

void codificar(char entrada[], char saida[], double **x, int N, int K)
{
    /******************************************************************
     * Função que codifica e decodifica os vetores do arquivo de nome  *
     * "entrada", de acordo com o dicionário "x" de dimensão N,K.      *
     ******************************************************************/
    
    FILE *entradap, *saidap;
    int i, j, indice;
    double *y;
    
    y = aloca_vetor_double(K);
    
    entradap = fopen(entrada,"r");
    if (!entradap)
    {
        printf("\nErro na abertura do arquivo dos dados a serem codificados\n");
        exit(1);
    }
    
    saidap = fopen(saida,"w");
    if (!saidap)
    {
        printf("\nErro na abertura do arquivo dos dados codificados\n");
        exit(1);
    }
    
    // Monta o vetor-codigo apartir do arquivo.dat original
    for(i=0; !feof(entradap); i+=K)
    {
        fscanf(entradap, "%lf", &y[0]);
        
        if(feof(entradap))
            break;
        
        for(j=1;j<K;j++)
        {
            fscanf(entradap,"%lf",&y[j]);
        }
        
        indice = buscar_vmp(x,y,N,K);
        fprintf(saidap,"%d\n",indice);
    }
    
    fclose(saidap);
    fclose(entradap);
    free(y);
}

void decodificar(char entrada[], char saida[], double **x, int N, int K)
{
    /******************************************************************
     * FunÁ„o que codifica e decodifica os vetores do arquivo de nome  *
     * "entrada", de acordo com o dicion·rio "x" de dimens„o N,K.      *
     ******************************************************************/
    
    FILE *entradap, *saidap;
    int i, j, indice;
    double *y;
    
    y = aloca_vetor_double(K);
    
    entradap = fopen(entrada,"r");
    if (!entradap)
    {
        printf("\nErro na abertura do arquivo dos dados a serem codificados\n");
        exit(1);
    }
    
    saidap = fopen(saida,"w");
    if (!saidap)
    {
        printf("\nErro na abertura do arquivo dos dados codificados\n");
        exit(1);
    }
    
    for(i = 0; !feof(entradap); i += K)
    {
        fscanf(entradap, "%d", &indice);
        
        for(j = 0; j < K; j++)
            fprintf(saidap, "%f\n", x[indice][j]);
        
        for(j = 0; j < K; j++)
            y[j] = 0.0;
    }
    
    fclose(saidap);
    fclose(entradap);
    free(y);
}

int *codificar(double **dic, double **treino, int N, int K, int nTreino)
{
    int *particao;
    
    particao = aloca_vetor_int(nTreino);
    
    for (int i = 0; i < nTreino; i++)
    {
        particao[i] = buscar_vmp(dic, treino[i], N, K);
    }
    
    return particao;
}

double **decodificar(double **dic, int *p, int K, int N, int n)
{
    double **data;
    int i, j;
    
    data = aloca_matrizd(n, K);
    
    for(i = 0; i < n; i++)
        for(j = 0; j < K; j++)
            data[i][j] = dic[p[i]][j];
    
    return data;
}

//Retorna o índice do vetor-código  que produz a menor distorção
//Dicionário (pta) - Vetor (y) - Número de vetores do dicionário (N)
int buscar_vmp(double **dic, double *vetor_treino, int N, int K)
{
    int i, cod;
    double dist, dist_atual;
    
    dist = calcular_distancia_euclidiana(dic[0], vetor_treino, K);
    cod = 0;
    
    for(i = 1; i < N; i++)
    {
        dist_atual = calcular_distancia_euclidiana(dic[i], vetor_treino, K);
        
        if(dist > dist_atual)
        {
            dist = dist_atual;
            cod = i;
        }
    }
    
    return cod;
}

// Converter matriz em blocos para imagem no formato PGM
//
// double **matriz -> Dados quantizados
// int x           -> Tamanho x do bloco
// int y           -> Tamanho y do bloco
// int Ximg        -> Dimensão X da imagem
// int Yimg        -> Dimensão Y da imagem
void dattopgm(double **matriz, int x, int y, const char *nome, int Ximg, int Yimg)
{
    int i,j,k,l;
    FILE *arquivo;
    int **saida, *aux;
    
    saida = aloca_matriz_int(Yimg, Ximg);
    aux = aloca_vetor_int(Ximg*Yimg);
    
    for(i = 0; i < Yimg; i++)
    {
        for (j = 0; j < Ximg; j++)
        {
            *(aux+i) = (int)round(matriz[i][j]);
        }
    }
    
    // reconstrução da imagem a partir dos vetores de entrada
    for(i=0;i<64;i++)
        for(j=0;j<64;j++)
            for(k=0;k<4;k++)
                for(l=0;l<4;l++)
                    *(saida[4*i+k] + (4*j+l)) = *(aux + 1024*i + 16*j + 4*k + l);
    
    //armazenamento da imagem reconstruida em arquivo
    arquivo = fopen(nome, "w");
    
    // Cabeçalho .pgm
    fprintf(arquivo,"P2 %i %i 255 ", Ximg, Yimg);
    
    for(i = 0; i < Yimg; i++)
    {
        for(j = 0; j < Ximg; j++)
        {
            fprintf(arquivo, "%d ", saida[i][j]);
        }
    }
    
    fclose(arquivo);
}

// Converter matriz em blocos para imagem no formato PGM
//
// double **matriz -> Dados quantizados
// int x           -> Tamanho x do bloco
// int y           -> Tamanho y do bloco
// int Ximg        -> Dimensão X da imagem
// int Yimg        -> Dimensão Y da imagem
void dat_to_pgm(double **matriz, int x, int y, const char *nome, int Ximg, int Yimg)
{
    int xaux, yaux, i, j, k, l, qtd;
    int **saida, *vet;
    FILE *arq;
    
    /*Alocação de memória e preenchimento da matriz imagem*/
    saida = aloca_matriz_int(Yimg, Ximg);
    
    /*Criando arquivo de saida*/
    arq = fopen(nome, "w");
    
    /*Inserindo cabeçalho*/
    fprintf(arq, "P2 %d %d 255 ", Ximg, Yimg);
    
    vet = new int[x*y];
    xaux = yaux = qtd = 0;
    
    /*Dividindo o arquivo em blocos x por y*/
    for(i = 0; i < Yimg; i++)
    {
        for(j = 0; j < Ximg; j++)
        {
            vet[qtd] = (int)round(matriz[i][j]);
            
            if(qtd == x*y - 1)
            {
                qtd = 0;
                
                for(k = 0; k < y; k++)
                {
                    for(l = 0; l < x; l++)
                    {
                        saida[yaux + k][xaux + l] = vet[qtd];
                        qtd++;
                    }
                }
                
                qtd = 0;
                xaux += x;
                
                if(xaux == Ximg)
                {
                    xaux = 0;
                    yaux += y;
                }
            }
            else
                qtd++;
        }
    }
    
    /*Escrevendo dados no arquivo .pgm*/
    for(i = 0; i < Yimg; i++)
        for(j = 0; j < Ximg; j++)
            fprintf(arq, "%d ", saida[i][j]);
    
    fclose(arq);
    delete vet;
    desaloca_matriz_int(saida, Yimg);
}

// Descricão: Gerar imagem erro (.PGM) com base na imagem original (.PGM) e quantizada (.PGM)
//
// Entrada:
// const char * img_original   -> nome da imagem original (.PGM)
// const char * img_quantizada -> nome da imagem quantizada (.PGM)
// const char * img_erro       -> nome da imagem erro gerada (.PGM)
//
// Saída:
// double *retorno
//      -> retorno[0] : Valor do píxel erro médio
//      -> retorno[1] : Valor do desvio padrão dos erros
double * gerar_imagem_erro(const char * img_original, const char * img_quantizada, const char * img_erro)
{
    double *retorno = new double[3];
    
    FILE *arq_original, *arq_quantizado, *arq_erro;
    int pixel_original, pixel_quantizado, pixel_erro;
    double pixel_medio = 0, desvio_padrao = 0;
    int *pixels_erro;
    int i, j, indice, max;
    char str[50];
    
    // Carregar os arquivos original e quantizado
    arq_original = fopen(img_original, "r");
    arq_quantizado = fopen(img_quantizada, "r");
    
    fscanf(arq_original, "%s %i %i %i ", str, &i, &j, &indice);
    fscanf(arq_quantizado, "%s %i %i %i ", str, &i, &j, &indice);
    
    arq_erro = fopen(img_erro, "w");
    
    fprintf(arq_erro, "P2 %i %i %i ", i, j, indice);
    
    pixels_erro = new int[i*j];
    indice = 0;
    
    while(!feof(arq_original) && !feof(arq_quantizado))
    {
        fscanf(arq_original, "%i ", &pixel_original);
        fscanf(arq_quantizado, "%i ", &pixel_quantizado);
        
        pixel_erro = abs(pixel_original - pixel_quantizado);
        
        if (pixel_erro > -1)
        {
            pixel_medio += pixel_erro;
            pixels_erro[indice] = pixel_erro;
            indice++;
        }
        	
        fprintf(arq_erro, "%i ", pixel_erro);
    }
    
    max = indice;

    pixel_medio = pixel_medio / (double)max;
    
    for (indice = 0; indice < max; indice++)
    {
        desvio_padrao += (pixels_erro[indice] - pixel_medio) * (pixels_erro[indice] - pixel_medio);
    }
    
    desvio_padrao = sqrt(desvio_padrao / (double)max);
    
    retorno[0] = pixel_medio;
    retorno[1] = desvio_padrao;
    retorno[2] = (desvio_padrao / pixel_medio) * 100;
    
    return retorno;
}

#endif /* _CODEC_H_ */

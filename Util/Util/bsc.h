//
//  bsc.h
//  Util
//
//  Finalidade: Faz a transmissão não codificada de índices intei-
//  ros de um arquivo de entrada. A transmissão é feita em um canal
//  binário simétrico (BSC) com taxa de cruzamento igual a p.
//
//  Created by Felipe Ferreira on 15/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef bsc_h
#define bsc_h

#include <iostream>
#include <math.h>
#include "rand.h"

using namespace std;

void bsc(char [], char [], double);
int *bsc(int *, int, double);
int procurar_maior(char []);
int calcular_qtd_bits(int);
int transmitir(int, int, double);

void bsc(char name1[], char name2[], double perro){
    /******************************************************************
     * Função que transmite os índices decimais do arquivo "name1" atra*
     * vés de um canal binário simétrico com probabilidade de cruzamen-*
     * to igual a perro. Os índices recebidos são armazenados no arqui-*
     * vo name2.                                                       *
     ******************************************************************/

    int maior, nbits, orig;
    FILE *fp1, *fp2;
    int soma, receb;

    maior = procurar_maior(name1);
    nbits = calcular_qtd_bits(maior);

    /* Abertura do arquivo de entrada */

    fp1 = fopen(name1,"r");
    if(!fp1){
        printf("\nErro na abertura do arquivo de entrada!\n");
        exit(1);
    }

    /* Abertura do arquivo de saída */

    fp2 = fopen(name2,"w");
    if(!fp2){
        printf("\nErro na abertura do arquivo de saída!\n");
        exit(1);
    }

    soma=0;

    for(;!feof(fp1);)
    {
        fscanf(fp1,"%d\n",&orig);

        receb=transmitir(orig,nbits,perro);

        fprintf(fp2,"%d\n",receb);

        if(orig!=receb) soma++;
    }

    printf("Número de símbolos errados: %d\n",soma);

    fclose(fp1);
    fclose(fp2);

    return;
}

int *bsc(int *particao, int length, double perro){
    /******************************************************************
     * Função que transmite os índices decimais do arquivo "name1" atra*
     * vés de um canal binário simétrico com probabilidade de cruzamen-*
     * to igual a perro. Os índices recebidos são armazenados no arqui-*
     * vo name2.                                                       *
     ******************************************************************/
    
    int nbits, maior = 0, soma = 0;
    int *transmitido = aloca_vetor_int(length);

    for(int i = 0; i < length; i++)
    {
        if(particao[i] > maior)
        {
            maior = particao[i];
        }
    }

    nbits = calcular_qtd_bits(maior);

    for(int i = 0; i < length; i++)
    {
        transmitido[i] = transmitir(particao[i], nbits, perro);
        if(transmitido[i] != particao[i]) soma++;
    }

    cout << "Número de símbolos errados: " << soma << endl;

    return transmitido;
}

int transmitir(int orig, int nbits, double perro){
    /******************************************************************
     * Função que transmite o inteiro orig através um canal binário si-*
     * métrico, com probabilidade de cruzamento perro, utilizando nbits*
     * para representar o inteiro nbits.                               *
     ******************************************************************/
    int i, aux, receiv;

    aux = 0;

    for(i=0;i<nbits;i++)
    {
        if(unif() < perro)
            aux += pow(2,i);
    }


    receiv = aux^orig;

    return(receiv);

}

int procurar_maior (char name[]){
    /******************************************************************
     * Função que retorna o maior valor inteiro armazenada no arquivo  *
     * de nome name.                                                   *
     ******************************************************************/

    int aux, maior;
    FILE *fp;

    fp = fopen(name,"r");

    if(!fp){
        printf("\nErro na abertura do arquivo de entrada!\n");
        exit(1);
    }

    fscanf(fp,"%d",&maior);

    for(;!feof(fp);){
        fscanf(fp,"%d\n",&aux);
        if(aux > maior)
            maior = aux;
    }

    fclose(fp);

    return(maior);
}

int calcular_qtd_bits(int NPT){
    /******************************************************************
     * Determina o número de bits necessários para representar o intei-*
     * ro NPT.                                                         *
     ******************************************************************/

    int i, nbits = 0;

    for(i=1;i<20;i++)
        if((NPT&((int)pow(2,i)))>>i)
            nbits = i+1;

    return nbits;
}

#endif /* bsc_h */

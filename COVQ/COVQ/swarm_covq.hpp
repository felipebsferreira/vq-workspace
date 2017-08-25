//
//  swarm_covq.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 09/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef SWARM_COVQ_H
#define SWARM_COVQ_H

#include <iostream>
#include <cmath>
#include "../../Util/Util/rand.h"
#include "../../Util/Util/inout.h"
#include "../../Util/Util/definicoes.h"
#include "../../Util/Util/general_vq.h"

using namespace std;

class SWARM_COVQ
{
public:
    
    int qtd_part;
    int N;
    int K;
    int nTreino;
    double erro;
    SCALING scaleType;
    double scale;
    SWARM_TYPE swarm_type;
    
    virtual double * Run(double **dic, double **treino, NNS nns);
    
    // Construtores
    SWARM_COVQ();
    SWARM_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
    SWARM_COVQ(int N, int K, int qtd_part, SWARM_TYPE swarm_type, double erro, SCALING scaleType, double scale);
    
    // Destrutores
    ~SWARM_COVQ();
    
protected:
    double **gbest;
    double dist_gbest;
    
    void AtualizarCentroides(double **treino, double **dic, int *particao, int iteracao);
    double CalcularDistorcaoMedia(double **treino, double **dic, int *particao);
    double CalcularDEQ(double *vetor_treino, double *vetor_dic);
    
    void Classificar(double **treino, double **dic, int *particao, NNS nns);
    
    void CalcularParamTreino(double **vector, NNS nns);
    void DesalocarParams(NNS nns);
    
    void OrdenarDicionarioPorMedias(double **dic);
    
private:
    Params params;
    double **matriz_prob;
    
    void ClassificarPorBT(double **treino, double **dic, int *particao);
    void ClassificarPorPDS(double **treino, double **dic, int *particao);
    void ClassificarPorENNS(double **treino, double **dic, int *particao);
    void ClassificarPorIENNS(double **treino, double **dic, int *particao);
    void ClassificarPorEENNS(double **treino, double **dic, int *particao);
    void ClassificarPorIEENNS(double **treino, double **dic, int *particao);
    void ClassificarPorEEENNS(double **treino, double **dic, int *particao);
    
    void CalcularParamENNS(double **vector, double *mean, int num_vector);
    void CalcularParamIENNS(double **vector, double *Mean, double *Sum_f, double *Sum_s, int num_vector);
    void CalcularParamEENNS(double **vector, double *Mean, double *Variance, int num_vector);
    void CalcularParamEEENNS(double **vector, double *Mean, double *Variance, double *Norm, int num_vector);
    void Ordenar(int *Sort_Mean_c, int *Sort_Mean_inv, double *Mean_c);
    void ClassificarMedias(int *Classifica_Dif, int *Sort_Mean_c, int *Sort_Mean_inv, double *mean_x, double *Mean_c);
};

SWARM_COVQ::SWARM_COVQ()
{
    this->gbest = NULL;
    this->matriz_prob = NULL;
}

SWARM_COVQ::SWARM_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale)
{
    this->N = N;
    this->K = K;
    this->qtd_part = qtd_part;
    this->erro = erro;
    this->scaleType = scaleType;
    this->scale = scale;
    this->matriz_prob = montar_matriz_probabilidades(N, erro);
    this->gbest = NULL;
}

SWARM_COVQ::SWARM_COVQ(int N, int K, int qtd_part, SWARM_TYPE swarm_type, double erro, SCALING scaleType, double scale)
{
    this->N = N;
    this->K = K;
    this->qtd_part = qtd_part;
    this->erro = erro;
    this->scaleType = scaleType;
    this->scale = scale;
    this->matriz_prob = montar_matriz_probabilidades(N, erro);
    this->gbest = NULL;
    this->swarm_type = swarm_type;
}

SWARM_COVQ::~SWARM_COVQ()
{
    if (matriz_prob != NULL)
        desaloca_matrizd(matriz_prob, N);
    
    if (gbest != NULL)
        desaloca_matrizd(gbest, N);
}

double * SWARM_COVQ::Run(double **dic, double **treino, NNS nns)
{
    cout << "Este metodo deve ser executado por uma classe herdada!" << endl;
    exit(1);
}

void SWARM_COVQ::AtualizarCentroides(double **treino, double **dic, int *particao, int iteracao)
{
    int i, n, j;
    int num_treino[N]; // Vetor que armazena a quantidade de vetores x de uma região Si
    double denominador, s = 0;
    double numerador[K];
    double soma_treino[N][K]; // Matriz que contém a soma dos vetores x pertencentes a região Si
    
    if(scaleType == FIXED)
    {
        s = scale;
    }
    else if(scaleType == VARIABLE)
    {
        s = 1 + (double)scale/(scale + iteracao);
    }
    
    // Zerando a matriz soma_x e o vetor num_x para o novo cálculo dos centróides
    for(i = 0; i < N; i++)
    {
        num_treino[i] = 0;
        
        for(j = 0; j < K; j++)
        {
            soma_treino[i][j] = 0;
        }
    }
    
    // Soma dos vetores x pertencentes a região Si
    // Soma da quantidade dos vetores x pertencentes a região Si
    for(n = 0; n < nTreino; n++)
    {
        for(i = 0; i < N; i++)
        {
            if(particao[n] == i)
            {
                num_treino[i]++;
                
                for(j = 0; j < K; j++)
                {
                    soma_treino[i][j] += treino[n][j];
                }
            }
        }
    }
    
    // Atualização dos Centróides
    for(j = 0; j < N; j++)
    {
        denominador = 0;
        
        for(n = 0; n < K; n++)
        {
            numerador[n] = 0; // Prepara o numerador para o próximo cálculo de centróide
        }
        
        for(i = 0; i < N; i++)
        {
            denominador += matriz_prob[i][j] * num_treino[i];
            
            for(n = 0; n < K; n++)
            {
                numerador[n] += (matriz_prob[i][j] * soma_treino[i][n]);
            }
        }
        
        for(n = 0; n < K; n++)
        {
            if(scaleType != NONE)
            {
                dic[j][n] = dic[j][n] + s * ((numerador[n] / denominador) - dic[j][n]);
            }
            else
            {
                dic[j][n] = (numerador[n] / denominador);
            }
        }
    }
}

double SWARM_COVQ::CalcularDistorcaoMedia(double **treino, double **dic, int *particao)
{
    int i, j;
    double dist_total, Dcovq;
    
    dist_total = 0.0;
    
    for(i = 0; i < nTreino; i++)
    {
        for(j = 0; j < N; j++)
        {
            dist_total += matriz_prob[particao[i]][j] * CalcularDEQ(treino[i], dic[j]);
        }
    }
    
    Dcovq = dist_total / (K * nTreino);
    
    return Dcovq;
}

// Calcular distância Euclidiana quadrática
double SWARM_COVQ::CalcularDEQ(double *vetor_treino, double *vetor_dic)
{
    double distq;
    int i;
    
    distq = 0.0;

    for(i = 0; i < K; i++)
    {
        distq += (vetor_treino[i] - vetor_dic[i]) * (vetor_treino[i] - vetor_dic[i]);
    }

    return distq;
}

void SWARM_COVQ::OrdenarDicionarioPorMedias(double **dic)
{
    double medias[N], aux;
    int i, j, k;
    
    for (i = 0; i < N; i++)
    {
        medias[i] = 0;
        
        for (j = 0; j < K; j++)
        {
            medias[i] += dic[i][j];
        }
        
        medias[i] = medias[i] / K;
    }
    
    for(i = 0; i < N - 1; i++)
    {
        for(j = i + 1; j < N; j++)
        {
            if(medias[i] > medias[j])
            {
                aux = medias[i];
                medias[i] = medias[j];
                medias[j] = aux;
                
                for (k = 0; k < K; k++)
                {
                    aux = dic[i][k];
                    dic[i][k] = dic[j][k];
                    dic[j][k] = aux;
                }
            }
        }
    }
}

void SWARM_COVQ::Classificar(double **treino, double **dic, int *particao, NNS nns)
{
    switch (nns)
    {
        case BT:
            ClassificarPorBT(treino, dic, particao);
            break;
            
        case PDS:
            ClassificarPorPDS(treino, dic, particao);
            break;
            
        case ENNS:
            ClassificarPorENNS(treino, dic, particao);
            break;
            
        case IENNS:
            ClassificarPorIENNS(treino, dic, particao);
            break;
            
        case EENNS:
            ClassificarPorEENNS(treino, dic, particao);
            break;
            
        case IEENNS:
            ClassificarPorIEENNS(treino, dic, particao);
            break;
            
        case EEENNS:
            ClassificarPorEEENNS(treino, dic, particao);
            break;
            
        default:
            break;
    }
}

void SWARM_COVQ::ClassificarPorBT(double **treino, double **dic, int *particao)
{
    int i, j, n;
    double dist, dist_aux;
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[0][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = 0;
        }
        
        for(i = 1; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                dist_aux += matriz_prob[i][j] * CalcularDEQ(treino[n], dic[j]);
            }
            
            if(dist > dist_aux)
            {
                dist = dist_aux;
                particao[n] = i;
            }
            
            dist_aux = 0;
        }
    }
    
}

void SWARM_COVQ::ClassificarPorPDS(double **treino, double **dic, int *particao)
{
    int i, j, n;
    double dist, dist_aux;
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 1e20;
        dist_aux = 0;
        
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                dist_aux += matriz_prob[i][j] * CalcularDEQ(treino[n], dic[j]);
                
                if(dist_aux > dist) break;
            }
            
            if(dist > dist_aux)
            {
                dist = dist_aux;
                particao[n] = i;
            }
            
            dist_aux = 0;
        }
    }
}

void SWARM_COVQ::ClassificarPorENNS(double **treino, double **dic, int *particao)
{
    int j, n, count_up, count_down;
    double dist, dist_aux;
    double mean_c[N];
    int sort_mean_c[N], sort_mean_inv[N], classifica_dif[nTreino];
    
    CalcularParamENNS(dic, mean_c, N); // Média dos vetores-código
    Ordenar(sort_mean_c, sort_mean_inv, mean_c);
    ClassificarMedias(classifica_dif, sort_mean_c, sort_mean_inv, params.mean_x, mean_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = classifica_dif[n];
        count_down = classifica_dif[n];
        
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[sort_mean_c[classifica_dif[n]]][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = sort_mean_c[classifica_dif[n]];
        }
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(mean_c[sort_mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    dist_aux = 0;
                    
                    for(j = 0; j < N; j++)
                    {
                        dist_aux += matriz_prob[sort_mean_c[count_down]][j] * CalcularDEQ(treino[n], dic[j]);
                        
                        if(dist_aux > dist) break;
                    }
                    
                    if(dist > dist_aux)
                    {
                        dist = dist_aux;
                        particao[n] = sort_mean_c[count_down];
                    }
                }
                else
                {
                    count_down = 0;
                }
            }
            
            if(count_up != N)
            {
                count_up++;
                
                if(count_up != N && mean_c[sort_mean_c[count_up]] < params.mean_x[n] + sqrt(dist/K))
                {
                    dist_aux = 0;
                    
                    for(j = 0; j < N; j++)
                    {
                        dist_aux += matriz_prob[sort_mean_c[count_up]][j] * CalcularDEQ(treino[n], dic[j]);
                        
                        if(dist_aux > dist) break;
                    }
                    
                    if(dist > dist_aux)
                    {
                        dist = dist_aux;
                        particao[n] = sort_mean_c[count_up];
                    }
                }
                else
                {
                    count_up = N;
                }
            }
            
            if(count_down == 0 && count_up == N)
            {
                break;
            }
        }
    }
}

void SWARM_COVQ::ClassificarPorIENNS(double **treino, double **dic, int *particao)
{
    int j, n, count_up, count_down;
    double dist, dist_aux;
    double Mean_c[N];
    int Sort_Mean_c[N], Sort_Mean_inv[N], Classifica_Dif[nTreino];
    double Sum_cf[N], Sum_cs[N], condition_sum_f, condition_sum_s;
    
    CalcularParamIENNS(dic, Mean_c, Sum_cf, Sum_cs, N); // Médias dos vetores-código
    Ordenar(Sort_Mean_c, Sort_Mean_inv, Mean_c);
    ClassificarMedias(Classifica_Dif, Sort_Mean_c, Sort_Mean_inv, params.mean_x, Mean_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = Classifica_Dif[n];
        count_down = Classifica_Dif[n];
        
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[Sort_Mean_c[Classifica_Dif[n]]][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        }
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(Mean_c[Sort_Mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    condition_sum_f = (params.sum_xf[n] - Sum_cf[Sort_Mean_c[count_down]]) * (params.sum_xf[n] - Sum_cf[Sort_Mean_c[count_down]]);
                    condition_sum_s = (params.sum_xs[n] - Sum_cs[Sort_Mean_c[count_down]]) * (params.sum_xs[n] - Sum_cs[Sort_Mean_c[count_down]]);
                    
                    if(condition_sum_f < K*0.5*dist && condition_sum_s < K*0.5*dist)
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_down]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_down];
                        }
                    }
                }
                else
                {
                    count_down = 0;
                }
            }
            
            if(count_up != N)
            {
                count_up++;
                
                if(count_up != N && Mean_c[Sort_Mean_c[count_up]] < params.mean_x[n] + sqrt(dist/K))
                {
                    condition_sum_f = (params.sum_xf[n] - Sum_cf[Sort_Mean_c[count_up]]) * (params.sum_xf[n] - Sum_cf[Sort_Mean_c[count_up]]);
                    condition_sum_s = (params.sum_xs[n] - Sum_cs[Sort_Mean_c[count_up]]) * (params.sum_xs[n] - Sum_cs[Sort_Mean_c[count_up]]);
                    
                    if(condition_sum_f < K*0.5*dist && condition_sum_s < K*0.5*dist)
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_up]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_up];
                        }
                    }
                }
                else
                {
                    count_up = N;
                }
            }
            
            if(count_down == 0 && count_up == N)
            {
                break;
            }
        }
    }
}

void SWARM_COVQ::ClassificarPorEENNS(double **treino, double **dic, int *particao)
{
    int j, n, count_up, count_down;
    double dist, dist_aux;
    double Mean_c[N], Variance_c[N];
    int Sort_Mean_c[N], Sort_Mean_inv[N], Classifica_Dif[nTreino];
    
    CalcularParamEENNS(dic, Mean_c, Variance_c, N); // Médias dos vetores-código
    Ordenar(Sort_Mean_c, Sort_Mean_inv, Mean_c);
    ClassificarMedias(Classifica_Dif, Sort_Mean_c, Sort_Mean_inv, params.mean_x, Mean_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = Classifica_Dif[n];
        count_down = Classifica_Dif[n];
        
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[Sort_Mean_c[Classifica_Dif[n]]][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        }
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(Mean_c[Sort_Mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    if(Variance_c[Sort_Mean_c[count_down]] > params.variance_x[n] - sqrt(dist))
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_down]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_down];
                        }
                    }
                }
                else
                {
                    count_down = 0;
                }
            }
            
            if(count_up != N)
            {
                count_up++;
                
                if(count_up != N && Mean_c[Sort_Mean_c[count_up]] < params.mean_x[n] + sqrt(dist/K))
                {
                    if(Variance_c[Sort_Mean_c[count_up]] < params.variance_x[n] + sqrt(dist))
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_up]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_up];
                        }
                    }
                }
                else
                {
                    count_up = N;
                }
            }
            
            if(count_down == 0 && count_up == N)
            {
                break;
            }
        }
    }
}

void SWARM_COVQ::ClassificarPorIEENNS(double **treino, double **dic, int *particao)
{
    int j, n, count_up, count_down;
    double dist, dist_aux, variance, mean;
    double Mean_c[N], Variance_c[N];
    int Sort_Mean_c[N], Sort_Mean_inv[N], Classifica_Dif[nTreino];
    
    CalcularParamEENNS(dic, Mean_c, Variance_c, N); // Médias dos vetores-código
    Ordenar(Sort_Mean_c, Sort_Mean_inv, Mean_c);
    ClassificarMedias(Classifica_Dif, Sort_Mean_c, Sort_Mean_inv, params.mean_x, Mean_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = Classifica_Dif[n];
        count_down = Classifica_Dif[n];
    
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[Sort_Mean_c[Classifica_Dif[n]]][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        }
    
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(Mean_c[Sort_Mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    variance = (params.variance_x[n] - Variance_c[Sort_Mean_c[count_down]])
                        * (params.variance_x[n] - Variance_c[Sort_Mean_c[count_down]]);
                    mean = (params.mean_x[n] - Mean_c[Sort_Mean_c[count_down]])
                        * (params.mean_x[n] - Mean_c[Sort_Mean_c[count_down]]);
                    
                    if( K * mean + variance < dist)
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_down]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_down];
                        }
                    }
                }
                else
                {
                    count_down = 0;
                }
            }
            
            if(count_up != N)
            {
                count_up++;
                
                if(count_up != N && Mean_c[Sort_Mean_c[count_up]] < params.mean_x[n] + sqrt(dist/K))
                {
                    variance = (params.variance_x[n] - Variance_c[Sort_Mean_c[count_up]])
                        * (params.variance_x[n] - Variance_c[Sort_Mean_c[count_up]]);
                    mean = (params.mean_x[n] - Mean_c[Sort_Mean_c[count_up]])
                        * (params.mean_x[n] - Mean_c[Sort_Mean_c[count_up]]);
                    
                    if( K * mean + variance < dist)
                    {
                        dist_aux = 0;
                        
                        for(j = 0; j < N; j++)
                        {
                            dist_aux += matriz_prob[Sort_Mean_c[count_up]][j] * CalcularDEQ(treino[n], dic[j]);
                            
                            if(dist_aux > dist) break;
                        }
                        
                        if(dist > dist_aux)
                        {
                            dist = dist_aux;
                            particao[n] = Sort_Mean_c[count_up];
                        }
                    }
                }
                else
                {
                    count_up = N;
                }
            }
            
            if(count_down == 0 && count_up == N)
            {
                break;
            }
        }
    }
}

void SWARM_COVQ::ClassificarPorEEENNS(double **treino, double **dic, int *particao)
{
    int j, n, count_up, count_down;
    double dist, dist_aux;
    double Mean_c[N], Variance_c[N], Norm_c[N];
    int Sort_Mean_c[N], Sort_Mean_inv[N], Classifica_Dif[nTreino];
    
    CalcularParamEEENNS(dic, Mean_c, Variance_c, Norm_c, N); // Médias dos vetores-código
    Ordenar(Sort_Mean_c, Sort_Mean_inv, Mean_c);
    ClassificarMedias(Classifica_Dif, Sort_Mean_c, Sort_Mean_inv, params.mean_x, Mean_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = Classifica_Dif[n];
        count_down = Classifica_Dif[n];
        
        for(j = 0; j < N; j++)
        {
            dist += matriz_prob[Sort_Mean_c[Classifica_Dif[n]]][j] * CalcularDEQ(treino[n], dic[j]);
            particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        }
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(Mean_c[Sort_Mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    if(Variance_c[Sort_Mean_c[count_down]] > params.variance_x[n] - sqrt(dist))
                    {
                        if(Norm_c[Sort_Mean_c[count_down]] > params.norm_x[n] - sqrt(dist))
                        {
                            dist_aux = 0;
                            
                            for(j = 0; j < N; j++)
                            {
                                dist_aux += matriz_prob[Sort_Mean_c[count_down]][j] * CalcularDEQ(treino[n], dic[j]);
                                
                                if(dist_aux > dist) break;
                            }
                            
                            if(dist > dist_aux)
                            {
                                dist = dist_aux;
                                particao[n] = Sort_Mean_c[count_down];
                            }
                        }
                    }
                }
                else
                {
                    count_down = 0;
                }
            }
            
            if(count_up != N)
            {
                count_up++;
                
                if(count_up != N && Mean_c[Sort_Mean_c[count_up]] < params.mean_x[n] + sqrt(dist/K))
                {
                    if(Variance_c[Sort_Mean_c[count_up]] < params.variance_x[n] + sqrt(dist))
                    {
                        if(Norm_c[Sort_Mean_c[count_up]] < params.norm_x[n] + sqrt(dist))
                        {
                            dist_aux = 0;
                            
                            for(j = 0; j < N; j++)
                            {
                                dist_aux += matriz_prob[Sort_Mean_c[count_up]][j] * CalcularDEQ(treino[n], dic[j]);
                                
                                if(dist_aux > dist) break;
                            }
                            
                            if(dist > dist_aux)
                            {
                                dist = dist_aux;
                                particao[n] = Sort_Mean_c[count_up];
                            }
                        }
                    }
                }
                else
                {
                    count_up = N;
                }
            }
            
            if(count_down == 0 && count_up == N)
            {
                break;
            }
        }
    }
}

void SWARM_COVQ::CalcularParamTreino(double **vector, NNS nns)
{
    params.mean_x = aloca_vetor_double(nTreino);
    
    switch (nns)
    {
        case ENNS:
            CalcularParamENNS(vector, params.mean_x, nTreino);
            break;
            
        case IENNS:
            params.sum_xf = aloca_vetor_double(nTreino);
            params.sum_xs = aloca_vetor_double(nTreino);
            
            CalcularParamIENNS(vector, params.mean_x, params.sum_xf, params.sum_xs, nTreino);
            break;
            
        case EENNS:
            params.variance_x = aloca_vetor_double(nTreino);
            
            CalcularParamEENNS(vector, params.mean_x, params.variance_x, nTreino);
            break;
            
        case IEENNS:
            params.variance_x = aloca_vetor_double(nTreino);
            
            CalcularParamEENNS(vector, params.mean_x, params.variance_x, nTreino);
            break;
            
        case EEENNS:
            params.variance_x = aloca_vetor_double(nTreino);
            params.norm_x = aloca_vetor_double(nTreino);
            
            CalcularParamEEENNS(vector, params.mean_x, params.variance_x, params.norm_x, nTreino);
            break;
            
        default:
            free(params.mean_x);
            break;
    }
}

void SWARM_COVQ::DesalocarParams(NNS nns)
{
    free(params.mean_x);
    
    switch (nns)
    {
        case IENNS:
            free(params.sum_xf);
            free(params.sum_xs);
            break;
            
        case EENNS:
            free(params.variance_x);
            break;
            
        case IEENNS:
            free(params.variance_x);
            break;
            
        case EEENNS:
            free(params.variance_x);
            free(params.norm_x);
            break;
            
        default:
            break;
    }
}

void SWARM_COVQ::CalcularParamENNS(double **vector, double *mean, int num_vector)
{
    int j, n;
    double value;
    
    for(n = 0; n < num_vector; n++)
    {
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += vector[n][j];
        }
        
        mean[n] = value/K;
    }
}

void SWARM_COVQ::CalcularParamIENNS(double **vector, double *Mean, double *Sum_f, double *Sum_s, int num_vector)
{
    int j, n;
    double value, value_f, value_s;
    
    for(n = 0; n < num_vector; n++)
    {
        value_f = 0;
        value_s = 0;
        
        // soma da primeira parte
        for(j = 0; j < (K/2) - 1; j++)
        {
            value_f += vector[n][j];
        }
        
        Sum_f[n] = value_f;
        
        // soma da segunda parte
        for(j = K/2; j < K; j++)
        {
            value_s += vector[n][j];
        }
        
        Sum_s[n] = value_s;
        value = value_f + value_s;
        Mean[n] = value/K;
    }
}

void SWARM_COVQ::CalcularParamEENNS(double **vector, double *Mean, double *Variance, int num_vector)
{
    int j, n;
    double value;
    
    for(n = 0; n < num_vector; n++)
    {
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += vector[n][j];
        }
        
        Mean[n] = value/K;
        
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += (vector[n][j] - Mean[n]) * (vector[n][j] - Mean[n]);
        }
        
        Variance[n] = sqrt(value);
    }
}

void SWARM_COVQ::CalcularParamEEENNS(double **vector, double *Mean, double *Variance, double *Norm, int num_vector)
{
    int j, n;
    double value;
    
    for(n = 0; n < num_vector; n++)
    {
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += vector[n][j];
        }
        
        Mean[n] = value/K;
        
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += (vector[n][j] - Mean[n]) * (vector[n][j] - Mean[n]);
        }
        
        Variance[n] = sqrt(value);
        
        value = 0;
        
        for(j = 0; j < K; j++)
        {
            value += vector[n][j] * vector[n][j];
        }
        
        Norm[n] = sqrt(value);
    }
}

void SWARM_COVQ::Ordenar(int *Sort_Mean_c, int *Sort_Mean_inv, double *Mean_c)
{
    int i, j, n, index;
    
    Sort_Mean_c[0] = 0;
    
    for(n = 1; n < N; n++)
    {
        if(Mean_c[n] > Mean_c[Sort_Mean_c[n-1]])
        {
            Sort_Mean_c[n] = n;
        }
        else if (Mean_c[n] == Mean_c[Sort_Mean_c[n-1]])
        {
            index = n - 2;
            
            while (index >= 0 && Mean_c[n] == Mean_c[Sort_Mean_c[index]])
            {
                index--;
            }
            
            if (index >= 0)
                Sort_Mean_c[n] = index;
            else
                Sort_Mean_c[n] = n + 1;
        }
        else
        {
            for(j = 0; j < n; j++)
            {
                if(Mean_c[n] < Mean_c[Sort_Mean_c[j]])
                {
                    for(i = n; i > j; i--)
                    {
                        Sort_Mean_c[i] = Sort_Mean_c[i-1];
                    }

                    Sort_Mean_c[j] = n;
                    break;
                }
            }
        }
    }
    
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(Sort_Mean_c[i] == j)
            {
                Sort_Mean_inv[j] = i;
            }
        }
    }
}

void SWARM_COVQ::ClassificarMedias(int *Classifica_Dif, int *Sort_Mean_c, int *Sort_Mean_inv, double *mean_x, double *Mean_c)
{
    int i, n;
    double value, current_value;
    
    for(n = 0; n < nTreino; n++)
    {
        value = 0;
        current_value = 0;
        
        for(i = 0; i < N; i++)
        {
            value = (mean_x[n] - Mean_c[i]) * (mean_x[n] - Mean_c[i]);
            
            if(i == 0)
            {
                current_value = value;
                Classifica_Dif[n] = Sort_Mean_inv[0]; // classificação do vetor-treino
            }
            else if(value < current_value)
            {
                current_value = value;
                Classifica_Dif[n] = Sort_Mean_inv[i]; // classificação do vetor-treino
            }
        }
    }
}

#endif /* SWARM_COVQ_H */

//
//  swarm_vq.hpp
//  VQ
//
//  Created by Felipe Ferreira on 18/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef swarm_vq_h
#define swarm_vq_h

#include <iostream>
#include <cmath>
#include "../../Util/Util/general.h"
#include "../../Util/Util/rand.h"
#include "../../Util/Util/inout.h"
#include "../../Util/Util/definicoes.h"

using namespace std;

class SWARM_VQ
{
public:
    int qtd_part;
    int N;
    int K;
    int nTreino;
    double complexidade;
    SCALING scaleType;
    double scale;
    NNS nns;
    
    virtual double * Run(double **dic, double **treino) = 0;
    
    // Construtores
    SWARM_VQ();
    SWARM_VQ(int N, int K);
    SWARM_VQ(int N, int K, NNS nns);
    SWARM_VQ(int N, int K, int qtd_part, NNS nns);
    SWARM_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale);
    
    // Destrutores
    virtual ~SWARM_VQ();
    
protected:
    double **gbest;
    double dist_gbest;
    
    void AtualizarCentroides(double **treino, double **dic, int *particao);
    void AtualizarCentroides(double **treino, double **dic, int *particao, int interacao);
    double CalcularDistorcaoMedia(double **treino, double **dic, int *particao);
    double CalcularDEQ(double *vetor_treino, double *vetor_dic, double distAtual);
    
    void Classificar(double **treino, double **dic, int *particao);
    
    void CalcularParamTreino(double **vector);
    void DesalocarParams();
    
private:
    Params params;
    
    void ClassificarPorBT(double **treino, double **dic, int *particao);
    void ClassificarPorPDS(double **treino, double **dic, int *particao);
    void ClassificarPorENNS(double **treino, double **dic, int *particao);
    void ClassificarPorIENNS(double **treino, double **dic, int *particao);
    void ClassificarPorEENNS(double **treino, double **dic, int *particao);
    void ClassificarPorIEENNS(double **treino, double **dic, int *particao);
    void ClassificarPorEEENNS(double **treino, double **dic, int *particao);
    void ClassificarPorDTA(double **treino, double **dic, int *particao);
    void ClassificarPorIDTA(double **treino, double **dic, int *particao);
    
    void CalcularParamENNS(double **vector, double *mean, int num_vector);
    void CalcularParamIENNS(double **vector, double *Mean, double *Sum_f, double *Sum_s, int num_vector);
    void CalcularParamEENNS(double **vector, double *Mean, double *Variance, int num_vector);
    void CalcularParamEEENNS(double **vector, double *Mean, double *Variance, double *Norm, int num_vector);
    void CalcularParamDTA(double **vector, double *sum_squared, double *sum, double *max, int num_vector);
    void CalcularParamIDTA(double **vector, double *sum_squared, double *sum, double *max, double *mean, int num_vector);
    void Ordenar(int *Sort_Mean_c, int *Sort_Mean_inv, double *Mean_c);
    void ClassificarMedias(int *Classifica_Dif, int *Sort_Mean_c, int *Sort_Mean_inv, double *mean_x, double *Mean_c);
};

SWARM_VQ::SWARM_VQ()
{
    this->gbest = NULL;
}

SWARM_VQ::SWARM_VQ(int N, int K)
{
    this->N = N;
    this->K = K;
    this->qtd_part = 1;
    this->gbest = nullptr;
    this->complexidade = 0;
}

SWARM_VQ::SWARM_VQ(int N, int K, NNS nns)
{
    this->N = N;
    this->K = K;
    this->qtd_part = 1;
    this->gbest = nullptr;
    this->complexidade = 0;
    this->nns = nns;
}

SWARM_VQ::SWARM_VQ(int N, int K, int qtd_part, NNS nns)
{
    this->N = N;
    this->K = K;
    this->qtd_part = qtd_part;
    this->gbest = nullptr;
    this->complexidade = 0;
    this->nns = nns;
}

SWARM_VQ::SWARM_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale)
{
    this->N = N;
    this->K = K;
    this->qtd_part = qtd_part;
    this->gbest = NULL;
    this->scaleType = scaleType;
    this->scale = scale;
    this->complexidade = 0;
    this->nns = nns;
}

SWARM_VQ::~SWARM_VQ()
{
    if (gbest != nullptr)
        desaloca_matrizd(gbest, N);
}

void SWARM_VQ::AtualizarCentroides(double **treino, double **dic, int *particao)
{
    int qtd[N];
    int i, j;
    
    for(i = 0; i < N; i++)
        qtd[i] = 0;
    
    for(i = 0; i < nTreino; i++)
        qtd[particao[i]]++;
    
    for(i = 0; i < N; i++)
        if(qtd[i] > 0)
            for(j = 0; j < K; j++)
                dic[i][j] = 0.0;
    
    for(i = 0; i < nTreino; i++)
        for(j = 0; j < K; j++)
            dic[particao[i]][j] += treino[i][j];
    
    for(i = 0; i < N; i++)
        if(qtd[i] > 0)
            for(j = 0; j < K; j++)
                dic[i][j] = dic[i][j]/qtd[i];
}

void SWARM_VQ::AtualizarCentroides(double **treino, double **dic, int *particao, int iteracao)
{
    int qtd[N];
    int i, j;
    double aux[N][K];
    double s = 0;
    
    if(scaleType == FIXED)
    {
        s = scale;
    }
    else if(scaleType == VARIABLE)
    {
        s = 1 + (double)scale/(scale + iteracao);
    }
    
    for(i = 0; i < N; i++)
        qtd[i] = 0;
    
    for(i = 0; i < nTreino; i++)
        qtd[particao[i]]++;
    
    for(i = 0; i < N; i++)
        if(qtd[i] > 0)
            for(j = 0; j < K; j++)
                aux[i][j] = 0.0;
    
    for(i = 0; i < nTreino; i++)
        for(j = 0; j < K; j++)
            aux[particao[i]][j] += treino[i][j];
    
    for(i = 0; i < N; i++)
        if(qtd[i] > 0)
            for(j = 0; j < K; j++)
            {
                if(scaleType != NONE)
                {
                    dic[i][j] = dic[i][j] + s * ((aux[i][j]/qtd[i]) - dic[i][j]);
                }
                else
                    dic[i][j] = aux[i][j]/qtd[i];
            }
}


double SWARM_VQ::CalcularDistorcaoMedia(double **treino, double **dic, int *particao)
{
    double dist = 0.0;
    int i, j;
    
    for(j = 0; j < nTreino; j++)
    {
        for(i = 0; i < K; i++)
        {
            dist += (dic[particao[j]][i] - treino[j][i]) * (dic[particao[j]][i] - treino[j][i]);
        }
    }
    
    return dist/(K*nTreino);
}

// Calcular distância Euclidiana quadrática
double SWARM_VQ::CalcularDEQ(double *vetor_treino, double *vetor_dic, double distAtual)
{
    double distq;
    int i;
    
    distq = 0.0;
    
    for(i = 0; i < K; i++)
    {
        distq += (vetor_treino[i] - vetor_dic[i]) * (vetor_treino[i] - vetor_dic[i]);
        
        if (distq > distAtual)
            break;
    }
    
    if (i == K)
    {
        complexidade++;
    }
    
    return distq;
}

void SWARM_VQ::Classificar(double **treino, double **dic, int *particao)
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
            
        case DTA:
            ClassificarPorDTA(treino, dic, particao);
            break;
            
        case IDTA:
            ClassificarPorIDTA(treino, dic, particao);
            break;
            
        default:
            break;
    }
}

void SWARM_VQ::ClassificarPorBT(double **treino, double **dic, int *particao)
{
    double dist, aux;
    int i, j, c;
    
    for(i = 0; i < nTreino; i++)
    {
        dist = 1e20;
        
        for(j = 0; j < N; j++)
        {
            aux = 0.0;
            
            for(c = 0; c < K; c++)
            {
                aux += (dic[j][c] - treino[i][c]) * (dic[j][c] - treino[i][c]);
            }
            
            complexidade++;
            
            if(aux < dist)
            {
                dist = aux;
                particao[i] = j;
            }
        }
    }
}

void SWARM_VQ::ClassificarPorPDS(double **treino, double **dic, int *particao)
{
    double dist, aux;
    int i, j, c;
    
    for(i = 0; i < nTreino; i++)
    {
        dist = 1e20;
        
        for(j = 0; j < N; j++)
        {
            aux = 0.0;
            
            for(c = 0; c < K; c++)
            {
                aux += (dic[j][c] - treino[i][c]) * (dic[j][c] - treino[i][c]);
                
                if (aux > dist)
                {
                    break;
                }
            }
            
            if (c == K)
            {
                complexidade++;
            }
            
            if(aux < dist)
            {
                dist = aux;
                particao[i] = j;
            }
        }
    }
}

void SWARM_VQ::ClassificarPorENNS(double **treino, double **dic, int *particao)
{
    int n, count_up, count_down;
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
        
        dist = CalcularDEQ(treino[n], dic[sort_mean_c[classifica_dif[n]]], 1e20);
        particao[n] = sort_mean_c[classifica_dif[n]];
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(mean_c[sort_mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    dist_aux = CalcularDEQ(treino[n], dic[sort_mean_c[count_down]], dist);
                    
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
                    dist_aux = CalcularDEQ(treino[n], dic[sort_mean_c[count_up]], dist);
                    
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

void SWARM_VQ::ClassificarPorIENNS(double **treino, double **dic, int *particao)
{
    int n, count_up, count_down;
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

        dist = CalcularDEQ(treino[n], dic[Sort_Mean_c[Classifica_Dif[n]]], 1e20);
        particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        
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
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_down]], dist);
                        
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
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_up]], dist);
                        
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

void SWARM_VQ::ClassificarPorEENNS(double **treino, double **dic, int *particao)
{
    int n, count_up, count_down;
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
        
        dist = CalcularDEQ(treino[n], dic[Sort_Mean_c[Classifica_Dif[n]]], 1e20);
        particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        
        while(true)
        {
            if(count_down != 0)
            {
                count_down--;
                
                if(Mean_c[Sort_Mean_c[count_down]] > params.mean_x[n] - sqrt(dist/K))
                {
                    if(Variance_c[Sort_Mean_c[count_down]] > params.variance_x[n] - sqrt(dist))
                    {
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_down]], dist);
                        
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
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_up]], dist);
                        
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

void SWARM_VQ::ClassificarPorIEENNS(double **treino, double **dic, int *particao)
{
    int n, count_up, count_down;
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
        
        dist = CalcularDEQ(treino[n], dic[Sort_Mean_c[Classifica_Dif[n]]], 1e20);
        particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        
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
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_down]], dist);
                        
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
                        dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_up]], dist);
                        
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

void SWARM_VQ::ClassificarPorEEENNS(double **treino, double **dic, int *particao)
{
    int n, count_up, count_down;
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
        
        dist = CalcularDEQ(treino[n], dic[Sort_Mean_c[Classifica_Dif[n]]], 1e20);
        particao[n] = Sort_Mean_c[Classifica_Dif[n]];
        
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
                            dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_down]], dist);
                            
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
                            dist_aux = CalcularDEQ(treino[n], dic[Sort_Mean_c[count_up]], dist);
                            
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

void SWARM_VQ::ClassificarPorDTA(double **treino, double **dic, int *particao)
{
    int n, c, k;
    double dist, dist_aux;
    double sum_squared_c[N], sum_c[N], max_c[N];
    
    CalcularParamDTA(dic, sum_squared_c, sum_c, max_c, N);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 1e20;
        
        for (c = 0; c < N; c++)
        {
            if((params.sum_squared_x[n] + sum_squared_c[c] - params.max_x[n]*sum_c[c]) <= dist)
            {
                if((params.sum_squared_x[n] + sum_squared_c[c] - max_c[c]*params.sum_x[n]) <= dist)
                {
                    dist_aux = 0.0;
                    
                    for(k = 0; k < K; k++)
                    {
                        dist_aux += (dic[c][k] - treino[n][k]) * (dic[c][k] - treino[n][k]);
                        
                        if (dist_aux > dist)
                        {
                            break;
                        }
                    }
                    
                    if(dist_aux < dist)
                    {
                        dist = dist_aux;
                        particao[n] = c;
                    }
                }
            }
        }
    }
}

void SWARM_VQ::ClassificarPorIDTA(double **treino, double **dic, int *particao)
{
    int n, k, count_up, count_down;
    double dist, dist_aux;
    double sum_squared_c[N], sum_c[N], max_c[N];
    double mean_c[N];
    int sort_mean_c[N], sort_mean_inv[N], classifica_dif[nTreino];
    
    CalcularParamIDTA(dic, sum_squared_c, sum_c, max_c, mean_c, N);
    Ordenar(sort_mean_c, sort_mean_inv, mean_c);
    ClassificarMedias(classifica_dif, sort_mean_c, sort_mean_inv, params.sum_x, sum_c);
    
    for(n = 0; n < nTreino; n++)
    {
        dist = 0;
        dist_aux = 0;
        count_up = classifica_dif[n];
        count_down = classifica_dif[n];
        
        dist = CalcularDEQ(treino[n], dic[sort_mean_c[classifica_dif[n]]], 1e20);
        particao[n] = sort_mean_c[classifica_dif[n]];
        
        while(count_up >= 0 || count_down < N)
        {
            if(count_up >= 0)
            {
                if(((params.sum_squared_x[n] + sum_squared_c[sort_mean_c[count_up]] - params.max_x[n]*sum_c[sort_mean_c[count_up]]) <= dist) &&
                   ((params.sum_squared_x[n] + sum_squared_c[sort_mean_c[count_up]] - max_c[sort_mean_c[count_up]]*params.sum_x[n]) <= dist))
                {
                    dist_aux = 0.0;
                    
                    for(k = 0; k < K; k++)
                    {
                        dist_aux += (dic[sort_mean_c[count_up]][k] - treino[n][k]) * (dic[sort_mean_c[count_up]][k] - treino[n][k]);
                        
                        if (dist_aux > dist)
                        {
                            break;
                        }
                    }
                    
                    if(dist_aux < dist)
                    {
                        dist = dist_aux;
                        particao[n] = sort_mean_c[count_up];
                    }
                    
                    
                }

                count_up--;
            }
            
            if(count_down < N)
            {
                if(((params.sum_squared_x[n] + sum_squared_c[sort_mean_c[count_down]] - params.max_x[n]*sum_c[sort_mean_c[count_down]]) <= dist) &&
                   ((params.sum_squared_x[n] + sum_squared_c[sort_mean_c[count_down]] - max_c[sort_mean_c[count_down]]*params.sum_x[n]) <= dist))
                {
                    dist_aux = 0.0;
                    
                    for(k = 0; k < K; k++)
                    {
                        dist_aux += (dic[sort_mean_c[count_down]][k] - treino[n][k]) * (dic[sort_mean_c[count_down]][k] - treino[n][k]);
                        
                        if (dist_aux > dist)
                        {
                            break;
                        }
                    }
                    
                    if(dist_aux < dist)
                    {
                        dist = dist_aux;
                        particao[n] = sort_mean_c[count_down];
                    }
                    
                    
                }

                count_down++;
            }
        }
    }
}

void SWARM_VQ::CalcularParamTreino(double **vector)
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
            
        case DTA:
            params.sum_squared_x = aloca_vetor_double(nTreino);
            params.sum_x = aloca_vetor_double(nTreino);
            params.max_x = aloca_vetor_double(nTreino);
            
            CalcularParamDTA(vector, params.sum_squared_x, params.sum_x, params.max_x, nTreino);
            break;
            
        case IDTA:
            params.sum_squared_x = aloca_vetor_double(nTreino);
            params.sum_x = aloca_vetor_double(nTreino);
            params.max_x = aloca_vetor_double(nTreino);
            
            CalcularParamIDTA(vector, params.sum_squared_x, params.sum_x, params.max_x, params.mean_x, nTreino);
            break;
            
        default:
            free(params.mean_x);
            break;
    }
}

void SWARM_VQ::DesalocarParams()
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
            
        case DTA:
            free(params.sum_squared_x);
            free(params.sum_x);
            free(params.max_x);
            break;
            
        case IDTA:
            free(params.sum_squared_x);
            free(params.sum_x);
            free(params.max_x);
            break;
            
        default:
            break;
    }
}

void SWARM_VQ::CalcularParamENNS(double **vector, double *mean, int num_vector)
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

void SWARM_VQ::CalcularParamIENNS(double **vector, double *Mean, double *Sum_f, double *Sum_s, int num_vector)
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

void SWARM_VQ::CalcularParamEENNS(double **vector, double *Mean, double *Variance, int num_vector)
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

void SWARM_VQ::CalcularParamEEENNS(double **vector, double *Mean, double *Variance, double *Norm, int num_vector)
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

void SWARM_VQ::CalcularParamDTA(double **vector, double *sum_squared, double *sum, double *max, int num_vector)
{
    int j, n;
    double value, value_squared, value_max;
    
    for(n = 0; n < num_vector; n++)
    {
        value = 0;
        value_squared = 0;
        value_max = 0;
        
        for(j = 0; j < K; j++)
        {
            value_squared += vector[n][j] * vector[n][j];
            value += vector[n][j];
            
            if(vector[n][j] > value_max)
            {
                value_max = vector[n][j];
            }
        }
        
        sum_squared[n] = value_squared;
        sum[n] = value;
        max[n] = 2*value_max;
    }
}

void SWARM_VQ::CalcularParamIDTA(double **vector, double *sum_squared, double *sum, double *max, double *mean, int num_vector)
{
    int j, n;
    double value, value_squared, value_max;
    
    for(n = 0; n < num_vector; n++)
    {
        value = 0;
        value_squared = 0;
        value_max = 0;
        
        for(j = 0; j < K; j++)
        {
            value_squared += vector[n][j] * vector[n][j];
            value += vector[n][j];
            
            if(vector[n][j] > value_max)
            {
                value_max = vector[n][j];
            }
        }
        
        sum_squared[n] = value_squared;
        sum[n] = value;
        max[n] = 2*value_max;
        mean[n] = value / K;
    }
}

void SWARM_VQ::Ordenar(int *Sort_Mean_c, int *Sort_Mean_inv, double *Mean_c)
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

void SWARM_VQ::ClassificarMedias(int *Classifica_Dif, int *Sort_Mean_c, int *Sort_Mean_inv, double *mean_x, double *Mean_c)
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

#endif /* swarm_vq_h */

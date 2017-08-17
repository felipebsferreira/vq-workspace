//
//  PSO_COVQ.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 09/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _PSO_COVQ_H_
#define _PSO_COVQ_H_

#include "swarm_covq.hpp"

using namespace std;

class PSO_COVQ : public SWARM_COVQ
{
public:
    double c1 = 0.8;          // Constante de aceleração cognitiva
    double c2 = 0.2;          // Constante de aceleração social
    double w = 0.7;
    
    double * Run(double **dic, double **treino, NNS nns);
    
    PSO_COVQ();
    PSO_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
    
private:
    Particle * InitiateSwarm(double **dic);
    void DealocateSwarm(Particle *swarm);
    
    void PSOUpdate(Particle *swarm, double **treino, int niter);
    
    void UpdateVelocity(Particle *swarm);
    void UpdatePosition(Particle *swarm);
};

PSO_COVQ::PSO_COVQ() : SWARM_COVQ()
{
    
}

PSO_COVQ::PSO_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, qtd_part, erro, scaleType, scale)
{
    
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * PSO_COVQ::Run(double **dic, double **treino, NNS nns)
{
    int niter;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;
    Particle *swarm;
    
    result = aloca_vetor_double(4);
    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino, nns);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio PSO" << endl;
    
    inicio = clock();
    
    while(true)
    {
        cout << "Iteracao = " << niter << endl;
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            swarm[i].dist_atual = CalcularDistorcaoMedia(treino, swarm[i].dicionario, swarm[i].particao);
            
            if(swarm[i].dist_atual < dist_gbest)
            {
                dist_gbest = swarm[i].dist_atual;
                copy(gbest, swarm[i].dicionario, N, K);
            }
        }
        
        dist_atual = dist_gbest;
        reducao = ((dist_anterior - dist_atual)/dist_anterior);
        
        cout << "Reducao Percentual = " << reducao << endl;
        cout << "Distorcao = " << dist_gbest << endl;
        
        if(dist_atual == 0 || fabs(reducao) < limiar) break;
        
        dist_anterior = dist_atual;
        
        PSOUpdate(swarm, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            AtualizarCentroides(treino, swarm[i].dicionario, swarm[i].particao, niter);
        }
        
        niter++;
    }
    
    fim = clock();
    
    result[0] = niter;
    result[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    result[2] = dist_gbest;
    
    DealocateSwarm(swarm);
    
    if (nns != BT && nns != PDS)
        DesalocarParams(nns);
    
    copy(dic, gbest, N, K);
    
    return result;
}

void PSO_COVQ::PSOUpdate(Particle *swarm, double **treino, int niter)
{
    int i;
    
    for (i = 0; i < qtd_part; i++)
    {
        swarm[i].dist_atual = CalcularDistorcaoMedia(treino, swarm[i].dicionario, swarm[i].particao);
        
        // Novo pbest?
        if (swarm[i].dist_atual < swarm[i].dist_pbest)
        {
            copy(swarm[i].pbest, swarm[i].dicionario, N, K);
            swarm[i].dist_pbest = swarm[i].dist_atual;
        }
        
        // Novo gbest?
        if (swarm[i].dist_atual < dist_gbest)
        {
            copy(gbest, swarm[i].dicionario, N, K);
            dist_gbest = swarm[i].dist_atual;
        }
    }
    
    UpdateVelocity(swarm);
    UpdatePosition(swarm);
}

Particle * PSO_COVQ::InitiateSwarm(double **dic)
{
    Particle *swarm = new Particle[qtd_part];
    
    for(int i = 0; i < qtd_part; i++)
    {
        swarm[i].dicionario = aloca_matrizd(N, K);
        swarm[i].pbest = aloca_matrizd(N, K);
        swarm[i].velocidade = aloca_matrizd(N, K);
        swarm[i].particao = aloca_vetor_int(nTreino);
        swarm[i].dist_pbest = 1e20;
        
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < K; k++)
            {
                swarm[i].dicionario[j][k] = dic[(i * N) + j][k];
            }
        }
    }
    
    gbest = aloca_matrizd(N, K);
    dist_gbest = 1e20;
    
    return swarm;
}

void PSO_COVQ::DealocateSwarm(Particle *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].dicionario, N);
        desaloca_matrizd(swarm[i].pbest, N);
        free(swarm[i].particao);
    }
    
    delete swarm;
}

void PSO_COVQ::UpdateVelocity(Particle *swarm)
{
    int i, j, k;
    double r1, r2;
    
    for (i = 0; i < qtd_part; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < K; k++)
            {
                r1 = unif();
                r2 = unif();
                
                swarm[i].velocidade[j][k] = w * swarm[i].velocidade[j][k]
                    + c1*r1*(swarm[i].pbest[j][k] - swarm[i].dicionario[j][k])
                    + c2*r2*(gbest[j][k] - swarm[i].dicionario[j][k]);
            }
        }
    }
}

void PSO_COVQ::UpdatePosition(Particle *swarm)
{
    int i, j, k;
    
    for (i = 0; i < qtd_part; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < K; k++)
            {
                swarm[i].dicionario[j][k] = swarm[i].dicionario[j][k] + swarm[i].velocidade[j][k];
                
                if(swarm[i].dicionario[j][k] > 255)
                {
                    swarm[i].dicionario[j][k] = 255;
                }
                else if(swarm[i].dicionario[j][k] < 0)
                {
                    swarm[i].dicionario[j][k] = 0;
                }
            }
        }
    }
}

#endif /* _PSO_COVQ_H_ */

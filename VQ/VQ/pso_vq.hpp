//
//  pso_vq.hpp
//  VQ
//
//  Created by Felipe Ferreira on 21/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef pso_vq_h
#define pso_vq_h

#include "swarm_vq.hpp"

using namespace std;

class PSO_VQ : public SWARM_VQ
{
public:
    double c1 = 0.1; //0.8          // Constante de aceleração cognitiva
    double c2 = 0.2; //0.2         // Constante de aceleração social
    double w = 0.7;           // Fator de inércia
    
    double * Run(double **dic, double **treino);
    
    PSO_VQ();
    PSO_VQ(int N, int K, int qtd_part, NNS nns);
    PSO_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale);
    
private:
    Particle * InitiateSwarm(double **dic);
    void DealocateSwarm(Particle *swarm);
    
    void PSOUpdate(Particle *swarm, double **treino, int niter);
    
    void UpdateVelocity(Particle *swarm);
    void UpdatePosition(Particle *swarm);
};

PSO_VQ::PSO_VQ() : SWARM_VQ()
{
    
}

PSO_VQ::PSO_VQ(int N, int K, int qtd_part, NNS nns) : SWARM_VQ(N, K, qtd_part, nns)
{
    
}

PSO_VQ::PSO_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale) : SWARM_VQ(N, K, qtd_part, nns, scaleType, scale)
{
    
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * PSO_VQ::Run(double **dic, double **treino)
{
    int niter;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;
    Particle *swarm;
    
    result = aloca_vetor_double(4);
    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio PSO" << endl;
    
    inicio = clock();
    
    while(true)
    {
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao);
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
        
        if(dist_atual == 0 || fabs(reducao) < limiar)
            break;
        
        dist_anterior = dist_atual;
        
        PSOUpdate(swarm, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao);
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
        DesalocarParams();
    
    copy(dic, gbest, N, K);
    
    return result;
}

void PSO_VQ::PSOUpdate(Particle *swarm, double **treino, int niter)
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

Particle * PSO_VQ::InitiateSwarm(double **dic)
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

void PSO_VQ::DealocateSwarm(Particle *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].dicionario, N);
        desaloca_matrizd(swarm[i].pbest, N);
        free(swarm[i].particao);
    }
    
    delete swarm;
}

void PSO_VQ::UpdateVelocity(Particle *swarm)
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

void PSO_VQ::UpdatePosition(Particle *swarm)
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

#endif /* pso_vq_h */

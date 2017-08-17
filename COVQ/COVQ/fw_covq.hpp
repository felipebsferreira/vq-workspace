//
//  FW_COVQ.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 09/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _FW_COVQ_H_
#define _FW_COVQ_H_

#include "swarm_covq.hpp"
#include <math.h>

using namespace std;

class FW_COVQ : public SWARM_COVQ
{
public:
    int m;              // total number of regular sparks for all fireworks
    int mHat;           // number of Gaussian sparks
    double a;           // controls min sparks per firework
    double b;           // controls max sparks per firework
    int A;              // max amplitude
    double minX;
    double maxX;
    
    double * Run(double **dic, double **treino, NNS nns);
    
    FW_COVQ();
    FW_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
    
private:
    Firework * InitiateSwarm(double **dic);
    void DealocateSwarm(Firework *swarm);
    
    void FWUpdate(Firework *swarm, double **treino, int niter);
    int NumberOfSparks(Firework *fireworks, int index);
    double YMax(Firework *fireworks);
    double YMin(Firework *fireworks);
};

FW_COVQ::FW_COVQ() : SWARM_COVQ()
{

}

FW_COVQ::FW_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, qtd_part, erro, scaleType, scale)
{
    m = qtd_part * 10;
    mHat = 5;             
    a = 0.04;
    b = 0.8; 
    A = 40;  
    minX = 0;
    maxX = 255;
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * FW_COVQ::Run(double **dic, double **treino, NNS nns)
{
    int niter;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;
    Firework *swarm;

    result = aloca_vetor_double(4);
    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino, nns);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio FW" << endl;
    
    inicio = clock();
    
    while(true)
    {
        cout << "Iteracao = " << niter << endl;
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].position, swarm[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            swarm[i].error = CalcularDistorcaoMedia(treino, swarm[i].position, swarm[i].particao);
            
            if(swarm[i].error < dist_gbest)
            {
                dist_gbest = swarm[i].error;
                copy(gbest, swarm[i].position, N, K);
            }
        }
        
        dist_atual = dist_gbest;
        reducao = ((dist_anterior - dist_atual)/dist_anterior);
        
        cout << "Reducao Percentual = " << reducao << endl;
        cout << "Distorcao = " << dist_gbest << endl;
        
        if(dist_atual == 0 || fabs(reducao) < limiar) break;
        
        dist_anterior = dist_atual;
        
        FWUpdate(swarm, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].position, swarm[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            AtualizarCentroides(treino, swarm[i].position, swarm[i].particao, niter);
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

void FW_COVQ::FWUpdate(Firework *swarm, double **treino, int niter)
{
    int i;
    
    //swarm[i].error = CalcularDistorcaoMedia(treino, swarm[i].position, swarm[i].particao);

    for (i = 0; i < qtd_part; i++)
    {
          
        
        
    }
}

int FW_COVQ::NumberOfSparks(Firework *fireworks, int index)
{
    int minSparks = round(a * m);
    
    if (minSparks < 1) 
        minSparks = 1;
  
    int maxSparks = round(b * m);
  
    if (maxSparks > m - (qtd_part - 1) * minSparks)
        maxSparks = m - (qtd_part - 1) * minSparks;

    double yMax = YMax(fireworks);
    double sumDeltas = 0.0; // sum diffs between yMax and each error
    
    for (int i = 0; i < qtd_part; ++i)
        sumDeltas += yMax - fireworks[i].error;

  int[] numSparks = new int[n]; // the result
  for (int i = 0; i < n; ++i)
  {
    numSparks[i] = (int)Math.Round(m * (yMax - fireworks[i].error + 1.0E-10) / (sumDeltas + 1.0E-10));
    if (numSparks[i] < minSparks)
      numSparks[i] = minSparks;
    else if (numSparks[i] > maxSparks)
      numSparks[i] = maxSparks;
  }
  return numSparks;
}

double FW_COVQ::YMax(Firework *fireworks)
{
    // largest (worst) error in any firework
    double result = fireworks[0].error;
    
    for (int i = 1; i < qtd_part; i++)
    {
        if (fireworks[i].error > result)
            result = fireworks[i].error;    
    }
    
    return result;
}

double FW_COVQ::YMin(Firework *fireworks)
{
    // smallest (best) error in any firework
    double result = fireworks[0].error;
    
    for (int i = 1; i < qtd_part; i++)
    {
        if (fireworks[i].error < result)
            result = fireworks[i].error;    
    }
    
    return result;
}

Firework * FW_COVQ::InitiateSwarm(double **dic)
{
    Firework *swarm = new Firework[qtd_part];
    
    for(int i = 0; i < qtd_part; i++)
    {
        swarm[i].position = aloca_matrizd(N, K);
        swarm[i].error = 1e20;
        
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < K; k++)
            {
                swarm[i].position[j][k] = dic[(i * N) + j][k];
            }
        }
    }
    
    gbest = aloca_matrizd(N, K);
    dist_gbest = 1e20;
    
    return swarm;
}

void FW_COVQ::DealocateSwarm(Firework *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].position, N);
        free(swarm[i].particao);
    }
    
    delete swarm;
}

#endif /* _FW_COVQ_H_ */

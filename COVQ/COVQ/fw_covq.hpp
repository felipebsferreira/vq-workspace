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
#include "../../Util/Util/general.h"
#include "../../Util/Util/rand.h"
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
    NNS _nns;
    
    double * Run(double **dic, double **treino, NNS nns);
    
    FW_COVQ();
    FW_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
    
private:
    Firework * InitiateSwarm(double **dic);
    void DealocateSwarm(Firework *swarm);
    
    void FWUpdate(Firework *swarm, double **treino, int niter);
    int * NumberOfSparks(Firework *fireworks);
    double * Amplitudes(Firework * fireworks, int epoch, int maxEpochs);
    double Error(Firework * firework, double ** treino);
    int ** PickDimensions(int z);
    double MinAmplitude(int epoch, int maxEpochs);
    double YMax(Firework *fireworks);
    double YMin(Firework *fireworks);
    Firework * AlocateSparks(int qty);
    void DealocateSparks(Firework *sparks, int qty_sparks);
    void AddGaussianSparks(Firework *fireworks, double **treino);
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
    _nns = nns;
    
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

void FW_COVQ::FWUpdate(Firework *fireworks, double **treino, int niter)
{
    int maxEpochs = 30; // <- Criar mecanismo para não depender da qtd máxima de iterações

    int * numberSparks = NumberOfSparks(fireworks);
    double * amplitudes = Amplitudes(fireworks, niter, maxEpochs);
    
    for(int i = 0; i < qtd_part; i++)
    {
        double amp = amplitudes[i];

        if(fireworks[i].qty_sparks > 0)
        {
            DealocateSparks(fireworks[i].sparks, m + mHat);
        }
        
        fireworks[i].qty_sparks = numberSparks[i];
        
        fireworks[i].sparks = AlocateSparks(m + mHat);
        
        for(int j = 0; j < fireworks[i].qty_sparks; j++)
        {
            Firework *spark = new Firework;
            spark->position = aloca_matrizd(N, K);
            spark->particao = aloca_vetor_int(nTreino);
            spark->error = 1e20;
            spark->qty_sparks = 0;

            copy(spark->position, fireworks[i].position, N, K);
            
            int ** dimensions = PickDimensions((int)round(rand() % (N*K)));
            
            for(int n = 0; n < N; n++)
            {
                for(int k = 0; k < K; k++)
                {
                    if(dimensions[n][k] == 1)
                    {
                        spark->position[n][k] += amp * unif_neg();

                        if(spark->position[n][k] < minX || spark->position[n][k] > maxX)
                        {
                            spark->position[n][k] = (maxX - minX) * unif() + minX;
                        }
                    }
                }
            }
            
            spark->error = Error(spark, treino);
            
            if(spark->error < dist_gbest)
            {
                dist_gbest = spark->error;
                copy(gbest, spark->position, N, K);
            }

            fireworks[i].sparks[j] = *spark;

            desaloca_matriz_int(dimensions, N);
        }
        
        AddGaussianSparks(fireworks, treino);
        
        double **bestSparkPos = aloca_matrizd(N, K);
        double bestSparkErr = 1e20; 

        double **worstSparkPos = aloca_matrizd(N, K);
        double worstSparkErr = 0;
        
        for (int i = 0; i < qtd_part; i++) // numFireworks
        {
            for (int j = 0; j < fireworks[i].qty_sparks; j++) // number sparks in each firework
            {
                if (fireworks[i].sparks[j].error < bestSparkErr)
                {
                    bestSparkErr = fireworks[i].sparks[j].error;
                    copy(bestSparkPos, fireworks[i].sparks[j].position, N, K);
                }

                if (fireworks[i].sparks[j].error > worstSparkErr)
                {
                    worstSparkErr = fireworks[i].sparks[j].error;
                    copy(worstSparkPos, fireworks[i].sparks[j].position, N, K);
                }
            } // each spark
        } // each firework
        
        copy(fireworks[0].position, bestSparkPos, N, K);
        fireworks[0].error = bestSparkErr;

        copy(fireworks[1].position, worstSparkPos, N, K);
        fireworks[1].error = worstSparkErr;
        
        for (int i = 2; i < qtd_part; i++) // n-2 random sparks 
        {
            int row = (int)round(rand() % qtd_part);
            
            if(fireworks[row].qty_sparks > 0)
            {
                int col = (int)round(rand() % fireworks[row].qty_sparks);

                copy(fireworks[i].position, fireworks[row].sparks[col].position, N, K);
                fireworks[i].error = fireworks[row].sparks[col].error;
            }
        }
    }
}

void FW_COVQ::AddGaussianSparks(Firework *fireworks, double **treino)
{
    for(int g = 0; g < mHat; g++)
    {
        Firework *gSpark = new Firework;

        gSpark->position = aloca_matrizd(N, K);
        gSpark->particao = aloca_vetor_int(nTreino);
        gSpark->error = 1e20;
        gSpark->qty_sparks = 0;

        int i = (int)round(rand() % qtd_part);

        copy(gSpark->position, fireworks[i].position, N, K);

        int ** dimensions = PickDimensions((int)round(rand() % (N*K)));

        double e = cos(2.0 * M_PI * unif()) * sqrt(-2.0 * log(unif()));

        for(int n = 0; n < N; n++)
        {
            for(int k = 0; k < K; k++)
            {
                if(dimensions[n][k] == 1)
                {
                    gSpark->position[n][k] = gSpark->position[n][k] 
                        + (gbest[n][k] - gSpark->position[n][k]) * e; 

                    if(gSpark->position[n][k] < minX || gSpark->position[n][k] > maxX)
                    {
                        gSpark->position[n][k] = (maxX - minX) * unif() + minX;
                    }
                }
            }
        }
        
        gSpark->error = Error(gSpark, treino);
        
        if(gSpark->error < dist_gbest)
        {
            dist_gbest = gSpark->error;
            copy(gbest, gSpark->position, N, K);
        }
        
        fireworks[i].sparks[fireworks[i].qty_sparks] = *gSpark;
        fireworks[i].qty_sparks++;
        
        desaloca_matriz_int(dimensions, N);
    }
}

int * FW_COVQ::NumberOfSparks(Firework *fireworks)
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

    int * numSparks = aloca_vetor_int(qtd_part);
    
    for (int i = 0; i < qtd_part; i++)
    {
        numSparks[i] = (int)round(m * (yMax - fireworks[i].error + 1.0E-10) / (sumDeltas + 1.0E-10));

        if(numSparks[i] < minSparks)
            numSparks[i] = minSparks;

        if(numSparks[i] > maxSparks)
            numSparks[i] = maxSparks;
    }

    return numSparks;
}

double * FW_COVQ::Amplitudes(Firework * fireworks, int epoch, int maxEpochs)
{
    double yMin = YMin(fireworks);
    double sumDeltas = 0.0; // sum  diffs between yMin and each error
    
    for (int i = 0; i < qtd_part; i++)
        sumDeltas += fireworks[i].error - yMin;

    double * result = aloca_vetor_double(qtd_part); // an amplitude for each firework

    double minAmplitude = MinAmplitude(epoch, maxEpochs);

    for (int i = 0; i < qtd_part; i++)
    {
        result[i] = A * (fireworks[i].error - yMin + 1.0E-10) / (sumDeltas + 1.0E-10);
        if (result[i] < minAmplitude)
            result[i] = minAmplitude;
    }
    
    return result;
}

double FW_COVQ::Error(Firework * firework, double ** treino)
{   
    Classificar(treino, firework->position, firework->particao, _nns);
    return CalcularDistorcaoMedia(treino, firework->position, firework->particao);
}


int ** FW_COVQ::PickDimensions(int z)
{
    // pick z random dimensions of a position
    int ** result = aloca_matriz_int(N, K);
    int n, k;

    for(int i = 0; i < z; i++)
    {
        n = rand() % N;
        k = rand() % K;

        result[n][k] = 1;
    }

    return result;
}

double FW_COVQ::MinAmplitude(int epoch, int maxEpochs)
{
    // minimum amplitude for any firework at curr epoch 
    double Ainit = (maxX - minX) * 0.02;
    double Afinal = (maxX - minX) * 0.001;
    return Ainit - (Ainit - Afinal) * (sqrt(2 * maxEpochs - epoch) * epoch) / maxEpochs;
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
        swarm[i].sparks = AlocateSparks(m + mHat);
        swarm[i].qty_sparks = 0;
        swarm[i].particao = aloca_vetor_int(nTreino);
        
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

Firework * FW_COVQ::AlocateSparks(int qty)
{
    Firework *sparks = new Firework[qty];
    
    for(int i = 0; i < qty; i++)
    {
        sparks[i].position = aloca_matrizd(N, K);
        sparks[i].particao = aloca_vetor_int(nTreino);
        sparks[i].error = 1e20;
        sparks[i].qty_sparks = 0;
    }
    
    return sparks;
}

void FW_COVQ::DealocateSparks(Firework *sparks, int qty_sparks)
{
    for(int i = 0; i < qty_sparks; i++)
    {
        desaloca_matrizd(sparks[i].position, N);
        free(sparks[i].particao);
    }
    
    delete sparks;
}

#endif /* _FW_COVQ_H_ */

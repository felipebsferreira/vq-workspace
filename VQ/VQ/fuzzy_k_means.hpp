//
//  fuzzy_k_means.hpp
//  VQ
//
//  A fuzzy vector quantization approach to image compression
//  George E. Tsekouras
//
//  Created by Felipe Ferreira on 7/18/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef fuzzy_k_means_h
#define fuzzy_k_means_h

#include "../../Util/Util/alocacao.h"
#include "../../Util/Util/general_vq.h"

#include "swarm_vq.hpp"

class FuzzyKMeans : public SWARM_VQ
{
public:
    double m;
    
    double * Run(double **centroids, double **treino);
    
    FuzzyKMeans();
    FuzzyKMeans(int N, int K, NNS nns);
    
private:
    double UpdateMembership(int index, int cardinality, int *hypersphere, double *treino, double **dicionario);
    void UpdateMembership(double **membership, double **centroids, double **treino);
    void UpdateCentroids(double **membership, double **centroids, double **treino);
    int UpdateHypersphere(int cardinality, int *hypersphere, double *treino, double **dicionario, double **membership, int index);
    double ComputeDistortion(double **membership, double **centroids, double **treino);
    double MembershipFunction(double membership);
    void InitMembership(double **membership, double **treino);
};

FuzzyKMeans::FuzzyKMeans() : SWARM_VQ()
{
    
}

FuzzyKMeans::FuzzyKMeans(int N, int K, NNS nns) : SWARM_VQ(N, K, nns)
{

}

double * FuzzyKMeans::Run(double **dicionario, double **treino)
{
    double **membership;
    double *resultado;
    double currentDistortion, previousDistortion;
    double fuzzyThreshold, crispThreshold;
    double coeff, numerator, denominator, reducao;
    int **hypersphere;
    int *cardinality, *particao, *hypersphere_aux;
    int niter, cardinality_aux;
    bool turnAllCrisp;
    clock_t inicio, fim;
    
    resultado = aloca_vetor_double(4);
    membership = aloca_matriz_double(N, nTreino);
    hypersphere = aloca_matriz_int(nTreino, N);
    cardinality = aloca_vetor_int(nTreino);
    particao = aloca_vetor_int(nTreino);
    hypersphere_aux = aloca_vetor_int(N);
    
    CalcularParamTreino(treino);
    
    turnAllCrisp = false;
    previousDistortion = 1e20;
    fuzzyThreshold = 0.1;
    crispThreshold = 0.001;
    niter = 1;
    
    cout << "Inicio Fuzzy K-Means" << endl;

    inicio = clock();
    
    for (int i = 0; i < nTreino; i++)
    {
        for (int j = 0 ; j < N; j++)
        {
            hypersphere[i][j] = j;
        }
        
        cardinality[i] = N;
    }
    
    while (true)
    {
        Classificar(treino, dicionario, particao);
        currentDistortion = CalcularDistorcaoMedia(treino, dicionario, particao);
        
        //TODO: If (|Dold - D|) / Dold < e' Then index = 1
        
        reducao = fabs(previousDistortion - currentDistortion) / previousDistortion;
        
        cout << "Iteracao " << niter << " / " << turnAllCrisp << endl;
        cout << reducao << endl;
        
        if (reducao < fuzzyThreshold)
            turnAllCrisp = true;
        
        //TODO: If (|Dold - D|) / Dold < e Then Stop
        
        if (reducao < crispThreshold)
            break;
        
        previousDistortion = currentDistortion;
        
        if(!turnAllCrisp)
        {
            for (int k = 0; k < nTreino; k++)
            {
                if(cardinality[k] != 1)
                {
                    cardinality_aux = 0;
                    
                    for (int i = 0; i < cardinality[k]; i++)
                    {
                        membership[hypersphere[k][i]][k] = UpdateMembership(i, cardinality[k], hypersphere[k], treino[k], dicionario);
                        
                        if(membership[hypersphere[k][i]][k] < 0)
                            membership[hypersphere[k][i]][k] = 0;
                        
                        if(membership[hypersphere[k][i]][k] > 0)
                        {
                            hypersphere_aux[cardinality_aux] = hypersphere[k][i];
                            cardinality_aux++;
                        }
                    }
                    
                    cardinality[k] = cardinality_aux;
                    
                    for (int i = 0; i < cardinality_aux; i++)
                    {
                        hypersphere[k][i] = hypersphere_aux[i];
                    }
                    
                    cardinality[k] = UpdateHypersphere(cardinality[k], hypersphere[k], treino[k], dicionario, membership, k);
                }
                
                if(cardinality[k] <= 1)
                {
                    for (int i = 0; i < N; i++)
                    {
                        if(i == particao[k])
                        {
                            membership[i][k] = 1;
                        }
                        else
                            membership[i][k] = 0;
                    }
                }
            }
        }
        else
        {
            for (int k = 0; k < nTreino; k++)
            {
                for (int i = 0; i < N; i++)
                {
                    if(i == particao[k])
                    {
                        membership[i][k] = 1;
                    }
                    else
                        membership[i][k] = 0;
                }
            }
        }
        
        //TODO: Calculate the normalized membership degrees wik (1 <= i <= c) using Eq. (18)
        
        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < nTreino; k++)
            {
                coeff = 1e-20;
                
                for (int j = 0; j < N; j++)
                {
                    coeff += membership[j][k];
                }
                
                membership[i][k] = membership[i][k] / coeff;
            }
        }
        
        //TODO: Update the codebook vectors vi (1 <= i <= c) using Eq. (19)
        
        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < K; k++)
            {
                numerator = 0;
                denominator = 1e-20;
                
                for (int j = 0; j < nTreino; j++)
                {
                    numerator += MembershipFunction(membership[i][j]) * treino[j][k];
                    denominator += MembershipFunction(membership[i][j]);
                }
                
                dicionario[i][k] = numerator / denominator;
            }
        }
        
        niter++;
    }
    
    fim = clock();
    
    resultado[0] = niter;
    resultado[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    resultado[2] = currentDistortion;
    
    desaloca_matrizd(membership, N);
    desaloca_matriz_int(hypersphere, N);
    free(cardinality);
    free(particao);
    
    if (nns != BT && nns != PDS)
        DesalocarParams();
    
    return resultado;
}

void FuzzyKMeans::UpdateCentroids(double **membership, double **centroids, double **treino)
{
    double centroidMembershipSum[K];
    double membershipSum;
    
    for(int j = 0; j < N; j++)
    {
        membershipSum = 0.0;
        
        for(int k = 0; k < K; k++)
        {
            centroidMembershipSum[k] = 0;
        }
        
        for(int i = 0; i < nTreino; i++)
        {
            if(membership[j][i] != 0)
            {
                membershipSum += pow(membership[j][i], (double)m);
                
                for(int k = 0; k < K; k++)
                {
                    centroidMembershipSum[k] += treino[i][k] * pow(membership[j][i], (double)m);
                }
            }
        }
        
        for(int k = 0; k < K; k++)
        {
            centroids[j][k] = centroidMembershipSum[k] / membershipSum;
        }
    }
}

double FuzzyKMeans::UpdateMembership(int index, int cardinality, int *hypersphere, double *treino, double **dicionario)
{
    double coeff, dist_i, dist_j;
    int j, k;
    
    dist_i = 0;
    
    for(k = 0; k < K; k++)
    {
        dist_i += (treino[k] - dicionario[hypersphere[index]][k]) * (treino[k] - dicionario[hypersphere[index]][k]);
    }
    
    coeff = 0;
    
    for (j = 0; j < cardinality; j++)
    {
        dist_j = 0;
        
        for(k = 0; k < K; k++)
        {
            dist_j += (treino[k] - dicionario[hypersphere[j]][k]) * (treino[k] - dicionario[hypersphere[j]][k]);
        }
        
        coeff += (dist_i / dist_j) * (dist_i / dist_j);
    }
    
    return ((cardinality + 2) / 2) * (1 / coeff) - 0.5;
}

void FuzzyKMeans::UpdateMembership(double **membership, double **centroids, double **treino)
{
    double coeff, aux, dist;
    double M = 1.0/(m - 1);
    int i, j, k, c;
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < nTreino; j++)
        {
            dist = 0;
            
            for(k = 0; k < K; k++)
            {
                dist += (treino[j][k] - centroids[i][k]) * (treino[j][k] - centroids[i][k]);
            }
            
            coeff = 0.0;
            
            for (c = 0; c < N; c++)
            {
                aux = 0;
                
                for(k = 0; k < K; k++)
                {
                    aux += (treino[j][k] - centroids[c][k]) * (treino[j][k] - centroids[c][k]);
                }
                
                if (aux != 0)
                {
                    coeff += dist / aux;
                }
            }
            
            membership[i][j] = coeff == 0 ? 0 : 1 / pow(coeff, M);
        }
    }
}

int FuzzyKMeans::UpdateHypersphere(int cardinality, int *hypersphere, double *treino, double **dicionario, double **membership, int index)
{
    double coeff, dist_i, dist_j;
    int hypersphereOld[cardinality];
    int i, j, k, cardinalityNew;
    
    cardinalityNew = 0;
    
    for(i = 0; i < cardinality; i++)
    {
        hypersphereOld[i] = hypersphere[i];
    }
    
    for(i = 0; i < cardinality; i++)
    {
        dist_i = 0;
        
        for(k = 0; k < K; k++)
        {
            dist_i += (treino[k] - dicionario[hypersphereOld[i]][k]) * (treino[k] - dicionario[hypersphereOld[i]][k]);
        }
        
        coeff = 0;
        
        for (j = 0; j < cardinality; j++)
        {
            dist_j = 0;
            
            for(k = 0; k < K; k++)
            {
                dist_j += (treino[k] - dicionario[hypersphereOld[j]][k]) * (treino[k] - dicionario[hypersphereOld[j]][k]);
            }
            
            coeff += (1 / dist_j) * (1 / dist_j);
        }
        
        coeff = (cardinality + 2) / coeff;
        
        if(dist_i < coeff)
        {
            hypersphere[cardinalityNew] = hypersphereOld[i];
            cardinalityNew++;
        }
        else
        {
            membership[hypersphereOld[i]][index] = 0;
        }
    }
    
    return cardinalityNew;
}

double FuzzyKMeans::MembershipFunction(double membership)
{
    return membership / 2 + (membership * membership) / 2;
}

double FuzzyKMeans::ComputeDistortion(double **membership, double **centroids, double **treino)
{
    double dist = 0.0, aux;
    
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < nTreino; i++)
        {
            if(membership[j][i] != 0)
            {
                aux = 0;
                
                for(int k = 0; k < K; k++)
                {
                    aux += (treino[i][k] - centroids[j][k]) * (treino[i][k] - centroids[j][k]);
                }
                
                dist += pow(membership[j][i], (double)m) * aux;
            }
        }
    }
    
    return dist/(K*nTreino);
}

#endif /* fuzzy_k_means_h */

//
//  vq.hpp
//  VQ
//
//  Created by Felipe Ferreira on 09/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef VQ_H
#define VQ_H

#include "swarm_vq.hpp"

using namespace std;

class VQ : public SWARM_VQ
{
public:
    double * Run(double **dicionario, double **treino);

    VQ();
    VQ(int N, int K);
    VQ(int N, int K, NNS nns, SCALING scaleType, double scale);
};

VQ::VQ() : SWARM_VQ()
{

}

VQ::VQ(int N, int K) : SWARM_VQ(N, K)
{

}

VQ::VQ(int N, int K, NNS nns, SCALING scaleType, double scale) : SWARM_VQ(N, K, 1, nns, scaleType, scale)
{

}

// Executar algoritmo LBG
double * VQ::Run(double **dicionario, double **treino)
{
    int *particao;
    int niter;
    double limiar = 0.001;
    double distAtual, distAnterior, reducao;
    double *resultado;
    clock_t inicio, fim;

    particao = aloca_vetor_int(nTreino);
    resultado = aloca_vetor_double(5);
    CalcularParamTreino(treino);

    distAnterior = 1e20;
    niter = 1;

    cout << "Inicio QV" << endl;

    inicio = clock();

    while(true)
    {
        Classificar(treino, dicionario, particao);
        distAtual = CalcularDistorcaoMedia(treino, dicionario, particao);

        reducao = ((distAnterior - distAtual)/distAnterior);
        cout << reducao << endl;
        if(distAtual == 0 || fabs(reducao) < limiar) break;
        distAnterior = distAtual;

        AtualizarCentroides(treino, dicionario, particao, niter);

        niter++;
    }

    fim = clock();

    resultado[0] = niter;
    resultado[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    resultado[2] = distAtual;
    resultado[4] = complexidade/(nTreino * niter);

    free(particao);

    if (nns != BT && nns != PDS)
        DesalocarParams();

    return resultado;
}

#endif /* VQ_H */

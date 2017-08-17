//
//  rand.h
//  Util
//
//  Created by Felipe Ferreira on 30/11/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _RAND_H_
#define _RAND_H_

#include <stdlib.h>
#include <time.h>
#include <math.h>

void reset_rand(void);
double unif(void);
double gauss(float, float);
int gera_bits(float);

void reset_rand(void)
{
    /* reseta o gerador de numeros aleatorios */
    
    int stime;
    long ltime;
    ltime = time(NULL);
    stime =  (unsigned)ltime/2;
    srand(stime);
}

double unif(void)
{
    /* Gera uma sequencia de numeros aleatorios uniformemente distribuidos
     entre 0 e 1 */
    double x;
    
    x = ((double)rand()/(double)RAND_MAX);
    
    return (x);
}

double unif_neg(void)
{
    /* Gera uma sequencia de numeros aleatorios uniformemente distribuidos
     entre -1 e 1 */
    double x, signal;
    
    x = ((double)rand()/(double)RAND_MAX);
    signal = ((double)rand()/(double)RAND_MAX) < 0.5 ? -1 : 1;
    
    return (signal*x);
}

double gauss(float mean, float var)
{
    // Gera uma distribuicao gaussiana com media zero e variancia var
    // (Metodo Polar)
    
    double v1,v2,r;
    do{
        v1 = 2.0*unif() - 1.0;
        v2 = 2.0*unif() - 1.0;
        r = v1*v1 + v2*v2;
    } while(r>=1.0 || r == 0);
    
    return (mean + sqrt(var) * v1 * sqrt(-2.0 * log(r)/r));
}

int gera_bits(float p)
{
    //Gera um bit aleatorio com p(0)= p e p(1)= 1-p
    
    if (unif()<p)
        return 0;
    else
        return 1;
}

#endif /* _RAND_H_ */

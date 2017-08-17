//
//  Sort.h
//  Util
//
//  Created by Felipe Ferreira on 5/4/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef Sort_h
#define Sort_h

namespace Util
{
    void BubbleSort(int *vector, int size);
}

void Util::BubbleSort(int *vector, int size)
{
    int i, j;
    int swap;
    
    for (i = 0; i < size - 1; i++)
    {
        for (j = i + 1; j < size; j++)
        {
            if (vector[i] > vector[j])
            {
                swap = vector[i];
                vector[i] = vector[j];
                vector[j] = swap;
            }
        }
    }
}

#endif /* Sort_h */

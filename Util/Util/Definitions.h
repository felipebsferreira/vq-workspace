//
//  Definitions.h
//  Util
//
//  Created by Felipe Ferreira on 5/5/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef Definitions_h
#define Definitions_h

namespace ImageEnhancement
{
    enum MASK
    {
        LAPLACIAN,
        LAPLACIANDIAGONAL
    };
    
    const int LaPlacianMaskSize = 3;
    const int LaPlacianDiagonalSize = 3;
    
    const int LaPlacianMask[3][3] =
    {
         0, -1,  0,
        -1,  5, -1,
         0, -1,  0
    };
    
    const int LaPlacianDiagonalMask[3][3] =
    {
        -1, -1, -1,
        -1,  9, -1,
        -1, -1, -1
    };
}

#endif /* Definitions_h */

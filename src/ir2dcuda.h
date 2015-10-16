/*
 * Copyright (c) 2013-2015, Florent Hédin, Pierre-André Cazade, Markus Meuwly, 
 * and the University of Basel, Switzerland. All rights reserved.
 * 
 * The 3-clause BSD license is applied to this software.
 * 
 * See LICENSE.txt
 */

#ifndef IR2DCUDA_H_INCLUDED
#define IR2DCUDA_H_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus
/**
 *  2D IR spectrum calculations using a cuda gpu with CUFFT
 */
int ir2d_cuda(FILE* input, FILE* output);

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //IR2DCUDA_H_INCLUDED


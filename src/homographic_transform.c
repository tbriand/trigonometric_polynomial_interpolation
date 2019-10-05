/* SPDX-License-Identifier: GPL-2.0+
 * 
 * Thibaud Briand <briand.thibaud@gmail.com>
 * 
 * Copyright (c) 2018-2019, Thibaud Briand
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>

#include "homography_core.h"
#include "fft_core.h"
#define NFFT_PRECISION_DOUBLE
#include "external/nfft-3.5.0/include/nfft3mp.h"

#define N_MULTIPL 2
#define M_POLYDEG 6

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

// Compute the correspondences between positions in [0,nx) x [0,ny) (DFT convention)
// and positions in [-1/2,1/2)^2 (NDFT convention).
// This corresponds to Line 2 of Algorithm 2 (or Equation (51)). 
static void init_position(int nx, int ny, double *x, double *y, int numPixels, nfft_plan *my_plan)
{
    // sanity check
    assert((int) my_plan->M_total == numPixels);
    
    // rescaling parameters 
    double scaleX = 1.0/((double) nx);  
    double scaleY = 1.0/((double) ny);
    
    double ex, ey;
    for (int i = 0; i < numPixels; i++) {
        //rescale
        ex = x[i]*scaleX;
        ey = y[i]*scaleY;
        
        // periodization
        while ( ex < 0 )
            ex++;
        while ( ey < 0 )
            ey++;
        while ( ex >= 1 )
            ex--;
        while ( ey >= 1 )
            ey--;
        
        // set to [-0.5,0.5] and do not forget the minus in the NDFT convention
        ex = (ex <= 0.5) ? - ex : 1 - ex;
        ey = (ey <= 0.5) ? - ey : 1 - ey;
        my_plan->x[2*i]   =  ey;
        my_plan->x[2*i+1] =  ex;
    }
}

// Compute the next power of 2.
static unsigned int next_power_of_2(unsigned int v)
{
   v--;
   v |= v >> 1;
   v |= v >> 2;
   v |= v >> 4;
   v |= v >> 8;
   v |= v >> 16;
   v++;
   return v;
}

// Initialization of the NFFT plan
// The values Xband, Yband are the sizes of the spectrum for the input function, 
// m: is the parameter for selection the interpolation function
//
// AFTER INITIALIZING THE KNOTS ARE FIXED, ONLY CAN BE CHANGED THE COORDINATES
static void irregular_sampling_init(int Xband, int Yband, int num_knots, double n_multiplier, int m, nfft_plan *my_plan) {
    int my_N[2], my_n[2];

    // Nasty workarround for the NFFT problem with odd bandwidths
    // THE SOLUTION: is to extend by 1 the spectral dimension that 
    // is odd by allocation zeros. It will not affect the result
    // but the results must be extracted carefully
    // This corresponds to Line 5 of Algorithm 2 for the particular case (or Equation (54))
    if (Yband == 1) 
        my_N[0] = Yband;
    else
        my_N[0] = (int) ceil((double) Yband/2)*2;
    if (Xband==1)
        my_N[1] = Xband;
    else 
        my_N[1] = (int) ceil((double) Xband/2)*2;

    my_n[0] = next_power_of_2((int)(Yband*n_multiplier));
    my_n[1] = next_power_of_2((int)(Xband*n_multiplier));

    // window function m
    // -----------------------
    // Kaiser-Bessel         6
    // $ {\rm sinc}^{2m}$    9
    // $ B$-spline          11
    // Gaussian             12

    // PRE_PHI_HUT| PRE_PSI| PRE_FULL_PSI| MALLOC_X| MALLOC_F_HAT|
    // nfft_init_specific (PLAN, dimension, N (number of fourier coefficients in each dimension),
    // M (irregular knots to evaluate),
    // n (number of fourier coefficients computed for the interpolation, one for each dimension) ,
    // m (cut off parameter in time domain) 
    nfft_init_guru(my_plan, 2, my_N, num_knots,  my_n, m,
                   MALLOC_X| MALLOC_F_HAT|
                   MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE,
                   FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
}

// Compute the irregular samples of f given in Equation (50) from fhat using the NFFT algorithm.
// When a dimension is odd an extra frequency with zeros is added.
// The coefficients in fhat cannot be used directly from the FFTW as the fftshift function must be applied.
// This corresponds to Line 5 to 7 of Algorihtm 2.
static void irregular_sampling_fourier(int nx, int ny, const fftw_complex *fhat, double *out, nfft_plan *my_plan)
{
    // nx and ny are the bandwith of the input transform fhat
    int numknots = my_plan->M_total;
    int Yband = my_plan->N[0];
    int Xband = my_plan->N[1];

    // the values of difx are 0 or 1, depending if we added or not
    // an extra frequency (with zeros)
    int difx = Xband - nx;  
    int dify = Yband - ny;
    
    // load the fourier coefficients (Line 5 and Line 6 of Algorithm 2)
    for (int i = 0; i < Xband; i++)
        for (int j = 0; j < Yband; j++)  {
            int pos = i + Xband *j;
            int posIn = i - difx + nx * (j-dify);
            if ( ((i-difx)>=0 ) && ((j-dify)>=0) )
                    my_plan->f_hat[pos] = fhat[posIn];
            else
                    my_plan->f_hat[pos] = 0.0;
    }
    
    // execute NFFT
    nfft_trafo(my_plan);
    
    // Extract the results and normalize the values
    for (int i = 0; i < numknots; i++)
            out[i] = creal(my_plan->f[i]) / (nx*ny);
}

// Transformation of an image using trigonometric polynomial interpolation (Line 2 to 7 of Algorithm 2).
static void interpolate_at_locations_nfft(double *out, double *in, int nx, int ny, int nz,
                                   double *x, double *y, int numPixels, int interp)
{
    // plan initialization
    NFFT(plan) my_plan;
    irregular_sampling_init(nx, ny, numPixels, N_MULTIPL, M_POLYDEG, &my_plan);
    init_position(nx, ny, x, y, numPixels, &my_plan);

    // allocate memory for fourier transform
    fftw_complex *fhat = fftw_malloc(nx*ny*nz*sizeof*fhat);
    fftw_complex *fshift = fftw_malloc(nx*ny*nz*sizeof*fshift);

    // compute DFT of the input
    do_fft_real(fhat, in, nx, ny, nz);

    // fftshift (to have the right ordering of polynomial coefficients)
    fftshift(fshift, fhat, nx, ny, nz);
    free(fhat);

    // evaluation of the interpolated values for each channel
    for(int l = 0; l < nz; l++) {
        irregular_sampling_fourier(nx, ny, fshift + l*nx*ny, out + l*numPixels, &my_plan);
        
        // real convention adjustment using Equation (27)
        if( interp && !(nx%2) && !(ny%2) ) {
            double hf = creal(fshift[l*nx*ny])/(nx*ny);
            for(int i = 0; i < numPixels; i++)
                out[i + l*numPixels] += hf*sin(M_PI*x[i])*sin(M_PI*y[i]);
        }
    }

    //cleanup
    nfft_finalize(&my_plan);

    //free memory
    free(fshift);
}

// Homographic transformation of an image using trigonometric polynomial interpolation.
// This corresponds to Algorithm 2 when $\varphi$ is an homography.
void interpolate_image_homography(double *out, double *in, int nx, int ny, int nz, const double H[9], int interp) {
    
    // invert the homography
    double iH[9];
    invert_homography(iH, H);
    
    // number of output pixels
    int numPixels = nx*ny;

    // compute pixel locations (x,y) = H^{-1}(i,j)
    double *x = malloc(numPixels*sizeof*x);
    double *y = malloc(numPixels*sizeof*y);
    double p[2], q[2];
    for (int j = 0; j < ny; j++) {
        p[1] = j;
        for (int i = 0; i < nx; i++) {
            p[0] = i;
            apply_homography(q, p, iH);
            x[j*nx + i] = q[0];
            y[j*nx + i] = q[1];
        }
    }

    // interpolation at locations (x[i], y[i]) using NFFT
    interpolate_at_locations_nfft(out, in, nx, ny, nz, x, y, numPixels, interp);
    
    // free memory
    free(x);
    free(y);
}

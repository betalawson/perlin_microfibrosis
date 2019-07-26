
/*
 * Octave2D.cpp
 *
 * Creates an 'octave noise' pattern (multiple octaves of basic Perlin
 * noise, evaluated at a set of points provided by the user. A grid with
 * unit spacing is used for the Perlin grid - input points should be
 * transformed appropriately to generate features of desired size and 
 * orientation. The user may also specify the number of octaves, and the
 * factor successively applied to each octave. This seeded version of the
 * code also takes as input a permutation table of the numbers 0 to 255
 * and an m x 2 set of offsets in 2D space to apply to the Perlin grid for
 * each octave. m must therefore be at least as large as N_freq.
 *
 * Calling syntax in MATLAB is
 *
 *		noisefield = Octave2D(points, N_freq, fade_factor, Ps, offset)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "math.h"

                    
//Define Perlin vectors (256 points distributed evenly around circle)
static double vecs_x[256] = 
    {1.0000, 0.9997, 0.9988, 0.9973, 0.9952, 0.9925, 0.9892, 0.9853, 0.9808,
     0.9757, 0.9700, 0.9638, 0.9569, 0.9495, 0.9415, 0.9330, 0.9239, 0.9142,
     0.9040, 0.8932, 0.8819, 0.8701, 0.8577, 0.8449, 0.8315, 0.8176, 0.8032,
     0.7883, 0.7730, 0.7572, 0.7410, 0.7242, 0.7071, 0.6895, 0.6716, 0.6532,
     0.6344, 0.6152, 0.5957, 0.5758, 0.5556, 0.5350, 0.5141, 0.4929, 0.4714,
     0.4496, 0.4276, 0.4052, 0.3827, 0.3599, 0.3369, 0.3137, 0.2903, 0.2667,
     0.2430, 0.2191, 0.1951, 0.1710, 0.1467, 0.1224, 0.0980, 0.0736, 0.0491,
     0.0245, 0.0000, -0.0245, -0.0491, -0.0736, -0.0980, -0.1224, -0.1467, -0.1710,
     -0.1951, -0.2191, -0.2430, -0.2667, -0.2903, -0.3137, -0.3369, -0.3599, -0.3827,
     -0.4052, -0.4276, -0.4496, -0.4714, -0.4929, -0.5141, -0.5350, -0.5556, -0.5758,
     -0.5957, -0.6152, -0.6344, -0.6532, -0.6716, -0.6895, -0.7071, -0.7242, -0.7410,
     -0.7572, -0.7730, -0.7883, -0.8032, -0.8176, -0.8315, -0.8449, -0.8577, -0.8701,
     -0.8819, -0.8932, -0.9040, -0.9142, -0.9239, -0.9330, -0.9415, -0.9495, -0.9569,
     -0.9638, -0.9700, -0.9757, -0.9808, -0.9853, -0.9892, -0.9925, -0.9952, -0.9973,
     -0.9988, -0.9997, -1.0000, -0.9997, -0.9988, -0.9973, -0.9952, -0.9925, -0.9892,
     -0.9853, -0.9808, -0.9757, -0.9700, -0.9638, -0.9569, -0.9495, -0.9415, -0.9330,
     -0.9239, -0.9142, -0.9040, -0.8932, -0.8819, -0.8701, -0.8577, -0.8449, -0.8315,
     -0.8176, -0.8032, -0.7883, -0.7730, -0.7572, -0.7410, -0.7242, -0.7071, -0.6895,
     -0.6716, -0.6532, -0.6344, -0.6152, -0.5957, -0.5758, -0.5556, -0.5350, -0.5141,
     -0.4929, -0.4714, -0.4496, -0.4276, -0.4052, -0.3827, -0.3599, -0.3369, -0.3137,
     -0.2903, -0.2667, -0.2430, -0.2191, -0.1951, -0.1710, -0.1467, -0.1224, -0.0980,
     -0.0736, -0.0491, -0.0245, -0.0000, 0.0245, 0.0491, 0.0736, 0.0980, 0.1224,
     0.1467, 0.1710, 0.1951, 0.2191, 0.2430, 0.2667, 0.2903, 0.3137, 0.3369,
     0.3599, 0.3827, 0.4052, 0.4276, 0.4496, 0.4714, 0.4929, 0.5141, 0.5350,
     0.5556, 0.5758, 0.5957, 0.6152, 0.6344, 0.6532, 0.6716, 0.6895, 0.7071,
     0.7242, 0.7410, 0.7572, 0.7730, 0.7883, 0.8032, 0.8176, 0.8315, 0.8449,
     0.8577, 0.8701, 0.8819, 0.8932, 0.9040, 0.9142, 0.9239, 0.9330, 0.9415,
     0.9495, 0.9569, 0.9638, 0.9700, 0.9757, 0.9808, 0.9853, 0.9892, 0.9925,
     0.9952, 0.9973, 0.9988, 0.9997};
static double vecs_y[256] = 
    {0.0000, 0.0245, 0.0491, 0.0736, 0.0980, 0.1224, 0.1467, 0.1710, 0.1951,
     0.2191, 0.2430, 0.2667, 0.2903, 0.3137, 0.3369, 0.3599, 0.3827, 0.4052,
     0.4276, 0.4496, 0.4714, 0.4929, 0.5141, 0.5350, 0.5556, 0.5758, 0.5957,
     0.6152, 0.6344, 0.6532, 0.6716, 0.6895, 0.7071, 0.7242, 0.7410, 0.7572,
     0.7730, 0.7883, 0.8032, 0.8176, 0.8315, 0.8449, 0.8577, 0.8701, 0.8819,
     0.8932, 0.9040, 0.9142, 0.9239, 0.9330, 0.9415, 0.9495, 0.9569, 0.9638,
     0.9700, 0.9757, 0.9808, 0.9853, 0.9892, 0.9925, 0.9952, 0.9973, 0.9988,
     0.9997, 1.0000, 0.9997, 0.9988, 0.9973, 0.9952, 0.9925, 0.9892, 0.9853,
     0.9808, 0.9757, 0.9700, 0.9638, 0.9569, 0.9495, 0.9415, 0.9330, 0.9239,
     0.9142, 0.9040, 0.8932, 0.8819, 0.8701, 0.8577, 0.8449, 0.8315, 0.8176,
     0.8032, 0.7883, 0.7730, 0.7572, 0.7410, 0.7242, 0.7071, 0.6895, 0.6716,
     0.6532, 0.6344, 0.6152, 0.5957, 0.5758, 0.5556, 0.5350, 0.5141, 0.4929,
     0.4714, 0.4496, 0.4276, 0.4052, 0.3827, 0.3599, 0.3369, 0.3137, 0.2903,
     0.2667, 0.2430, 0.2191, 0.1951, 0.1710, 0.1467, 0.1224, 0.0980, 0.0736,
     0.0491, 0.0245, 0.0000, -0.0245, -0.0491, -0.0736, -0.0980, -0.1224, -0.1467,
     -0.1710, -0.1951, -0.2191, -0.2430, -0.2667, -0.2903, -0.3137, -0.3369, -0.3599,
     -0.3827, -0.4052, -0.4276, -0.4496, -0.4714, -0.4929, -0.5141, -0.5350, -0.5556,
     -0.5758, -0.5957, -0.6152, -0.6344, -0.6532, -0.6716, -0.6895, -0.7071, -0.7242,
     -0.7410, -0.7572, -0.7730, -0.7883, -0.8032, -0.8176, -0.8315, -0.8449, -0.8577,
     -0.8701, -0.8819, -0.8932, -0.9040, -0.9142, -0.9239, -0.9330, -0.9415, -0.9495,
     -0.9569, -0.9638, -0.9700, -0.9757, -0.9808, -0.9853, -0.9892, -0.9925, -0.9952,
     -0.9973, -0.9988, -0.9997, -1.0000, -0.9997, -0.9988, -0.9973, -0.9952, -0.9925,
     -0.9892, -0.9853, -0.9808, -0.9757, -0.9700, -0.9638, -0.9569, -0.9495, -0.9415,
     -0.9330, -0.9239, -0.9142, -0.9040, -0.8932, -0.8819, -0.8701, -0.8577, -0.8449,
     -0.8315, -0.8176, -0.8032, -0.7883, -0.7730, -0.7572, -0.7410, -0.7242, -0.7071,
     -0.6895, -0.6716, -0.6532, -0.6344, -0.6152, -0.5957, -0.5758, -0.5556, -0.5350,
     -0.5141, -0.4929, -0.4714, -0.4496, -0.4276, -0.4052, -0.3827, -0.3599, -0.3369,
     -0.3137, -0.2903, -0.2667, -0.2430, -0.2191, -0.1951, -0.1710, -0.1467, -0.1224,
     -0.0980, -0.0736, -0.0491, -0.0245, 
    };
    
typedef struct Vector2_s {
   double x, y;
} Vector2;
                     
/* Linear interpolation function. alpha is the position (alpha = 0 takes
   only starting component, alpha = 1 takes only ending component) */
inline double linearInterpolate(double fstart, double fend, double alpha)
{
    return fstart + alpha * (fend - fstart);
}

/* Smoothing function. Takes an input co-ordinate and smoothes it, in the
 * sense that proximity to 0 or 1 pushes the value closer to that extreme
 * limit. The function used is:
     x_smooth = 6 x^5 - 15 x^4 + 10 x^3                  */
inline double smoothInterpolate(double fstart, double fend, double alpha)
{
    double alpha_smoothed = alpha * alpha * alpha * (10 - alpha * (15 - 6 * alpha) );
    return linearInterpolate(fstart, fend, alpha_smoothed);
}

/* Dot product function. Takes the dot product of two 2D vectors */
inline double dotProduct2D(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.x + v1.y * v2.y);
}

/* Random vector selectio nfunction. Takes as input two integers, and uses
   hashing to convert these into a 'random' selection from 0-255 */
inline int selectVector(unsigned int x, unsigned int y, int *P)
{

    // Hash with bitwise XOR
    return P[(x % 256)^P[(y % 256)]];
}

/* Noise calculation function. Takes an input point in 2D space and
 calculates the dot product with the four vectors surrounding this point */
double noise2D(double x, double y, int *P)
{

    // Convert the input x and y into the corresponding integers
    int int_x = (int)floor(x);
    int int_y = (int)floor(y);
    // Use these to find the "box co-ordinates (co-ords relative to bottom corner of box)
    double box_x = x - int_x;
    double box_y = y - int_y;
    
    // Define vectors to the corners
    Vector2 v_dl = {-box_x, -box_y};
    Vector2 v_dr = {1-box_x, -box_y};
    Vector2 v_ul = {-box_x, 1-box_y};
    Vector2 v_ur = {1-box_x, 1-box_y};
    
    // Use hash function to assign gradient vectors to the corners
    unsigned int k_dl = selectVector(int_x, int_y, P);
    unsigned int k_dr = selectVector(int_x+1, int_y, P);
    unsigned int k_ul = selectVector(int_x, int_y+1, P);
    unsigned int k_ur = selectVector(int_x+1, int_y+1, P);
    
    Vector2 g_dl = {vecs_x[k_dl], vecs_y[k_dl]}; 
    Vector2 g_dr = {vecs_x[k_dr], vecs_y[k_dr]};
    Vector2 g_ul = {vecs_x[k_ul], vecs_y[k_ul]};
    Vector2 g_ur = {vecs_x[k_ur], vecs_y[k_ur]};
    
    // Calculate the values of the dot products
    double dl_val = dotProduct2D(v_dl, g_dl);
    double dr_val = dotProduct2D(v_dr, g_dr);
    double ul_val = dotProduct2D(v_ul, g_ul);
    double ur_val = dotProduct2D(v_ur, g_ur);

    // Perform linear interpolation of these dot products and return the result
    double bottom_noise = smoothInterpolate(dl_val, dr_val, box_x);
    double top_noise = smoothInterpolate(ul_val, ur_val, box_x);
    return smoothInterpolate(bottom_noise, top_noise, box_y);
    
}


/* Main runner function. Takes as input the points where Perlin noise is to
 * be evaluated, and calculates the value of octave noise at each point,
 * which is the combination of multiple 'octaves' of Perlin noise. The
 * result is scaled back to lie within [0,1] according to the theoretical
 * maximum and minimum values of 2D octave noise. */
void OctaveNoise2D(double *points, double *noisefield, mwSize n, unsigned int N_freq, double fade_factor, int *Ps, double *offsets)
{
    
    // Initialise loop variables
    mwSize i;
    unsigned int j;
    
    // Initialise multipliers
    int freq_mult;
    double scale_mult;
    
    // Initialise individual permutation tables
    int *P;
    
    // Calculate the scaling factor for normalisation
    // This is calculated using a geometric progression, with starting value sqrt(2)/2
    double normalisation_factor = 0.70710678 * (1 - (double)pow((double)fade_factor, (int)N_freq) ) / (1 - fade_factor);
    double normalisation_denominator = 0.5 / normalisation_factor;
    
    // Perform a loop over all input points
    for (i=0; i<n; i++) {
        
        // Initialise value of noisefield here as zero
        noisefield[i] = 0;
        
        // Loop over frequencies
        for (j=0; j<N_freq; j++) {
         
            // Set up multiplier for this frequency
            freq_mult = (int) pow((double)2 ,int(j));
            scale_mult = (double)pow((double)fade_factor, int(j));
            
            // Grab out this octave's permutation table (pointer pointing to shifted position along the full array Ps)
            P = &Ps[256*j];
                        
            // Add value of noisefield at this location
            noisefield[i] = noisefield[i] + scale_mult*noise2D(points[2*i]*freq_mult - offsets[2*j], points[2*i+1]*freq_mult - offsets[2*j+1], P);
            
        }
        
        // Convert this value of the noisefield to its [0,1] normalised equivalent
        noisefield[i] = ( noisefield[i] + normalisation_factor ) * normalisation_denominator;
        
    }
}

/* "Main" function used for interface between C and MATLAB */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    // Read in or otherwise initialise variables
    size_t n_points = mxGetN(prhs[0]);         // points supplied as rows
    double *points = mxGetPr(prhs[0]);
    unsigned int N_freq = mxGetScalar(prhs[1]);
    double fade_factor = mxGetScalar(prhs[2]);
    int *Ps = (int *)mxGetData(prhs[3]);
    double *offsets = mxGetPr(prhs[4]);
    double *noisefield;
        
    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)n_points,mxREAL);

    // Create a pointer to this output matrix
    noisefield = mxGetPr(plhs[0]);
    
    // Call the PerlinNoise2D function, which will update the values of noisefield
    OctaveNoise2D(points, noisefield, (mwSize)n_points, N_freq, fade_factor, Ps, offsets);
    
}
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <random>
#include <omp.h>

using namespace std;

//the version for SMLM analysis
void iter(double* newImg, const double* hiImg, const size_t *dims, const uint32_t *locs, const size_t numLocs, const double * psfs,const size_t psfSize, int begin=0, int finish=-1)
{
    uniform_real_distribution<> rand(0, 1);
    
    size_t psfHSize = (psfSize-1)/2;
    if (finish == -1) finish = numLocs;

    fill(newImg, newImg + dims[0] * dims[1], 0.5);

    #pragma omp parallel for
    for (int i = begin; i < finish; i ++) {
        double bins[psfSize*psfSize];
        
        uint32_t row = locs[i*3+1] - psfHSize;
        uint32_t col = locs[i*3] - psfHSize;
        const double * psf = psfs + locs[i*3+2] * psfSize * psfSize;
        
        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1]) continue;

        unsigned int idx = 0;
        double s = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx] = hiImg[c * dims[0] + r] * psf[idx];
                s += bins[idx];
                idx ++;
            }
        }

        idx = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                int tmp = c * dims[0] + r;
                #pragma omp atomic
                    newImg[tmp] += bins[idx] / s;
                idx ++;
            }
        }
    }
}

// Input -
// #1 - UINT32 Data. Either a low resolution image, or a 3XN array of SMLM data. 
//      If it's latter, row 1 is x, row 2 is y, row 3 is a index to psf array (see input #3)
// #2 - Image from last itration
// #3 - PSF. For image analysis, it's 2D double array. For SMLM, it's a SxSxN 3D array, representing
//      N different kinds of PSFs with varying localization accuracy. This is because in SMLM different
//      molecules are localized with different accuracy.
// Output -
// #1 - Image result from new iteration.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mwSize d = mxGetNumberOfDimensions(prhs[2]);
    
    const mwSize *datadims = mxGetDimensions(prhs[0]);
    const mwSize *imgdims = mxGetDimensions(prhs[1]);
    const mwSize *psfdims = mxGetDimensions(prhs[2]);
    
    const mxUint32 * data = mxGetUint32s(prhs[0]); 
    const mxDouble * img = mxGetDoubles(prhs[1]);
    const mxDouble * psf = mxGetDoubles(prhs[2]);

    // if psf is 2 2D array, assuming input data is a image
    if ( d == 2 ) {
        //int scale = imgdims[0] / datadims[0];

        //plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
        //mxDouble * newImg = mxGetDoubles(plhs[0]);
        //uint32_t * tmpImg1 = new uint32_t[imgdims[0] * imgdims[1]];
        //uint32_t * tmpImg2 = new uint32_t[imgdims[0] * imgdims[1]];

        // sas(tmpImg1, data, datadims, scale);
        // sample_r(tmpImg2, img, imgdims, tmpImg1, psf, psfdims[0]);
        // sample_t(newImg, tmpImg2, imgdims);

        //delete [] tmpImg1;
        //delete [] tmpImg2;
    } else { // else psf is 3D arrray, input data is single molecule localizations
        plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
        mxDouble * newImg = mxGetDoubles(plhs[0]);

        iter(newImg, img, imgdims, data, datadims[1], psf, psfdims[0]);
    }
}

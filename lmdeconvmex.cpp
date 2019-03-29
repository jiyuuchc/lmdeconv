#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdint>
#include <algorithm>
#include <omp.h>

using namespace std;

//For images
void iter(double* newImg, const double* hiImg, const uint32_t * loImg, const size_t *dims , const double * psf, const size_t psfSize, int beginCol=0, int endCol = -1)
{
    size_t numel = dims[0] * dims[1];
    size_t psfHSize = (psfSize - 1) / 2;

    if (endCol == -1) endCol = dims[1];

    #pragma omp parallel for
    for (int i = beginCol * dims[0]; i < endCol * dims[0]; i ++) {

        double bins[psfSize*psfSize];
        bins[0] = 0;
    
        if (loImg[i] == 0) continue;

        size_t row = i % dims[0] - psfHSize;
        size_t col = i / dims[0] - psfHSize;
        
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

        if (s > 0) {
            idx = 0;
            for (int c = col ; c < col + psfSize; c++) {
                for (int r = row ; r < row + psfSize; r++) {
                    int tmp = c * dims[0] + r;
                    #pragma omp atomic
                        newImg[tmp] += double(loImg[i]) * bins[idx] / s;
                    idx ++;
                }
            }
        }
    }    
}

//the version for SMLM analysis
void iter(double* newImg, const double* hiImg, const size_t *dims, const uint32_t *locs, const size_t numLocs, const double * psfs,const size_t psfSize, int begin=0, int finish=-1)
{
    size_t psfHSize = (psfSize-1)/2;
    if (finish == -1) finish = numLocs;

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
    if (nrhs != 3 && nrhs != 4) {
        mexErrMsgTxt("lmdeconvmex: incorrect number of inputs");
    }

    plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
    mxDouble * newImg = mxGetDoubles(plhs[0]);

    fill(newImg, newImg + mxGetNumberOfElements(plhs[0]), 0);

    // if psf is 2 2D array, assuming input data is a image
    if ( d == 2 ) {
        iter(newImg, img, data, imgdims, psf, psfdims[0]);
    } else { // else psf is 3D arrray, input data is single molecule localizations
        iter(newImg, img, imgdims, data, datadims[1], psf, psfdims[0]);
    }

    if (nrhs == 4) {
        if (mxIsScalar(prhs[3])) {
            for (int i = 0 ; i < mxGetNumberOfElements(plhs[0]); i++) {
                newImg[i] += mxGetScalar(prhs[3]) - 1.0;
                if (newImg[i] < 0) {
                    newImg[i] = 0;
                }
            }
        } else {
            mxDouble * prior = mxGetDoubles(prhs[3]);
            for (int i = 0; i < mxGetNumberOfElements(plhs[0]); i ++) {
                newImg[i] += prior[i] - 1.0;
                if (newImg[i] < 0) {
                    newImg[i] = 0;
                }
            }
        }
    }
}

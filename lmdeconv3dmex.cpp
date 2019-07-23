#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdint>
#include <algorithm>
#include <omp.h>

using namespace std;

void iter3d(double* newImg, const double* hiImg, const size_t *dims, const uint32_t *locs, const size_t numLocs, const double * psfs,const size_t psfSize, const double *psfZs, const size_t psfZSize, int begin=0, int finish=-1)
{
    //mexPrintf("%d , %d, %d\n", dims[0], dims[1], dims[2]);
    size_t psfHSize = (psfSize-1)/2;
    size_t psfZHSize = (psfZSize-1)/2;
    if (finish == -1) finish = numLocs;

    /*#pragma omp parallel for*/
    for (int i = begin; i < finish; i ++) {
        double bins[psfSize*psfSize*psfZSize];

        uint32_t row = locs[i*5+1] - psfHSize;
        uint32_t col = locs[i*5] - psfHSize;
        uint32_t slice = locs[i*5+3] - psfZHSize;
        uint32_t psfInd = locs[i*5 + 2];
        uint32_t psfZInd = locs[i*5 + 4];

        const double * psf = psfs + psfInd * psfSize * psfSize;
        const double * psfz = psfZs + psfZInd * psfZSize;

        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1] || slice < 0 || slice + psfZSize >= dims[2]) {
            continue;
        }

        unsigned int idx0 = 0, idx1;
        double s = 0;
        for (int z = slice; z < slice + psfZSize; z ++) {
            idx1 = 0;
            for (int c = col ; c < col + psfSize; c++) {
                for (int r = row ; r < row + psfSize; r++) {
                    //mexPrintf("%d - %d - %d\n", z * dims[0] * dims[1] + c * dims[0] + r, idx0, idx1);
                    bins[idx0] = hiImg[z * dims[0] * dims[1] + c * dims[0] + r] * psf[idx1] * psfz[z-slice];
                    s += bins[idx0];
                    idx0 ++;
                    idx1 ++;
                }
            }
        }
        idx0 = 0;
        if ( s > 0 ) {
            for (int z = slice; z < slice + psfZSize; z ++) {
                for (int c = col ; c < col + psfSize; c++) {
                    for (int r = row ; r < row + psfSize; r++) {
                        int tmp = z * dims[0] * dims[1] + c * dims[0] + r;
                        double cnt = bins[idx0] / s;
                        /*#pragma omp atomic*/
                            newImg[tmp] += cnt;
                        idx0 ++;
                    }
                }
            }
        }
    }
}

// Input -
// #1 - UINT32 Data. 5XN array of SMLM data. 
//      row 1 is x, row 2 is y, row 3 is z, row 4 is a index to xypsf, row 5 is a index to zpsf
// #2 - Image from last itration
// #3 - xy-PSF. 3D array, representing N different kinds of PSFs with varying localization accuracy.
// #4 - z-PSF. 2D array.
// #5 - optional prior 
// Output -
// #1 - Image result from new iteration.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mwSize imgND = mxGetNumberOfDimensions(prhs[1]); // image dimension

    const mwSize *datadims = mxGetDimensions(prhs[0]);
    const mwSize *imgdims = mxGetDimensions(prhs[1]);
    const mwSize *psfxydims = mxGetDimensions(prhs[2]);
    const mwSize *psfzdims = mxGetDimensions(prhs[3]);

    const mxUint32 * data = mxGetUint32s(prhs[0]);
    const mxDouble * img = mxGetDoubles(prhs[1]);
    const mxDouble * psfxy = mxGetDoubles(prhs[2]);
    const mxDouble * psfz = mxGetDoubles(prhs[3]);
    
    plhs[0] = mxCreateNumericArray(imgND, imgdims, mxDOUBLE_CLASS, mxREAL);
    mxDouble * newImg = mxGetDoubles(plhs[0]);
    fill(newImg, newImg + mxGetNumberOfElements(plhs[0]), 0);

    //mexPrintf("%d , %d, %d\n", imgdims[0], imgdims[1], imgdims[2]);
    iter3d(newImg, img, imgdims, data, datadims[1], psfxy, psfxydims[0], psfz, psfzdims[0]);

    if (nrhs == 5) {
        if (mxIsScalar(prhs[4])) {
            for (int i = 0 ; i < mxGetNumberOfElements(plhs[0]); i++) {
                newImg[i] += mxGetScalar(prhs[4]) - 1.0;
                if (newImg[i] < 0) {
                    newImg[i] = 0;
                }
            }
        } else {
            mxDouble * prior = mxGetDoubles(prhs[4]);
            for (int i = 0; i < mxGetNumberOfElements(plhs[0]); i ++) {
                newImg[i] += prior[i] - 1.0;
                if (newImg[i] < 0) {
                    newImg[i] = 0;
                }
            }
        }
    }
}

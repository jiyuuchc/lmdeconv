#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdint>
#include <algorithm>
#include <omp.h>

using namespace std;
static void init() __attribute__ ((constructor));

void init()
{
    mexPrintf("LMDeconvmex:init with %d threads\n",omp_get_max_threads());
}

void iter(double* newImg, const mxArray * data, const mxArray * img, const mxArray * psfs)
{
    const mwSize * imgDims = mxGetDimensions(img);
    double * imgBuf = mxGetDoubles(img);
    uint32_t * dataBuf = mxGetUint32s(data);
    size_t nLocs = mxGetDimensions(data)[1];

    //size_t psfHSize = (psfSize-1)/2;

    #pragma omp parallel for
    for (int i = 0; i < nLocs; i ++) {
        uint32_t psfIdx = dataBuf[i*3+2];
        mxArray * psf = mxGetCell(psfs, psfIdx);
        double * psfBuf = mxGetDoubles(psf);
        size_t psfSize = mxGetDimensions(psf)[0];
        size_t psfHSize = (psfSize-1)/2;

        double * bins = new double[psfSize*psfSize];
        
        int32_t row = dataBuf[i*3+1] - psfHSize;
        int32_t col = dataBuf[i*3] - psfHSize;

        if (row < 0 || col < 0 || row + psfSize >= imgDims[0] || col + psfSize >= imgDims[1]) continue;

        unsigned int idx = 0;
        double s = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx] = imgBuf[c * imgDims[0] + r] * psfBuf[idx];
                s += bins[idx];
                idx ++;
            }
        }

        idx = 0;
        if ( s > 0 ) {
            for (int c = col ; c < col + psfSize; c++) {
                for (int r = row ; r < row + psfSize; r++) {
                    int tmp = c * imgDims[0] + r;
                    double cnt = bins[idx] / s;
                    #pragma omp atomic
                        newImg[tmp] += cnt;
                    idx ++;
                }
            }
        }
        
        delete [] bins;
    }
}

// Input -
// #1 - UINT32 Data. 3XN array of SMLM data. 
//      row 1 is x, row 2 is y, row 3 is a index to psf array (see input #3)
// #2 - Image from last itration
// #3 - PSFs. Cell array. Eepresenting N different kinds of PSFs with varying localization accuracy. 
// #4 - Prior. Scalar or array. Optional. 
// Output -
// #1 - Image result from new iteration.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mwSize *imgDims = mxGetDimensions(prhs[1]);
    if (nrhs != 3 && nrhs != 4) {
        mexErrMsgTxt("lmdeconvmex: incorrect number of inputs");
    }

    plhs[0] = mxCreateDoubleMatrix(imgDims[0], imgDims[1], mxREAL);
    mxDouble * newImg = mxGetDoubles(plhs[0]);
    fill(newImg, newImg + mxGetNumberOfElements(plhs[0]), 0);

    iter(newImg, prhs[0], prhs[1], prhs[2]);

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

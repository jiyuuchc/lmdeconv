#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <random>
#include <omp.h>

using namespace std;

static void init() __attribute__ ((constructor));
static void fini() __attribute__ ((destructor));

mt19937_64 * mts;
const uint32_t INTENSITY_FOR_POISS_APPROX = 5000;
        
void init() 
{
    random_device rd;
    mts = new mt19937_64[omp_get_max_threads()];
    for (int i = 0; i < omp_get_max_threads(); i++) 
        mts[i].seed(rd());
    mexPrintf("LMSampler:init with %d threads\n",omp_get_max_threads());
}

void fini()
{
    delete [] mts;
    mexPrintf("LMSampler:exit\n");
}

void sample_r(uint32_t * newImg, const mxArray * data, const mxArray * sample, const mxArray * psfs) {
    uniform_real_distribution<> rand(0, 1);
    const mwSize * imgDims = mxGetDimensions(sample);
    const mwSize * dataDims = mxGetDimensions(data);
    size_t nLocs = dataDims[1];

    uint32_t * dataBuf = mxGetUint32s(data);
    double * sampleBuf = mxGetDoubles(sample);
        
    // fill(newImg, newImg + dims[0] * dims[1], 0);

    #pragma omp parallel for
    for (int i = 0; i < nLocs; i ++) {
        uint32_t psfIdx = dataBuf[i*3+2];
        mxArray * psf = mxGetCell(psfs, psfIdx);
        double * psfBuf = mxGetDoubles(psf);
        size_t psfSize = mxGetDimensions(psf)[0];
        size_t psfHSize = (psfSize-1)/2;

        double * bins = new double[psfSize*psfSize + 1];
        
        int32_t row = dataBuf[i*3+1] - psfHSize;
        int32_t col = dataBuf[i*3] - psfHSize;
        
        if (row < 0 || col < 0 || row + psfSize >= imgDims[0] || col + psfSize >= imgDims[1]) continue;

        int idx = 0;
        bins[0] = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx + 1] = bins[idx] + sampleBuf[c * imgDims[0] + r] * psfBuf[idx];
                idx ++;
            }
        }
        
        double r = rand(mts[omp_get_thread_num()]) * bins[psfSize * psfSize];
        
        idx = psfSize * psfSize / 2;
        int llimit = 1;
        int rlimit = psfSize * psfSize;
        while (1) {
            if (r >= bins[idx]) {
                llimit = idx + 1;                    
                idx += (rlimit - idx + 1) / 2;
            } else if ( r < bins[idx-1]) {
                rlimit = idx - 1;                    
                idx -= (idx - llimit + 1) / 2;
            } else {
                idx --;
                break;
            }
        }
        delete [] bins;

        int rr = idx % psfSize;
        int cc = idx / psfSize;
        int tmp = row + rr + (col+cc) * imgDims[0];

        #pragma omp atomic
        newImg[tmp] ++;
    }
}

void sample_t(mxArray *img, uint32_t * imgIn, double * prior)
{
    size_t numPixels = mxGetNumberOfElements(img);
    double * imgBuf = mxGetDoubles(img);

    double s = 0;
    
    #pragma omp parallel for reduction(+:s)
    for (int i = 0; i < numPixels; i++) {
        gamma_distribution<double> randgam(imgIn[i] + prior[i]);
        imgBuf[i] = randgam(mts[omp_get_thread_num()]);
        s += imgBuf[i];
    }

    #pragma omp parallel for
    for (int i = 0; i < numPixels; i++) {
        imgBuf[i] /= s;
    }    
}

// Input -
// #1 - UINT32 Data. 3XN array of SMLM data. 
//      row 1 is x, row 2 is y, row 3 is a index to psf array (see input #3)
// #2 - The image sample from last itration
// #3 - Cached PSFs. cell array, representing N different kinds of PSFs with 
//      varying localization accuracy
// #4 - Optional Prior. Either a scalar or a arrray/matrix. 
// Output -
// #1 - The sample. A double 2D array. It is normalized (summed to be 1.0).
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs != 3 && nrhs !=4) || nlhs != 1) {
        mexErrMsgTxt("lmsamplermex: Wrong numer of input/output arguments");
    }

    mxDouble * prior = new mxDouble[mxGetNumberOfElements(prhs[1])];
    if (nrhs == 3) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), 1.0);
    } else if (mxIsScalar(prhs[3])) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), mxGetScalar(prhs[3]));
    } else {
        mxDouble * prior_in = mxGetDoubles(prhs[3]);
        copy(prior_in, prior_in + mxGetNumberOfElements(prhs[1]), prior);
    }

    const mwSize * imgDims = mxGetDimensions(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(imgDims[0], imgDims[1], mxREAL);
    uint32_t * tmpImg = new uint32_t[imgDims[0] * imgDims[1]];
    fill(tmpImg, tmpImg + imgDims[0] * imgDims[1], 0);

    sample_r(tmpImg, prhs[0], prhs[1], prhs[2]);
    sample_t(plhs[0], tmpImg, prior);

    delete [] tmpImg;
    delete [] prior;
}

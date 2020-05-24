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
    // mexPrintf("LMSampler:init with %d threads\n",omp_get_max_threads());
}

void fini()
{
    delete [] mts;
    mexPrintf("SDS:exit\n");
}

void sample_r(uint32_t * newImg, const double* hiImg, const size_t *dims, const uint32_t *locs, const size_t numLocs, const double * psfs,const size_t psfSize, int begin=0, int finish=-1)
{
    uniform_real_distribution<> rand(0, 1);
    
    size_t psfHSize = (psfSize-1)/2;
    if (finish == -1) finish = numLocs;

    //mexPrintf("Sample_R \n");
    //mexPrintf("Img dimension %d - %d\n", dims[0], dims[1]);
    //mexPrintf("Psf dimension %d\n", psfSize);
    //mexPrintf("Num of locs %d\n", numLocs);
    
    fill(newImg, newImg + dims[0] * dims[1], 0);

    #pragma omp parallel for
    for (int i = begin; i < finish; i ++) {
        double bins[psfSize*psfSize + 1];
        
        int32_t row = locs[i*3+1] - psfHSize;
        int32_t col = locs[i*3] - psfHSize;
        const double * psf = psfs + locs[i*3+2] * psfSize * psfSize;
        
        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1]) continue;

        int idx = 0;
        bins[0] = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx + 1] = bins[idx] + hiImg[c * dims[0] + r] * psf[idx];
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
        int rr = idx % psfSize;
        int cc = idx / psfSize;
        int tmp = row + rr + (col+cc) * dims[0];

        #pragma omp atomic
        newImg[tmp] ++;
    }
}

void sample_t(double *img, uint32_t * imgIn, const size_t* dims, double * prior, int begin = 0, int finish = -1)
{
    int numel = dims[0] * dims[1];
    if (finish == -1) finish = numel;
    
    double s = 0;
    #pragma omp parallel for reduction(+:s)
    for (int i = begin; i < finish; i++) {
        gamma_distribution<double> randgam(imgIn[i] + prior[i]);
        img[i] = randgam(mts[omp_get_thread_num()]);
        s += img[i];
    }

    #pragma omp parallel for
    for (int i = begin; i < finish; i++) {
        img[i] /= s;
    }    
}

// Input -
// #1 - UINT32 Data. Either a low resolution image, or a 3XN array of SMLM data. 
//      If it's latter, row 1 is x, row 2 is y, row 3 is a index to psf array (see input #3)
// #2 - The image sample from last itration
// #3 - PSF. SxSxN 3D array, representing N different kinds of PSFs with 
//      varying localization accuracy. This is because in SMLM different
//      molecules are localized with different accuracy.
// #4 - Optional Prior. Either a scalar or a arrray/matrix. 
// Output -
// #1 - The sample. A double 2D array. It is normalized (summed to be 1.0).
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs != 3 && nrhs !=4) || nlhs != 1) {
        mexErrMsgTxt("lmsamplermex: Wrong numer of input/output arguments");
    }

    const mwSize d = mxGetNumberOfDimensions(prhs[2]);
    
    const mwSize *datadims = mxGetDimensions(prhs[0]);
    const mwSize *imgdims = mxGetDimensions(prhs[1]);
    const mwSize *psfdims = mxGetDimensions(prhs[2]);
    
    const mxUint32 * data = mxGetUint32s(prhs[0]);
    const mxDouble * img = mxGetDoubles(prhs[1]);
    const mxDouble * psf = mxGetDoubles(prhs[2]);

    mxDouble * prior = new mxDouble[mxGetNumberOfElements(prhs[1])];

    if (nrhs == 3) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), 1.0);
    } else if (mxIsScalar(prhs[3])) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), mxGetScalar(prhs[3]));
    } else {
        mxDouble * prior_in = mxGetDoubles(prhs[3]);
        copy(prior_in, prior_in + mxGetNumberOfElements(prhs[1]), prior);
    }

    plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
    mxDouble * newImg = mxGetDoubles(plhs[0]);
    uint32_t * tmpImg = new uint32_t[imgdims[0] * imgdims[1]];

    sample_r(tmpImg, img, imgdims, data, datadims[1], psf, psfdims[0]);
    sample_t(newImg, tmpImg, imgdims, prior);

    delete [] tmpImg;
    delete [] prior;
}

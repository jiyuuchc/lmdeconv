#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <random>
#include <omp.h>

using namespace std;

mt19937_64 mt;

static void init() __attribute__ ((constructor));
static void fini() __attribute__ ((destructor));

//mt19937_64 * mts;
//const uint32_t INTENSITY_FOR_POISS_APPROX = 5000;


void init() 
{
    random_device rd;
    /*mts = new mt19937_64[omp_get_max_threads()];
    for (int i = 0; i < omp_get_max_threads(); i++) 
        mts[i] = mt19937_64(rd());
    mexPrintf("LMSampler:init with %d threads\n",omp_get_max_threads());
     */
    mt = mt19937_64(rd());
}

void fini()
{
    //delete [] mts;
    //mexPrintf("SDS:exit\n");
}


void sample_r(uint32_t * newImg, const double* hiImg, const size_t *dims, const uint32_t *locs, const size_t numLocs, const double * psfs,const size_t psfSize, const double * psfZs, const size_t psfZSize, int begin=0, int finish=-1)
{
    uniform_real_distribution<> rand(0, 1);
    
    size_t psfHSize = (psfSize-1)/2;
    size_t psfZHsize = (psfZSize-1)/2;

    if (finish == -1) finish = numLocs;

    //mexPrintf("Sample_R \n");
    //mexPrintf("Img dimension %d - %d\n", dims[0], dims[1]);
    //mexPrintf("Psf dimension %d\n", psfSize);
    //mexPrintf("Num of locs %d\n", numLocs);
    
    fill(newImg, newImg + dims[0] * dims[1] * dims[2], 0);

    /* #pragma omp parallel for */
    for (int i = begin; i < finish; i ++) {
        double bins[psfSize*psfSize*psfZSize + 1];

        uint32_t row = locs[i*5+1] - psfHSize;
        uint32_t col = locs[i*5] - psfHSize;
        const double * psf = psfs + locs[i*5+2] * psfSize * psfSize;
        uint32_t slice = locs[i*5+3] - psfZHsize;
        const double * psfz = psfZs + locs[i*5+4] * psfZSize;

        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1] || slice < 0 || slice + psfZSize >= dims[2]) continue;

        unsigned int idx = 0, psfidx;
        bins[0] = 0;
        for (int z = slice; z < slice + psfZSize; z++) {
            psfidx = 0;
            for (int c = col ; c < col + psfSize; c++) {
                for (int r = row ; r < row + psfSize; r++) {
                    bins[idx + 1] = bins[idx] + hiImg[z * dims[0] * dims[1] + c * dims[0] + r] * psf[psfidx] * psfz[z-slice];
                    psfidx ++;
                    idx ++;
                }
            }
        }
        
        //double r = rand(mts[omp_get_thread_num()]) * bins[psfSize * psfSize * psfZSize];
        double r = rand(mt) * bins[psfSize * psfSize * psfZSize];
        
        idx = psfSize * psfSize * psfZSize / 2;
        unsigned int llimit = 1;
        unsigned int rlimit = psfSize * psfSize * psfZSize;
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
        int zz = idx / (psfSize * psfSize);
        idx = idx - zz * psfSize * psfSize;
        int rr = idx % psfSize;
        int cc = idx / psfSize;
        int tmp = (zz + slice) * dims[0]*dims[1] + row + rr + (col+cc) * dims[0];

        #pragma omp atomic
        newImg[tmp] ++;
    }
}

void sample_t(double *img, uint32_t * imgIn, const size_t* dims, double * prior, int begin = 0, int finish = -1)
{
    int numel = dims[0] * dims[1] * dims[2];
    if (finish == -1) finish = numel;
    
    double s = 0;
    /*#pragma omp parallel for reduction(+:s)*/
    for (int i = begin; i < finish; i++) {
        gamma_distribution<double> randgam(imgIn[i] + prior[i]);
        //img[i] = randgam(mts[omp_get_thread_num()]);
        img[i] = randgam(mt);
        s += img[i];
    }

    /*#pragma omp parallel for*/
    for (int i = begin; i < finish; i++) {
        img[i] /= s;
    }    
}

// Input -
// #1 - UINT32 Data. 5XN array of SMLM data. 
// #2 - The image sample from last itration
// #3 - xyPSF. SxSxN 3D array, representing N different kinds of PSFs with varying localization accuracy. 
// #4 - zPSF, 2D array
// #5 - Optional Prior. Either a scalar or a arrray/matrix. 
// Output -
// #1 - The sample. A double 3D array. It is normalized (summed to be 1.0).
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs != 4 && nrhs != 5) || nlhs != 1) {
        mexErrMsgTxt("lmsampler3dmex: Wrong numer of input/output arguments");
    }

    const mwSize imgND = mxGetNumberOfDimensions(prhs[1]); // image dimension
    if (imgND != 3) {
        mexErrMsgTxt("lmsampler3dmex: image is not 3d");
    }

    const mwSize *datadims = mxGetDimensions(prhs[0]);
    const mwSize *imgdims = mxGetDimensions(prhs[1]);
    const mwSize *psfxydims = mxGetDimensions(prhs[2]);
    const mwSize *psfzdims = mxGetDimensions(prhs[3]);
    
    const mxUint32 * data = mxGetUint32s(prhs[0]);
    const mxDouble * img = mxGetDoubles(prhs[1]);
    const mxDouble * psfxy = mxGetDoubles(prhs[2]);
    const mxDouble * psfz = mxGetDoubles(prhs[3]);

    
    mxDouble * prior = new mxDouble[mxGetNumberOfElements(prhs[1])];
    if (nrhs == 4) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), 1.0);
    } else if (mxIsScalar(prhs[4])) {
        fill(prior, prior + mxGetNumberOfElements(prhs[1]), mxGetScalar(prhs[4]));
    } else {
        mxDouble * prior_in = mxGetDoubles(prhs[4]);
        copy(prior_in, prior_in + mxGetNumberOfElements(prhs[1]), prior);
    }

    plhs[0] = mxCreateNumericArray(imgND, imgdims, mxDOUBLE_CLASS, mxREAL);
    mxDouble * newImg = mxGetDoubles(plhs[0]);
    uint32_t * tmpImg = new uint32_t[mxGetNumberOfElements(plhs[0])];

    sample_r(tmpImg, img, imgdims, data, datadims[1], psfxy, psfxydims[0], psfz, psfzdims[0]);
    sample_t(newImg, tmpImg, imgdims, prior);

    delete [] tmpImg;
    delete [] prior;
}

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
        mts[i] = mt19937_64(rd());
    mexPrintf("SDS:init with %d threads\n",omp_get_max_threads());
}

void fini()
{
    delete [] mts;
    mexPrintf("SDS:exit\n");
}

void sas(uint32_t* newImg, const uint32_t * img, const size_t * dims, const uint32_t scale, int begin = 0, int finish = -1) 
{
    uniform_int_distribution<uint32_t> randi(0, scale * scale-1);
    size_t numel = dims[0] * dims[1];
    
    if (finish == -1) finish = numel;
    
    #pragma omp parallel for
    for (int i = begin; i < finish; i ++) {
        uint32_t counts[scale*scale] = {};
        for (uint32_t c = 0; c < img[i]; c++) {
            counts[randi(mts[omp_get_thread_num()])] ++;
        }
        int idx = 0;
        int x0 = (int)(i / dims[0]) * scale;
        int y0 = (i % dims[0]) * scale;
        for (int x = 0; x < scale; x++) {
            for (int y = 0; y < scale; y++) {
                newImg[(x0 + x) * dims[0] * scale + y0 + y] = counts[idx++];
            }
        }
    }
}

// the version for image analysis
void sample_r(uint32_t* newImg, const double* hiImg, const size_t *dims, const uint32_t *loImg, const double * psf, const size_t psfSize, int beginCol=0, int endCol=-1)
{
    uniform_real_distribution<> rand(0, 1);
    
    size_t numel = dims[0] * dims[1];
    size_t psfHSize = (psfSize - 1) / 2;
    
    if (endCol == -1) endCol = dims[1];
    fill(newImg, newImg + numel, 0);

    #pragma omp parallel for
    for (int i = beginCol * dims[0]; i < endCol * dims[0]; i ++) {

        double bins[psfSize*psfSize + 1];
        bins[0] = 0;
    
        if (loImg[i] == 0) continue;

        size_t row = i % dims[0] - psfHSize;
        size_t col = i / dims[0] - psfHSize;
        
        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1]) continue;

        unsigned int idx = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx + 1] = bins[idx] + hiImg[c * dims[0] + r] * psf[idx];
                idx ++;
            }
        }
        
//        if (loImg[i] < INTENSITY_FOR_POISS_APPROX) {
            for (int k = 0; k < loImg[i]; k++) {
                double r = rand(mts[omp_get_thread_num()]) * bins[psfSize * psfSize];
                unsigned int idx = psfSize * psfSize / 2;
                unsigned int llimit = 1;
                unsigned int rlimit = psfSize * psfSize;
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
//            }
/*        } else { // for very bright pixles we do a poisson-distribution-based approximation
            idx = 0;
            for (int c  = col;c < col + psfSize; c++ ) {
                for (int r = row; r < row + psfSize; r++) {
                    poisson_distribution<uint32_t> randpoiss((bins[idx+1] - bins[idx]) * loImg[i] / bins[psfSize * psfSize]);
                    uint32_t d = randpoiss(mts[omp_get_thread_num()]);
                    int tmp = c * dims[0] + r;
                    idx ++;
                    #pragma omp atomic
                    newImg[tmp] += d;
                }
            } */
        }
    }
}

//the version for SMLM analysis
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
        
        uint32_t row = locs[i*3+1] - psfHSize;
        uint32_t col = locs[i*3] - psfHSize;
        const double * psf = psfs + locs[i*3+2] * psfSize * psfSize;
        
        if (row < 0 || col < 0 || row + psfSize >= dims[0] || col + psfSize >= dims[1]) continue;

        unsigned int idx = 0;
        bins[0] = 0;
        for (int c = col ; c < col + psfSize; c++) {
            for (int r = row ; r < row + psfSize; r++) {
                bins[idx + 1] = bins[idx] + hiImg[c * dims[0] + r] * psf[idx];
                idx ++;
            }
        }
        
        double r = rand(mts[omp_get_thread_num()]) * bins[psfSize * psfSize];
        
        idx = psfSize * psfSize / 2;
        unsigned int llimit = 1;
        unsigned int rlimit = psfSize * psfSize;
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

void sample_t(double *img, uint32_t * imgIn, const size_t* dims, int begin = 0, int finish = -1)
{
    int numel = dims[0] * dims[1];
    if (finish == -1) finish = numel;
    
    double s = 0;
    #pragma omp parallel for reduction(+:s)
    for (int i = begin; i < finish; i++) {
        gamma_distribution<double> randgam(imgIn[i] + 1);
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
// #3 - PSF. For image analysis, it's 2D double array. For SMLM, it's a SxSxN 3D array, representing
//      N different kinds of PSFs with varying localization accuracy. This is because in SMLM different
//      molecules are localized with different accuracy.
// Output -
// #1 - The sample. A double 2D array. It is normalized (summed to be 1.0).
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
        int scale = imgdims[0] / datadims[0];

        plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
        mxDouble * newImg = mxGetDoubles(plhs[0]);
        uint32_t * tmpImg1 = new uint32_t[imgdims[0] * imgdims[1]];
        uint32_t * tmpImg2 = new uint32_t[imgdims[0] * imgdims[1]];
	//uint32_t tmpImg1[imgdims[0] * imgdims[1]];
	//uint32_t tmpImg2[imgdims[0] * imgdims[1]];

        sas(tmpImg1, data, datadims, scale);
        sample_r(tmpImg2, img, imgdims, tmpImg1, psf, psfdims[0]);
        sample_t(newImg, tmpImg2, imgdims);

        delete [] tmpImg1;
        delete [] tmpImg2;
    } else { // else psf is 3D arrray, input data is single molecule localizations
        plhs[0] = mxCreateDoubleMatrix(imgdims[0], imgdims[1], mxREAL);
        mxDouble * newImg = mxGetDoubles(plhs[0]);
        uint32_t * tmpImg = new uint32_t[imgdims[0] * imgdims[1]];

        sample_r(tmpImg, img, imgdims, data, datadims[1], psf, psfdims[0]);
        sample_t(newImg, tmpImg, imgdims);
        
        delete [] tmpImg;
    }
}

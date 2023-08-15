/* Copyright (c) 2009-2013 Kyle Gorman
 *
 * Permission is hereby granted, free of charge, to any person obtaining a 
 * copy of this software and associated documentation files (the 
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to 
 * permit persons to whom the Software is furnished to do so, subject to 
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * swipe.c: primary functions
 * Kyle Gorman <gormanky@ohsu.edu>
 */

#define VNUM    1.5 // current version

#include <eigen3/Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <climits>
#include <cstdbool>
#include <cstring>
#include <ctime>
#include <vector>

#include <fftw3.h>   // http://www.fftw.org/

#include "swipe.h"


// #define DERBS    .01 
#define POLYV    (1./768.) //  1 / 12 / 64 = 1 / 768
// #define DLOG2P   (1./25.) // 1/25

// feel free to change these defaults
#define ST       .1
#define DT       .1
#define MIN      4.0
#define MAX      43.0

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef NAN
    #define NAN sqrt(-1.)
#endif

#ifndef isnan
int isnan(double x) { 
    return(x != x);
}
#endif

#ifndef log2
// a base-2 log function
double log2(double x) { 
    return log(x) / log(2.);
}
#endif

#ifndef round
// rounds a double to the nearest integer value
double round(double x) { 
    return(x >= 0. ? floor(x + .5) : floor(x - .5));
}
#endif

// converts from hertz to Mel frequency
double hz2mel(double hz) { 
    return(1127.01048 * log(1. + hz / 700.));
}

// converts from hertz to ERBs
double hz2erb(double hz) { 
    return(21.4 * log10(1. + hz / 229.));
}

// converts from ERBs to hertz 
double erb2hz(double erb) { 
    return((pow(10, erb / 21.4) - 1.) * 229.);
}

// a silly function that treats NaNs as 0.
double fixnan(double x) { 
    return(isnan(x) ? 0. : x);
}

// MAC WISDOM
#ifdef __APPLE__
const char* wisdom = "(fftw-3.3.10 fftw_wisdom #x4be12fff #x7b2df9b2 #xa5975329 #x385b0041\n"
"  (fftw_codelet_hc2cf_8 0 #x1040 #x1040 #x0 #x46284a7e #x5d7c8ac8 #x67a40edd #x701ec2f6)\n"
"  (fftw_codelet_hc2cf_4 0 #x1040 #x1040 #x0 #xf400d043 #xadbcc975 #x4bac40cb #x7833f0aa)\n"
"  (fftw_rdft_vrank_geq1_register 0 #x1040 #x1040 #x0 #x082776fc #x4b9b1d15 #xf8f6c54e #xed4e03fe)\n"
"  (fftw_codelet_r2cf_8 2 #x1040 #x1040 #x0 #xfeb8a7d9 #x0bb74e7a #x59436ea4 #x21d5fb09)\n"
"  (fftw_codelet_r2cf_8 2 #x1040 #x1040 #x0 #x5a87b94f #x630d3e77 #xa2e2443a #x3b0bc7d0)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #xa3e07b96 #x4a81ecf1 #xc26d608c #x8c15863d)\n"
"  (fftw_codelet_hc2cf2_4 0 #x1040 #x1040 #x0 #xaf51018d #x49bad7da #x8193787f #x1bb50304)\n"
"  (fftw_codelet_r2cfII_16 0 #x1040 #x1040 #x0 #x5c9657a7 #xb583909a #x3ee7b98a #x923ef24f)\n"
"  (fftw_codelet_r2cf_64 0 #x1040 #x1040 #x0 #x12a44a20 #x5dd60b71 #x21e0a37e #x23d64023)\n"
"  (fftw_codelet_r2cf_64 2 #x1040 #x1040 #x0 #x5ce6486a #xed9080e3 #x3bb7119c #xa780bec1)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #xd68c1b0d #xdec5fb3c #xea68bdea #xdd20ce78)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #xdc2d71d8 #xa8de4c88 #x8e4960e8 #x53d2af90)\n"
"  (fftw_codelet_hf_16 0 #x1040 #x1040 #x0 #x27839d50 #x0aae9d24 #x722b47c5 #x1b0ecefd)\n"
"  (fftw_codelet_r2cf_16 0 #x1040 #x1040 #x0 #x70ef7d20 #x19ebf46c #xd621f90d #x79951d6a)\n"
"  (fftw_codelet_r2cfII_8 2 #x1040 #x1040 #x0 #x1022d070 #x988061ba #xcd44ef4e #x6145c585)\n"
"  (fftw_codelet_hc2cf_8 0 #x1040 #x1040 #x0 #x4e9725b5 #x248e8217 #x00b1fddf #x59463c98)\n"
"  (fftw_codelet_r2cfII_8 2 #x1040 #x1040 #x0 #xc03df26b #xacaf3f9b #x5578608b #x5a1acfad)\n"
"  (fftw_codelet_r2cfII_8 2 #x1040 #x1040 #x0 #xdcb3eadc #x7aae4296 #x086ff51f #x41e592f6)\n"
"  (fftw_codelet_r2cf_64 0 #x1040 #x1040 #x0 #x1efbf9b7 #x4b507218 #x5daae168 #x61b8dd12)\n"
"  (fftw_codelet_r2cf_64 0 #x1040 #x1040 #x0 #xc9220c71 #x2b3a4a2a #xb4c31cac #x09d63773)\n"
"  (fftw_codelet_r2cf_4 2 #x1040 #x1040 #x0 #xb664f3c3 #xa36f4781 #x53e13a01 #x47de2854)\n"
"  (fftw_codelet_hc2cf_8 0 #x1040 #x1040 #x0 #x650927ce #xad2953b7 #x0d47a3c2 #x4d06139d)\n"
"  (fftw_codelet_r2cf_8 2 #x1040 #x1040 #x0 #xd57d4c22 #x4e1c4494 #x99472055 #x063b74e7)\n"
"  (fftw_codelet_hf_16 0 #x1040 #x1040 #x0 #xbe5138c4 #xc409e5b3 #xf0c4aa0b #x547d24ba)\n"
"  (fftw_codelet_hc2cf_4 0 #x1040 #x1040 #x0 #xcf2dc8f6 #xeaf35734 #x85e24035 #x769f96cf)\n"
"  (fftw_codelet_r2cf_16 0 #x1040 #x1040 #x0 #xafb941dc #xe9671516 #x1f560985 #x21b6d110)\n"
"  (fftw_codelet_r2cf_16 0 #x1040 #x1040 #x0 #x4de25858 #x88fbf594 #x99ddcdff #xcb89de26)\n"
"  (fftw_rdft_vrank_geq1_register 0 #x1040 #x1040 #x0 #xaf98a4e0 #xfcddea34 #x980b2d82 #xa57cca64)\n"
"  (fftw_codelet_r2cf_8 2 #x1040 #x1040 #x0 #xf34e66e1 #xe08b968c #xc34fa4dd #xe516a4e0)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #x5409a8a1 #x46be8189 #xad3795a2 #x77d8c9b9)\n"
"  (fftw_codelet_r2cfII_16 0 #x1040 #x1040 #x0 #xddb07789 #x5207474f #x82896f70 #xb790ff89)\n"
"  (fftw_codelet_r2cfII_4 2 #x1040 #x1040 #x0 #x184802f9 #xbba93e31 #x22f14020 #x93936099)\n"
"  (fftw_codelet_r2cf_16 0 #x1040 #x1040 #x0 #xee6fda48 #xf507a36d #x8f1189f7 #x16a0f19c)\n"
"  (fftw_codelet_r2cfII_4 2 #x1040 #x1040 #x0 #x90b835a1 #x0a5899b6 #x161bfa35 #x34436abd)\n"
"  (fftw_codelet_r2cf_32 2 #x1040 #x1040 #x0 #x1d2f3122 #x498102e6 #x68c7333f #x5ad939cc)\n"
"  (fftw_codelet_r2cf_4 2 #x1040 #x1040 #x0 #xafe2dd13 #x71b28566 #xf8fe01b0 #x581b6dd3)\n"
"  (fftw_rdft_vrank_geq1_register 0 #x1040 #x1040 #x0 #x8035d568 #x3014e232 #x12fd671f #x32ed03c0)\n"
"  (fftw_codelet_r2cf_8 2 #x1040 #x1040 #x0 #x4220eee2 #x68deaeed #x5efa6627 #x6f600652)\n"
"  (fftw_codelet_r2cf_4 2 #x1040 #x1040 #x0 #x20ce64c7 #xef0bd268 #x37826661 #xd29d0ed7)\n"
"  (fftw_codelet_hf_16 0 #x1040 #x1040 #x0 #xbe5a1e7f #xa439dc17 #xa27226b1 #x80295e83)\n"
"  (fftw_codelet_r2cf_4 2 #x1040 #x1040 #x0 #x29418a59 #xbb5e0c2f #xb3f1dadf #x2bfce015)\n"
"  (fftw_rdft_vrank_geq1_register 0 #x1040 #x1040 #x0 #x4580130f #xc7c5663e #xd6f21992 #x87ff4375)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #xde749c1e #x3e70a672 #xe8059160 #x2e308709)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #xef27b9b7 #xf73d2d36 #x7dbc85fb #x37bcd2da)\n"
"  (fftw_codelet_r2cf_16 2 #x1040 #x1040 #x0 #xc2956c83 #xf9e2390d #x332a3774 #x860920ff)\n"
"  (fftw_codelet_r2cf_64 0 #x1040 #x1040 #x0 #x9b2aa334 #x4e70617a #xcc6654d2 #x3d5e8a48)\n"
"  (fftw_codelet_r2cfII_8 2 #x1040 #x1040 #x0 #xdb9741c2 #xd78299ac #xa004c96e #xfeb498ab)\n"
"  (fftw_codelet_r2cfII_4 2 #x1040 #x1040 #x0 #xea121310 #x7812288f #x1cf96287 #xc6e9a81b)\n"
"  (fftw_codelet_hc2cf_8 0 #x1040 #x1040 #x0 #x587138ca #xa1bbe0b0 #x4572c05c #xe7acc963)\n"
"  (fftw_codelet_hf_16 0 #x1040 #x1040 #x0 #x361f62ab #x997c0f61 #x54a64057 #x1f35aef5)\n"
"  (fftw_codelet_r2cf_16 0 #x1040 #x1040 #x0 #x4db470f5 #x6d9348ae #x682f66b3 #xe516c1d4)\n"
"  (fftw_rdft_vrank_geq1_register 1 #x1040 #x1040 #x0 #x25724df0 #xba9bb091 #xb4ca44ff #xa583d488)\n"
")\n"
"";
#endif

// UNIX WISDOM
#ifdef __linux__
const char* wisdom = "(fftw-3.3.8 fftw_wisdom #x458a31c8 #x92381c4c #x4f974889 #xcd46f97e\n"
"  (fftw_codelet_n1fv_128_avx 0 #x31bff #x31bff #x0 #x37647671 #xa231fde8 #xcd1493d8 #x9f0a3dac)\n"
"  (fftw_codelet_hc2cfdftv_16_avx 0 #x31bff #x31bff #x0 #x7366d115 #xbc457627 #x02c92fe0 #xa85a47d3)\n"
"  (fftw_codelet_hc2cfdftv_32_avx 0 #x31bff #x31bff #x0 #xdc38cc14 #x0d1c6e3c #xc0ae48bf #xa0e0eb5a)\n"
"  (fftw_codelet_r2cfII_16 2 #x31bff #x31bff #x0 #x9fa443f9 #x0132a381 #xa9894cfa #x490e1490)\n"
"  (fftw_codelet_hc2cfdftv_2_avx 0 #x31bff #x31bff #x0 #x7da7817a #xbe20e7cd #x001a0a1b #x7f9f0731)\n"
"  (fftw_codelet_n1fv_16_avx 0 #x31bff #x31bff #x0 #x083037ce #x1237799f #x394f37a1 #x6cf21ac5)\n"
"  (fftw_codelet_r2cf_16 2 #x31bff #x31bff #x0 #x57f9650f #x9a38745d #x7c332073 #x7e66aa89)\n"
"  (fftw_codelet_n2fv_64_avx 0 #x31bff #x31bff #x0 #x4fa8b2f3 #xd93af2ba #xec419f5a #xae14e9fd)\n"
"  (fftw_codelet_n1fv_16_avx 0 #x31bff #x31bff #x0 #xfaab6250 #x004c256b #x78db2567 #xb6804f43)\n"
"  (fftw_codelet_n1fv_16_avx 0 #x31bff #x31bff #x0 #xb6ffb6dd #x8b384252 #xae9f2d40 #x43cb879f)\n"
"  (fftw_codelet_hc2cfdftv_2_avx 0 #x31bff #x31bff #x0 #x429d1761 #xc79f9ecc #xdb982251 #x10f2a2eb)\n"
"  (fftw_codelet_n2fv_32_avx 0 #x31bff #x31bff #x0 #x21bf8747 #x375e1866 #x1b9ee1b5 #xa0a9b5c3)\n"
"  (fftw_codelet_t3fv_16_avx 0 #x31bff #x31bff #x0 #x791ccda1 #xb979152f #x871ea4a6 #x2e274828)\n"
"  (fftw_codelet_r2cfII_2 2 #x31bff #x31bff #x0 #xfcde69de #xb320f2f0 #xf9c513b9 #x0388aa6b)\n"
"  (fftw_codelet_r2cfII_2 2 #x31bff #x31bff #x0 #x326e1847 #x46d0d5f5 #xb2f15e22 #xf96e89af)\n"
"  (fftw_codelet_n1fv_32_avx 0 #x31bff #x31bff #x0 #x11e5ffc0 #xcdd528bc #x8da1fac2 #x4a17575d)\n"
"  (fftw_dft_vrank_geq1_register 0 #x31bff #x31bff #x0 #x730e35f0 #x80430f15 #xe928d6bf #xcc157970)\n"
"  (fftw_codelet_r2cfII_32 2 #x31bff #x31bff #x0 #x0d723485 #x067d5676 #x69a91896 #x7274711f)\n"
"  (fftw_codelet_r2cfII_2 2 #x31bff #x31bff #x0 #x6e5132eb #x80cf1f87 #x62d8f14f #x169bb56f)\n"
"  (fftw_codelet_r2cf_2 2 #x31bff #x31bff #x0 #x5e097b05 #x6e0d6250 #xbd8a700a #xb202d582)\n"
"  (fftw_codelet_r2cf_2 2 #x31bff #x31bff #x0 #x6f5d40e3 #x02501563 #x153b8f39 #xffdc1a66)\n"
"  (fftw_codelet_hc2cfdftv_32_avx 0 #x31bff #x31bff #x0 #x2b2224fd #xc6d31fc1 #xf08ddb32 #xd3726977)\n"
"  (fftw_codelet_hc2cfdftv_32_avx 0 #x31bff #x31bff #x0 #x4fc99822 #xf85b48fc #x7d0a0313 #x81834ad0)\n"
"  (fftw_codelet_hc2cfdftv_2_avx 0 #x31bff #x31bff #x0 #xb9ca5008 #xdaec3213 #x06c7b7c1 #x9db6000f)\n"
"  (fftw_codelet_r2cfII_32 2 #x31bff #x31bff #x0 #x11bed0d2 #x53e7db34 #xa88d8340 #x0a29b429)\n"
"  (fftw_codelet_r2cf_2 2 #x31bff #x31bff #x0 #x96022fcd #x29ef0675 #x1b46036f #x81e6cbd5)\n"
"  (fftw_codelet_r2cfII_32 2 #x31bff #x31bff #x0 #x9854524d #x45e23508 #xe087f9f8 #xdf099e07)\n"
")";
#endif


// TODO: Hardcoded primes may need more depending on use case
const int primes[169] = {1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
941,947,953,967,971,977,983,991,997};

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

extern "C" void splev_(double* t, int* n, double* c, int* nc, int* k, double* x, double* y, int* m, int* e, int* ier);

extern "C" void splder_(double* t, int* n, double* c, int* nc, int* k, int* nu, double* x, double* y, int* m, int* e, double* wrk, int* ier);

extern "C" void curfit_(int* iopt, int* m, double* x, double* y, double* w, double* xb, double*  xe, int* k, 
    double* s, int* nest, int* n, double* t, double* c, double* fp, double* wrk, int* lwrk, int* iwrk, int* ier);

extern "C" int dgels_(char* trans, int* m, int* n, int* nrhs,
                           double* A, int* lda, double* B, int* ldb,
                           double* work, int* lwork, int* info);



// given a vector of coefficients and a value for x, evaluate the polynomial
double polyval(const std::vector<double>& coefs, double val) {
    double sum = 0.;
    size_t i;
    for (i = 0; i < coefs.size(); i++)
        sum += coefs[i] * pow(val, coefs.size() - i - 1);
    return(sum);
}


// polynomial fitting with CLAPACK: solves poly(A, m) * X = B
std::vector<double> polyfit(std::vector<double>& A, std::vector<double>& B, int degree) {
    int info;
    degree++; // I find it intuitive this way...
    double* Ap = (double *)malloc(sizeof(double) * degree * A.size());
    int i, j;
    for (i = 0; i < degree; i++)
        for (j = 0; j < (int)A.size(); j++)
            Ap[i * A.size() + j] = pow(A[j], degree - i - 1); // mimics MATLAB
    std::vector<double> Bp(degree >= static_cast<int>(B.size()) ? degree : B.size());
    for (size_t i = 0; i < B.size(); i++)
        Bp[i] = B[i];
    i = 1; // nrhs, j is info
    j = A.size() + degree; // lwork
    double* work = (double*) malloc(sizeof(double) * j);

    int a_size = (int)A.size();
    int b_size = (int)B.size();

    char n = 'N';

    dgels_(&n, &a_size, &degree, &i, Ap, &b_size, Bp.data(), &degree, work, &j,
                                                              &info);
    free(Ap);
    free(work);
    if (info < 0) {
        fprintf(stderr, "LAPACK routine dgels() returned error: %d\n",
                                                                info);
        exit(EXIT_FAILURE);
    }
    return(Bp);
}

/* like bisectv(), but the minimum starting value is passed as an argument.
 * This is good for multiple bisection calls for forming a new vector when
 * the queries are a non-constant interval; but make sure to use bisectv()
 * the first time.
 */
int bilookv(std::vector<double>& yr_vector, double key, int lo) {
    int md;
    size_t hi = yr_vector.size();
    lo--;
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vector[md] > key)
            hi = md;
        else
            lo = md;
    }
    return(hi);
}

std::vector<double> fbspline(std::vector<double>& f, std::vector<double>& a, int k, std::vector<double>& ferbs) {
    int iopt = 0;
    int m = f.size();
    std::vector<double> w(m, 1);
    double s = 0.0;
    double xb = f[0];
    double xe = f.back();
    int nest = m+k+1 > 2*k+3 ? m+k+1 : 2*k+3;
    int n = nest;

    std::vector<double> t(nest);
    std::vector<int> iwrk(nest);
    std::vector<double> wrk(m*(k + 1) + nest*(7 + 3*k));

    int lwrk = m*(k + 1) + nest*(7 + 3*k);
    int ier;
    double fp;

    std::vector<double> c(nest);

    curfit_(&iopt, &m, f.data(), a.data(), w.data(), &xb, &xe, &k, &s, &nest, &n, t.data(), c.data(), &fp, wrk.data(), &lwrk, iwrk.data(), &ier);


    int nu = 0;
    int e = 0;
    ier = 0;
    m = ferbs.size();
    n = nest;
    int nc = c.size();

    std::vector<double> y(m, 0);

    splev_(t.data(), &n, c.data(), &nc, &k, ferbs.data(), y.data(), &m, &e, &ier);

    return y;
}

// find the bisection index of the std::vector<double> for key
int bisectv(const std::vector<double>& yr_vector, double key) {
    int md;
    int lo = 1;
    int hi = yr_vector.size();
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vector[md] > key)
            hi = md;
        else
            lo = md;
    }
    return(hi);
}

// a helper function for loudness() for individual fft slices
void La(Matrix& L, std::vector<double>& f, std::vector<double>& fERBs, fftw_plan plan, 
                            fftw_complex* fo, int w2, int hi, int i) {
    int j;
    fftw_execute(plan);
    std::vector<double> a(w2+1);

    for (j = 0; j < w2+1; j++) // this iterates over only the first half
        a[j] = sqrt(fo[j][0] * fo[j][0] + fo[j][1] * fo[j][1]); // spectrum in python

    // std::vector<double> a2 = spline(f, a); // a2 is now the result of the cubic spline
    std::vector<double> y = fbspline(f, a, 3, fERBs);

    for (j = 0; j < L.cols(); j++) {
        L(i,j) = fixnan(sqrt(y[j]));
    }
}

// a function for populating the loudness matrix with a signal x
Matrix loudness(std::vector<double>& x, std::vector<double>& fERBs, double nyquist, int w, int w2) {
    fftw_plan plan;
    int i, j, hi; 
    int offset = 0;
    double td = nyquist / w2; // this is equivalent to fstep
    // testing showed this configuration of fftw to be fastest
    double* fi = (double*) fftw_malloc(sizeof(double) * w); 
    fftw_complex* fo = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w);

    plan = fftw_plan_dft_r2c_1d(w, fi, fo, FFTW_ESTIMATE);
 
    std::vector<double> hann(w); // this defines the Hann[ing] window
    for (i = 0; i < w; i++) 
        hann[i] = .5 - (.5 * cos(2. * M_PI * ((double) i / w)));

    std::vector<double> f(w2+1);

    for (i = 0; i < (int)f.size(); i++) 
        f[i] = i * td;

    hi = bisectv(f, fERBs[0]); // all calls to La() will begin here

    Matrix L((int)(ceil((double) x.size() / w2) + 1), fERBs.size());

    for (j = 0; j < w2; j++) // left boundary case
        fi[j] = 0.; // more explicitly, 0. * hanna[j]
    for (; j < w; j++) 
        fi[j] = x[j-w2] * hann[j];
    La(L, f, fERBs, plan, fo, w2, hi, 0); 

    for (i = 1; i < L.rows() - 2; i++) {
        for (j=0; j < w; j++) 
            fi[j] = x[j + offset] * hann[j];
        La(L, f, fERBs, plan, fo, w2, hi, i); 
        offset += w2;
    }

    for (/* i = L.size() - 2; */; i < L.rows(); i++) { // right two boundary cases
        for (j = 0; j < (int)x.size() - offset; j++) // this dies at x.size() + w2
            fi[j] = x[j + offset] * hann[j];
        for (/* j = x.size() - offset */; j < w; j++) 
            fi[j] = 0.; // once again, 0. * hanna[j] 

        La(L, f, fERBs, plan, fo, w2, hi, i);
        offset += w2;
    } // now L is fully valued

    // L must now be normalized
    for (i = 0; i < L.rows(); i++) { 
        td = 0.; // td is the value of the normalization factor

        td = (L.row(i).array() * L.row(i).array()).sum();

        if (td != 0.) { // catches zero-division
            td = sqrt(td);
            L.row(i) /= td;
        } // otherwise, it is already 0.
    } 

    fftw_destroy_plan(plan); 
    fftw_free(fi); 
    fftw_free(fo); 
    return(L);
}


// populates the strength matrix using the loudness matrix
void Sadd(Matrix& S, const Matrix& L, std::vector<double>& fERBs, std::vector<double>& pci, std::vector<double>& mu, 
                                            int ps, double dt, 
                                            double nyquist2, int lo, 
                                            int psz, int w2) {
    int i, j, k;
    double t = 0.;
    double tp = 0.;
    double td;
    double dtp = w2 / nyquist2;

    Matrix Slocal(psz, L.rows());

    for (i = 0; i < Slocal.rows(); i++) {
        std::vector<double> q(fERBs.size());

        for (j = 0; j < (int)q.size(); j++) q[j] = fERBs[j] / pci[i];

        Eigen::VectorXd kernel = Eigen::VectorXd::Zero(fERBs.size()); // a zero-filled kernel vector

        int klast = 0;
        for (j = 0; primes[j] <= ps; j++) {    
            int prime = primes[j];

            double min = prime - .75;
            k = bilookv(q, min, klast);
            klast = k;

            for (; k < kernel.size(); k++) {
                td = fabs(q[k] - prime); 
                if (td < .25) // peaks
                    kernel[k] = cos(2. * M_PI * q[k]);
                else if (td < .75)  // valleys
                    kernel[k] += cos(2. * M_PI * q[k]) / 2.;
                else {
                    break;
                }
            }
        }

        td = 0.; 
        for (j = 0; j < kernel.size(); j++) {
            kernel[j] *= sqrt(1. / fERBs[j]); // applying the envelope
            if (kernel[j] > 0.) 
                td += kernel[j] * kernel[j];
        }
        td = sqrt(td); // now, td is the p=2 norm factor
        kernel /= td;

        Slocal.row(i) = L * kernel;
    } // Slocal is filled out; time to interpolate

    k = 0; 
    for (j = 0; j < S.cols(); j++) { // determine the interpolation params 
        td = t - tp; 
        while (td >= 0.) {
            k++;
            tp += dtp;
            td -= dtp;
        } // td now equals the time difference
        for (i = 0; i < psz; i++) {
            S(lo + i,j) += (Slocal(i,k) + (td * (Slocal(i,k) - Slocal(i,k-1))) / dtp) * mu[i];
        }
        t += dt;
    }
}

// helper function for populating the strength matrix on left boundary
void Sfirst(Matrix& S, std::vector<double>& x, std::vector<double>& pc, std::vector<double>& fERBs, std::vector<double>& d, 
                                           std::vector<int>& ws, int ps, 
                                           double nyquist, double nyquist2,
                                           double dt, int n) {
    size_t i; 
    int w2 = ws[n] / 2;
    std::vector<double> x_pad(x.size() + ws[n], 0);
    for (i = 0; i < x.size(); i++) {
        x_pad[i] = x[i];
    }

    Matrix L = loudness(x_pad, fERBs, nyquist, ws[n], w2);

    int lo = 0; // the start of Sfirst-specific code
    int hi = bisectv(d, 2.);
    int psz = hi - lo;
    std::vector<double> mu(psz);
    std::vector<double> pci(psz);
    for (i = 0; (int)i < hi; i++) {
        pci[i] = pc[i];
        mu[i] =  d[i]-1. <= 0 ? 1. : 1. - fabs(d[i]-1);
    } // end of Sfirst-specific code

    Sadd(S, L, fERBs, pci, mu, ps, dt, nyquist2, lo, psz, w2); 
}

// generic helper function for populating the strength matrix
void Snth(Matrix& S, std::vector<double>& x, std::vector<double>& pc, std::vector<double>& fERBs, std::vector<double>& d,
                              std::vector<int>& ws, int ps, double nyquist, 
                              double nyquist2, double dt, int n) {
    int i;
    int w2 = ws[n] / 2;
    std::vector<double> x_pad(x.size() + ws[n], 0);
    for (int i = 0; i < (int)x.size(); i++) {
        x_pad[i] = x[i];
    }
    Matrix L = loudness(x_pad, fERBs, nyquist, ws[n], w2);

    int lo = bisectv(d, n); // start of Snth-specific code
    int hi = bisectv(d, n + 2);
    int psz = hi - lo;
    std::vector<double> mu(psz);
    std::vector<double> pci(psz);
    int ti = 0;
    for (i = lo; i < hi; i++) {
        pci[ti] = pc[i];
        mu[ti] = 1. - fabs(d[i] - (n+1));
        ti++;
    }// end of Snth-specific code

    Sadd(S, L, fERBs, pci, mu, ps, dt, nyquist2, lo, psz, w2); 
}

// helper function for populating the strength matrix from the right boundary
void Slast(Matrix& S, std::vector<double>& x, std::vector<double>& pc, std::vector<double>& fERBs, std::vector<double>& d, 
                                          std::vector<int>& ws, int ps, 
                                          double nyquist, double nyquist2, 
                                          double dt, int n) {
    size_t i;
    int w2 = ws[n] / 2;
    std::vector<double> x_pad(x.size() + ws[n], 0);
    for (i = 0; i < x.size(); i++) {
        x_pad[i] = x[i];
    }

    Matrix L = loudness(x_pad, fERBs, nyquist, ws[n], w2);
    int lo = bisectv(d, n); // start of Slast-specific code
    size_t hi = d.size();
    int psz = hi - lo;
    std::vector<double> mu(psz);
    std::vector<double> pci(psz);
    int ti = 0;
    for (i = lo; i < hi; i++) {
        pci[ti] = pc[i];
        mu[ti] = d[i]-(n+1) >= 0 ? 1. : 1. - fabs(d[i] - (n+1));
        ti++;
    } // end of Slast-specific code

    Sadd(S, L, fERBs, pci, mu, ps, dt, nyquist2, lo, psz, w2); 
}

// performs polynomial tuning on the strength matrix to determine the pitch
PitchVector pitch(Matrix& S, std::vector<double>& pc, double st) {
    int i, j;
    int maxi = -1;
    int search = (int) round((log2(pc[2]) - log2(pc[0])) / POLYV + 1.);
    double nftc, maxv, log2pc;
    double tc2 = 1. / pc[1];
    std::vector<double> coefs;
    std::vector<double> s(3);
    std::vector<double> ntc(3);
    ntc[0] = ((1. / pc[0]) / tc2 - 1.) * 2. * M_PI; 
    ntc[1] = (tc2 / tc2 - 1.) * 2. * M_PI; 
    ntc[2] = ((1. / pc[2]) / tc2 - 1.) * 2. * M_PI;

    PitchVector results;
    results.strength = std::vector<double>(S.cols());
    results.pitch = std::vector<double>(S.cols());

    for (j = 0; j < S.cols(); j++) {
        maxv = SHRT_MIN;  
        for (i = 0; i < S.rows(); i++) {
            if (S(i,j) > maxv) {
                maxv = S(i,j);
                maxi = i;
            }
        }
        if (maxv > st) { // make sure it's big enough
            log2pc = log2(pc[maxi - 1]);
            if (maxi == 0 || maxi == S.rows() - 1) {
                results.pitch[j] = pc[0];
                // results.pitch[j] = 0;
                results.strength[j] = maxv;
            } 
            else { // general case
                tc2 = 1. / pc[maxi]; 
                s[0] = S(maxi - 1,j);
                s[1] = S(maxi,j);
                s[2] = S(maxi + 1,j); 
                
                coefs = polyfit(ntc, s, 2); 
                maxv = SHRT_MIN; 
                
                for (i = 0; i < search; i++) { // check the nftc space
                    nftc = polyval(coefs, ((1. / pow(2, i * POLYV + 
                                   log2pc)) / tc2 - 1) * 2 * M_PI);
                    if (nftc > maxv) {
                        maxv = nftc;
                        maxi = i;
                    }
                } // now we've got the pitch numbers we need
                results.pitch[j] = pow(2, log2pc + (maxi * POLYV));
                results.strength[j] = maxv;
            }
        }
        else{
			results.pitch[j] = 0;
            results.strength[j] = maxv;
		}
    }

    return(results);
}

// primary utility function for each pitch extraction
PitchVector swipe(std::vector<double>& x, int32_t fs, int32_t step, double min, double max, double dlog2p, double derbs, double st) {
    int i; 
    double td = 0.;
    double dt = (double)step/fs;

    if (fftw_import_wisdom_from_string(wisdom) == 0) {
        printf("Failed to load wisdom!\n");
        printf("Will create wisdom using FFTW_ESTIMATE flag\n");
        // exit(EXIT_FAILURE);
    }

    double nyquist = fs / 2.;
    double nyquist2 = fs;
    double nyquist16 = fs * 8.;
    if (max > nyquist) {
        fprintf(stderr, "Max pitch exceeds Nyquist frequency...\n");
        exit(EXIT_FAILURE);
    }
    if (dt > nyquist2) {
        fprintf(stderr, "Max time step exceeds nyquist frequency\n");
        exit(EXIT_FAILURE);
    }


    std::vector<int> ws(round(log2((nyquist16) / min) -  
                                log2((nyquist16) / max)) + 1); 
    for (i = 0; i < (int) ws.size(); i++)
        ws[i] = pow(2, round(log2(nyquist16 / min))) / pow(2, i);

    std::vector<double> pc = std::vector<double>(ceil((log2(max) - log2(min)) / dlog2p));

    std::vector<double> d(pc.size());

    for (i = (int)pc.size() - 1; i >= 0; i--) { 
        td = log2(min) + (i * dlog2p);
        pc[i] = pow(2, td);
        d[i] = 1. + td - log2(nyquist16 / ws[0]); 
    } // td now equals log2(min)

    std::vector<double> fERBs(ceil((hz2erb(nyquist) - 
                               hz2erb(pow(2, td) / 4)) / derbs));
	// std::vector<double> fERBs(ceil((ceil(hz2erb(nyquist)) - hz2erb(pc[0]/4))/derbs));
    td = hz2erb(min / 4.);
    for (i = 0; i < (int)fERBs.size(); i++) 
        fERBs[i] = erb2hz(td + (i * derbs));



    int ps = floor(fERBs[fERBs.size() - 1] / pc[0] - .75);

    Matrix S((int)pc.size(), (int)ceil(((double) x.size() / nyquist2) / dt));
    S.setZero();

    Sfirst(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, dt, 0); 
    for (i = 1; i < (int)ws.size() - 1; i++) // S is updated inline here
        Snth(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, dt, i);
    // i is now (ws.size() - 1)
    Slast(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, dt, i);

    PitchVector p = pitch(S, pc, st); // find pitch using strength matrix

    fftw_cleanup();
    return(p);
}
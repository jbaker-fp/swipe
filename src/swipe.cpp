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

// #include <fftw3.h>   // http://www.fftw.org/

#include "pffft_double.h"

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
void La(Matrix& L, std::vector<double>& f, std::vector<double>& fERBs, PFFFTD_Setup* plan, 
                            double* fo, double* fi, int w2, int hi, int i) {
    int j;
    pffftd_transform_ordered(plan, fi, fo, NULL, PFFFT_FORWARD);
    std::vector<double> a(w2+1);

    for (j = 0; j < w2+1; j+=2) // this iterates over only the first half
        a[j] = sqrt(fo[j] * fo[j] + fo[j+1] * fo[j+1]); // spectrum in python

    // std::vector<double> a2 = spline(f, a); // a2 is now the result of the cubic spline
    std::vector<double> y = fbspline(f, a, 3, fERBs);

    for (j = 0; j < L.cols(); j++) {
        L(i,j) = fixnan(sqrt(y[j]));
    }
}

// a function for populating the loudness matrix with a signal x
Matrix loudness(std::vector<double>& x, std::vector<double>& fERBs, double nyquist, int w, int w2) {
    PFFFTD_Setup* plan = pffftd_new_setup(w, PFFFT_COMPLEX);
    int i, j, hi; 
    int offset = 0;
    double td = nyquist / w2; // this is equivalent to fstep

    // testing showed this configuration of fftw to be fastest
    double* fi = (double*) pffftd_aligned_malloc(sizeof(double) * w); 
    double* fo = (double*) pffftd_aligned_malloc(sizeof(double) * (w * 2));

    // plan = fftw_plan_dft_r2c_1d(w, fi, fo, FFTW_ESTIMATE);
 
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
    La(L, f, fERBs, plan, fo, fi, w2, hi, 0); 

    for (i = 1; i < L.rows() - 2; i++) {
        for (j=0; j < w; j++) 
            fi[j] = x[j + offset] * hann[j];
        La(L, f, fERBs, plan, fo, fi, w2, hi, i); 
        offset += w2;
    }

    for (/* i = L.size() - 2; */; i < L.rows(); i++) { // right two boundary cases
        for (j = 0; j < (int)x.size() - offset; j++) // this dies at x.size() + w2
            fi[j] = x[j + offset] * hann[j];
        for (/* j = x.size() - offset */; j < w; j++) 
            fi[j] = 0.; // once again, 0. * hanna[j] 

        La(L, f, fERBs, plan, fo, fi, w2, hi, i);
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

    pffftd_destroy_setup(plan); 
    free(fi); 
    free(fo); 
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
PitchVector swipe(std::vector<double>& x, int32_t fs, int32_t step, double min, double max, double dlog2p, double derbs, double st, const std::string& wisdom) {
    int i; 
    double td = 0.;
    double dt = (double)step/fs;

    if (wisdom != "") {
        if (fftw_import_wisdom_from_string(wisdom.data()) == 0) {
            printf("Failed to load wisdom!\n");
            printf("Will create wisdom using FFTW_ESTIMATE flag\n");
        }
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
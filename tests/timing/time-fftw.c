// Build with:  gcc time-fftw.c -std=c99 -lfftw3
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char**argv) {
    if (argc < 2) {
        printf("Need size argument.\n");
        return 1;
    }

    int n = atoi(argv[1]);

    int numIters = 1000 * 10;
    //printf("Running dft %d times for size %d.\n", numIters, n);

    fftw_complex *in, *out;

    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    for (int j=0; j<n; j++) {
            in[j] = 17;
   }
    p = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
   for (int i=0; i<numIters; i++) {
        fftw_execute(p);
        /*
        double s = (1 / sqrt ((double)n));
        for (int j=0; j<n; j++) {
            out[j] *= s;
        }
        */
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}


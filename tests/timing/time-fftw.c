#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>

int main() {
    int n = 1024;

    fftw_complex *in, *out;

    fftw_plan p;

    int repeats = 10000;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    for (int i=0; i<repeats; i++) {
        p = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        for (int j=0; j<n; j++) {
            in[j] = i;
        }
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

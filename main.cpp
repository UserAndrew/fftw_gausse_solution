#include <complex>
#include <fftw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>

int main()
{
    const int N = 1024;
    double *x_ = new double[N];
    double *p_ = new double[N];
    fftw_complex *in = new fftw_complex[N];
    fftw_complex *out = new fftw_complex[N];
    double Xmin = -10.;
    double Xmax = 10.;
    double h = (Xmax - Xmin)/N;
    double hp = 2.0*M_PI/(h*N);

    for (int i = 0; i < N/2; i++)
    {
        p_[i]=hp*(double)i;
    }
    for (int i = N/2; i < N; i++)
    {
        p_[i]= -hp*(double)(N-i);
    }
    int i = 0;
    int count = 0;
    double X = Xmin;
    while(X < Xmax)
    {
        x_[i] = X;
        in[i][0] = exp(((-1.)*X*X)/2.);
        in[i][1] = 0.;
        X += h;
        i++;
        count++;
    }
    std::cout << count << std::endl;
    std::ofstream _out("fft_gausse_test.txt");
    _out.precision(4);
    for(int i=0; i < N; ++i)
    {
        _out<<std::fixed<<in[i][0]<<"\t"<<x_[i]<<std::endl;
    }


    fftw_plan my_plan;
    my_plan=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(my_plan);

    for (int i = 0; i < N; i++)
    {
        double re = out[i][0];
        double im = out[i][1];
        double phase = -double(i)*M_PI;
        out[i][0] = re*cos(phase) - im*sin(phase);
        out[i][1] = re*sin(phase) + im*cos(phase);
    }

    std::ofstream _out_("fft_gausse_our.txt");
    _out_.precision(4);
    for (int i = N/2; i < N; i++)
    {
        _out_ << std::fixed << p_[i] << "\t"<< out[i][0]*h/(sqrt(2*M_PI)) <<
                 "\t"<< out[i][1]*h/(sqrt(2*M_PI)) << std::endl;
    }
    for (int i = 0; i < N/2; i++)
    {
        _out_ << std::fixed << p_[i] << "\t"<< out[i][0]*h/(sqrt(2*M_PI)) <<
                 "\t"<< out[i][1]*h/(sqrt(2*M_PI)) << std::endl;
    }

    fftw_destroy_plan(my_plan);

    delete [] out;
    delete [] in;
    delete [] p_;
    delete [] x_;

    return 0;
}

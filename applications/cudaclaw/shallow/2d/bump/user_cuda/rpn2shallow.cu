#include "../bump_user.h"

#include <fc2d_cudaclaw.h>

#include <fc2d_cudaclaw_check.h>

__constant__ double s_grav;


void setprob_cuda()
{
    double grav;
    FILE *f = fopen("setprob.data","r");
    fscanf(f,"%lf",&grav);
    fclose(f);
    
    CHECK(cudaMemcpyToSymbol(s_grav,  &grav, sizeof(double)));
}

__device__ void bump_compute_speeds(int idir, int meqn, int mwaves, int maux,
                                    double ql[], double  qr[],
                                    double auxl[], double auxr[],
                                    double s[])
{
    int mu = 1+idir;
    //int mv = 2-idir;    

    double h = (qr[0] + ql[0])/2.0;
    double hsqrtr = sqrt(qr[0]);
    double hsqrtl = sqrt(ql[0]);
    double hsq2 = hsqrtl + hsqrtr;

    double u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
    // double v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
    double a = sqrt(s_grav*h);    

    s[0] = u-a;
    s[1] = u;
    s[2] = u+a;
}

__device__ cudaclaw_cuda_speeds_t bump_speeds = bump_compute_speeds;

void bump_assign_speeds(cudaclaw_cuda_speeds_t *speeds)
{
    cudaError_t ce = cudaMemcpyFromSymbol(speeds, bump_speeds, sizeof(cudaclaw_cuda_speeds_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (bump_compute_speeds): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void bump_rpn2shallow(int idir, int meqn, int mwaves, 
                                 int maux, double ql[], double qr[], 
                                 double auxl[], double auxr[],
                                 double wave[], double s[], 
                                 double amdq[], double apdq[])
{
    //assert(mwaves == 3);
    //assert(meqn == 3);

    int mu = 1+idir;
    int mv = 2-idir;    

    /* This should be true ... */
    //assert(qr[0] > 0);
    //assert(ql[0] > 0);

	double h = (qr[0] + ql[0])/2.0;
	double hsqrtr = sqrt(qr[0]);
	double hsqrtl = sqrt(ql[0]);
	double hsq2 = hsqrtl + hsqrtr;

	double u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
	double v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
	double a = sqrt(s_grav*h);

	double delta1 = qr[0 ] - ql[0 ];
	double delta2 = qr[mu] - ql[mu];
	double delta3 = qr[mv] - ql[mv];
	double a1 = ( (u+a)*delta1 - delta2 )/(2*a);
	double a2 = -v*delta1 + delta3;
	double a3 = (-(u-a)*delta1 + delta2 )/(2*a);
    
    wave[0]  = a1;
    wave[mu] = a1*(u-a);
    wave[mv] = a1*v;
    s[0] = u-a;
    
    wave[0  + meqn] = 0.0;
    wave[mu + meqn] = 0.0;
    wave[mv + meqn] = a2;
    s[1] = u;
    
    wave[0  + meqn*2] = a3;
    wave[mu + meqn*2] = a3*(u+a);
    wave[mv + meqn*2] = a3*v;
    s[2] = u+a;
    
    double smin[3],smax[3];
    for(int mw = 0; mw < mwaves; mw++)
    {
        smin[mw] = (s[mw] < 0) ? s[mw] : 0.;
        smax[mw] = (s[mw] >= 0) ? s[mw] : 0.;
    }
    for(int mq = 0; mq < meqn; mq++)
    {
        /* Loop-unrolling! loop over mwaves=3*/
        amdq[mq]  = smin[0]*wave[mq];
        amdq[mq] += smin[1]*wave[meqn + mq];
        amdq[mq] += smin[2]*wave[2*meqn + mq];

        apdq[mq]  = smax[0]*wave[mq];
        apdq[mq] += smax[1]*wave[meqn + mq];
        apdq[mq] += smax[2]*wave[2*meqn + mq];
   }    

}

__device__ cudaclaw_cuda_rpn2_t bump_rpn2 = bump_rpn2shallow;

void bump_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, bump_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (bump_rpn2shallow): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void bump_rpt2shallow(int idir, int meqn, int mwaves, int maux,
                                 double ql[], double qr[], 
                                 double aux1[], double aux2[], double aux3[],
                                 int imp, double asdq[],
                                 double bmasdq[], double bpasdq[])
{
    int mu = 1+idir;
    int mv = 2-idir;    

    double h = (qr[0] + ql[0])/2.0;
    double hsqrtr = sqrt(qr[0]);
    double hsqrtl = sqrt(ql[0]);
    double hsq2 = hsqrtl + hsqrtr;

    double u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
    double v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
    double a = sqrt(s_grav*h);


    double alpha1 = ((v+a)*asdq[0] - asdq[mv])/(2.0*a);
    double alpha2 = asdq[mu] - u*asdq[0];
    double alpha3 = (-(v-a)*asdq[0] + asdq[mv])/(2.0*a);

    double waveb[9], sb[3];
    waveb[0]  = alpha1;
    waveb[mu] = alpha1*u;
    waveb[mv] = alpha1*(v-a);
    sb[0] = v - a;

    waveb[meqn + 0]  = 0.0;
    waveb[meqn + mu] = alpha2;
    waveb[meqn + mv] = 0.0;
    sb[1] = v;

    waveb[2*meqn + 0]  = alpha3;
    waveb[2*meqn + mu] = alpha3*u;
    waveb[2*meqn + mv] = alpha3*(v+a);
    sb[2] = v + a;

    double smin[3], smax[3];
    for(int mw = 0; mw < mwaves; mw++)
    {
        smin[mw] = (sb[mw] < 0) ? sb[mw] : 0.0;
        smax[mw] = (sb[mw] >= 0) ? sb[mw] : 0.0;        
    }

    for(int mq = 0; mq < meqn; mq++)
    {        
        /* Loop-unrolling! loop over mwaves=3*/
        bmasdq[mq]  = smin[0]*waveb[mq];
        bmasdq[mq] += smin[1]*waveb[meqn + mq];
        bmasdq[mq] += smin[2]*waveb[2*meqn + mq];

        bpasdq[mq]  = smax[0]*waveb[mq];
        bpasdq[mq] += smax[1]*waveb[meqn + mq];
        bpasdq[mq] += smax[2]*waveb[2*meqn + mq];
    }
}


__device__ cudaclaw_cuda_rpt2_t bump_rpt2 = bump_rpt2shallow;

void bump_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, bump_rpt2, 
                                          sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (bump_rpt2shallow): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


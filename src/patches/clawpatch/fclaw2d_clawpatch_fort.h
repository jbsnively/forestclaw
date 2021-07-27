/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FCLAW2D_CLAWPATCH_FORT_H
#define FCLAW2D_CLAWPATCH_FORT_H

#include <fclaw2d_defs.h>

#if FCLAW2D_PATCHDIM == 2
#include "fclaw2d_clawpatch_fort2.h"
#else
#include "fclaw2d_clawpatch_fort3.h"
#endif


#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

struct fclaw2d_patch_transform_data;  /* Should be replaced by long int?  */

/* Functions defined here are implemented in individual solvers (clawpack 4.6 and 
   clawpack 5.0) */

#if 0
/* Fix syntax highlighting */
#endif


#if 0
/* --------------------------- Ghost filling - patch specific ------------------------- */

typedef void (*clawpatch_fort_copy_face_t)(const int* mx, const int* my, const int* mbc, 
										   const int* meqn,
										   double qthis[],double qneighbor[], 
										   const int* a_idir,
										   struct fclaw2d_patch_transform_data** transform_cptr);

typedef void (*clawpatch_fort_average_face_t)(const int* mx, const int* my, const int* mbc,
											  const int* meqn,
											  double qcoarse[],double qfine[],
											  double areacoarse[], double areafine[],
											  const int* idir, const int* iside,
											  const int* num_neighbors,
											  const int* refratio, const int* igrid,
											  const int* manifold, 
											  struct fclaw2d_patch_transform_data** transform_cptr);
	
typedef void (*clawpatch_fort_interpolate_face_t)(const int* mx, const int* my, const int* mbc,
												  const int* meqn,
												  double qcoarse[],double qfine[],
												  const int* idir, const int* iside,
												  const int* num_neighbors,
												  const int* refratio, const int* igrid,
												  struct fclaw2d_patch_transform_data** transform_cptr);
	
	

typedef void (*clawpatch_fort_copy_corner_t)(const int* mx, const int* my, const int* mbc,
											 const int* meqn, double this_q[],double neighbor_q[],
											 const int* a_corner,
											 struct fclaw2d_patch_transform_data** transform_cptr);

typedef void (*clawpatch_fort_average_corner_t)(const int* mx, const int* my, const int* mbc,
												const int* meqn, const int* a_refratio,
												double qcoarse[], double qfine[],
												double areacoarse[], double areafine[],
												const int* manifold,
												const int* a_corner, 
												struct fclaw2d_patch_transform_data** transform_cptr);

typedef void (*clawpatch_fort_interpolate_corner_t)(const int* mx, const int* my, const int* mbc,
													const int* meqn, const int* a_refratio, 
													double this_q[],
													double neighbor_q[], const int* a_corner,
													struct fclaw2d_patch_transform_data** transform_cptr);
	

/* --------------------------------- Regridding functions ----------------------------- */

typedef void (*clawpatch_fort_tag4refinement_t)(const int* mx,const int* my,
												const int* mbc,const int* meqn,
												const double* xlower, const double* ylower,
												const double* dx, const double* dy,
												const int* blockno,
												double q[],
												const double* tag_threshold,
												const int* init_flag,
												int* tag_patch);

typedef void (*clawpatch_fort_tag4coarsening_t)(const int* mx, const int* my,
												const int* mbc, const int* meqn,
												double xlower[], double ylower[],
												const double* dx, const double* dy,
												const int* blockno,
												double q0[],double q1[],
												double q2[],double q3[],
												const double* tag_threshold,
                                                const int* init_flag,
												int* tag_patch);


typedef int (*clawpatch_fort_exceeds_threshold_t)(int* blockno,
                                                   double qval[], 
                                                   double* qmin, double *qmax,
                                                   double quad[], 
                                                   double *dx, double *dy, 
                                                   double *xc, double *yc, 
                                                   int* tag_threshold,
                                                   int* init_flag,
                                                   int* is_ghost);
    

typedef void (*clawpatch_fort_interpolate2fine_t)(const int* mx, const int* my,
												  const int* mbc, const int* meqn,
												  double qcoarse[], double qfine[],
												  double areacoarse[], double areafine[],
												  const int* igrid, const int* manifold);
	
typedef void (*clawpatch_fort_average2coarse_t)(const int* mx, const int* my,
												const int* mbc, const int* meqn,
												double qcoarse[],double qfine[],
												double areacoarse[],double areafine[],
												const int* igrid, const int* manifold);
	

/* ----------------------------------- time stepping ---------------------------------- */

typedef void (*clawpatch_fort_timeinterp_t)(const int *mx, const int* my, const int* mbc,
											const int *meqn, const int* psize,
											double qcurr[], double qlast[],
											double qinterp[],const double* alpha,
											const int* ierror);
	
/* ------------------------------- Parallel ghost patches ----------------------------- */

typedef void (*clawpatch_fort_local_ghost_pack_t)(int *mx, int *my, int *mbc,
												  int *meqn, int *mint,
												  double qdata[], double area[],
												  double qpack[], int *psize,
												  int *packmode, int *ierror);
	
typedef void (*clawpatch_fort_local_ghost_pack_aux_t)(struct fclaw2d_global *glob,
													  struct fclaw2d_patch *patch,
													  int mint,
													  double qpack[], int extrasize,
													  int packmode, int* ierror);
	
typedef void (*clawpatch_fort_local_ghost_pack_registers_t)(struct fclaw2d_global *glob,
													  struct fclaw2d_patch *patch,
													  double qpack[], int frsize,
													  int* ierror);

/* ---------------------------------- Output functions -------------------------------- */

typedef void  (*clawpatch_fort_header_ascii_t)(char* matname1,char* matname2,
											   double* time, int* meqn, int* maux,
											   int* ngrids);

/* Write out data */
typedef void (*clawpatch_fort_output_ascii_t)(char* matname1,
											  int* mx,        int* my,
											  int* meqn,      int* mbc,
											  double* xlower, double* ylower,
											  double* dx,     double* dy,
											  double q[],
											  int* patch_num, int* level,
											  int* blockno,   int* mpirank);


/* ----------------------------- Diagnostic functions --------------------------------- */

typedef void (*clawpatch_fort_error_t)(int* blockno, int *mx, int *my, int *mbc,int *meqn,
									   double *dx, double *dy, double *xlower,
									   double *ylower, double *t, double q[],
									   double error[], double soln[]);

typedef void (*clawpatch_fort_conscheck_t)(int *mx, int *my, int* mbc, int* meqn,
										   double *dx, double *dy,
										   double area[], double q[], double sum[],
                                           double *c_kahan);

typedef double (*clawpatch_fort_area_t)(int *mx, int* my, int*mbc, double* dx,
										double* dy, double area[]);

typedef void (*clawpatch_fort_norm_t)(int* blockno, int *mx, int *my, int *mbc,int *meqn,
									  double *dx, double *dy, double area[],
									  double error[], double error_norm[]);



/* -------------------------- User convenience headers -------------------------------- */

#define TAG4REFINEMENT FCLAW_F77_FUNC(tag4refinement,TAG4REFINEMENT)
void TAG4REFINEMENT(const int* mx,const int* my,
					const int* mbc,const int* meqn,
					const double* xlower, const double* ylower,
					const double* dx, const double* dy,
					const int* blockno,
					double q[],
					const double* tag_threshold,
					const int* init_flag,
					int* tag_patch);

#define TAG4COARSENING FCLAW_F77_FUNC(tag4coarsening,TAG4COARSENING)
void TAG4COARSENING(const int* mx, const int* my,
					const int* mbc, const int* meqn,
					const double* xlower, const double* ylower,
					const double* dx, const double* dy,
					const int* blockno,
					double q0[],double q1[],
					double q2[],double q3[],
					const double* tag_threshold,
                    const int* initflag,
					int* tag_patch);


#define FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_exceeds_threshold, \
                                 FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD)

int FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD(int* blockno,
                                        double qval[], 
                                        double* qmin, double *qmax,
                                        double quad[], 
                                        double *dx, double *dy, 
                                        double *xc, double *yc, 
                                        int* tag_threshold,
                                        int* init_flag,
                                        int* is_ghost);
#endif
/* --------------------- Dimension independent defs and routines ---------------------- */

/* These headers are independent of dimension and clawpack version */
typedef void  (*clawpatch_fort_header_ascii_t)(const char* matname1,const char* matname2,
                                               const double* time, const int* meqn, 
                                               const int* maux, const int* ngrids);

/* Even though this is for 3d patches, we still assume that tagging can only be
   dependent on two dimensions */
typedef int (*clawpatch_fort_exceeds_threshold_t)(const int *blockno,
                                                  const double *qval, 
                                                  const double *qmin, 
                                                  const double *qmax,
                                                  const double quad[], 
                                                  const double *dx, 
                                                  const double *dy, 
                                                  const double *xc, 
                                                  const double *yc, 
                                                  const double *tag_threshold,
                                                  const int    *init_flag,
                                                  const int    *is_ghost);

/* ----------------------------- Fortran headers ---------------------------------------*/

#define FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_get_refinement_criteria, \
                                 FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA)
int FCLAW2D_CLAWPATCH_GET_REFINEMENT_CRITERIA();


/* ------------------------------- General threshold ---------------------------------- */
#define FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_exceeds_threshold, \
                                  FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD)

int FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD(const int *blockno,
                                        const double *qval, 
                                        const double *qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *tag_threshold,
                                        const int *init_flag,
                                        const int *is_ghost);

/* ----------------------------- Value threshold -------------------------------------- */
#define FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_value_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH)
    
int FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH(const int* blockno,
                                       const double *qval, 
                                       const double* qmin, 
                                       const double *qmax,
                                       const double quad[], 
                                       const double *dx, 
                                       const double *dy, 
                                       const double *xc, 
                                       const double *yc, 
                                       const double* tag_threshold,
                                       const int* init_flag,
                                       const int* is_ghost);
    
/* ----------------------------- difference threshold --------------------------------- */

#define FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_difference_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(const int    *blockno,
                                            const double *qval, 
                                            const double *qmin, 
                                            const double *qmax,
                                            const double quad[], 
                                            const double *dx, 
                                            const double *dy, 
                                            const double *xc, 
                                            const double *yc, 
                                            const double *tag_threshold,
                                            const int *init_flag,
                                            const int *is_ghost);

/* --------------------------------- minmax threshold --------------------------------- */

#define FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_minmax_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH(const int *blockno,
                                        const double *qval, 
                                        const double* qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *tag_threshold,                        
                                        const int *init_flag,
                                        const int *is_ghost);

/* ------------------------------- gradient threshold --------------------------------- */
#define FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH \
                  FCLAW_F77_FUNC(fclaw2d_clawpatch_gradient_exceeds_th, \
                                 FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH)

int FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH(const int *blockno,
                                          const double *qval, 
                                          const double* qmin, 
                                          const double *qmax,
                                          const double quad[], 
                                          const double *dx, 
                                          const double *dy, 
                                          const double *xc, 
                                          const double *yc, 
                                          const double *tag_threshold,
                                          const int *init_flag,
                                          const int *is_ghost);


/* ------------------------------- user threshold --------------------------------- */
#define USER_EXCEEDS_TH \
                  FCLAW_F77_FUNC(user_exceeds_th, \
                                 USER_EXCEEDS_TH)

int USER_EXCEEDS_TH(const int *blockno,
                    const double *qval, 
                    const double* qmin, 
                    const double *qmax,
                    const double quad[], 
                    const double *dx, 
                    const double *dy, 
                    const double *xc, 
                    const double *yc, 
                    const double *tag_threshold,
                    const int *init_flag,
                    const int *is_ghost);

    
#ifdef __cplusplus
}
#endif

#endif

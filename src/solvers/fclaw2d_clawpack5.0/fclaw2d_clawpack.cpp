/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#include "amr_includes.H"

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "clawpack_fort.H"
#include "fclaw2d_solvers.H"
#include "fclaw2d_clawpack.H"

void set_clawpack_parms(fclaw2d_domain_t* domain,fclaw2d_clawpack_parms_t* clawpack_parms)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    ddata->clawpack_parms = (void*) clawpack_parms;
}

static
fclaw2d_clawpack_parms_t* get_clawpack_parms(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    fclaw2d_clawpack_parms_t *clawpack_parms = (fclaw2d_clawpack_parms_t*) ddata->clawpack_parms;
    return clawpack_parms;
}

static
fclaw2d_clawpack_patch_data_t* get_clawpack_patch_data(ClawPatch *cp)
{
    /* We need the this cast here, because ClawPatch::clawpack_patch_data() only returns a
       void* object */
    fclaw2d_clawpack_patch_data_t *wp = (fclaw2d_clawpack_patch_data_t*) cp->clawpack_patch_data();
    return wp;
}


void fclaw2d_clawpack_define_auxarray(fclaw2d_domain_t* domain, ClawPatch *cp)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;

    fclaw2d_clawpack_parms_t *clawpack_parms = get_clawpack_parms(domain);
    int maux = clawpack_parms->maux;

    // Construct an index box to hold the aux array
    int ll[2], ur[2];
    ll[0] = 1-mbc;
    ll[1] = 1-mbc;
    ur[0] = mx + mbc;
    ur[1] = my + mbc;
    Box box(ll,ur);

    // get solver specific data stored on this patch
    fclaw2d_clawpack_patch_data_t *clawpack_patch_data = get_clawpack_patch_data(cp);
    clawpack_patch_data->auxarray.define(box,maux);
    clawpack_patch_data->maux = maux; // set maux in solver specific patch data (for completeness)
}

void fclaw2d_clawpack_get_auxarray(fclaw2d_domain_t* domain,
                               ClawPatch *cp, double **aux, int* maux)
{
    fclaw2d_clawpack_parms_t *clawpack_parms = get_clawpack_parms(domain);
    *maux = clawpack_parms->maux;

    fclaw2d_clawpack_patch_data_t *clawpack_patch_data = get_clawpack_patch_data(cp);
    *aux = clawpack_patch_data->auxarray.dataPtr();
}


void fclaw2d_clawpack_setprob(fclaw2d_domain_t* domain)
{
    fclaw2d_clawpack_parms_t * clawpack_parms   = get_clawpack_parms(domain);

    setprob_();

    set_amr_module_(clawpack_parms->mwaves,
                    clawpack_parms->mcapa,
                    clawpack_parms->mthlim,
                    clawpack_parms->method);
}


/* This should only be called when a new ClawPatch is created. */
void fclaw2d_clawpack_setaux(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ----------------------------------------------------------- */
    fclaw2d_clawpack_define_auxarray(domain,cp);

    /* ----------------------------------------------------------- */
    // Pointers needed to pass to class setaux call, and other setaux
    // specific arguments
    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    /* ----------------------------------------------------------- */
    // Classic setaux call
    setaux_(mbc,mx,my,xlower,ylower,dx,dy,maux,aux);
}

void fclaw2d_clawpack_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------- */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    /* ------------------------------------------------------- */
    // Call to classic Clawpack 'qinit' routine.
    qinit_(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux);
}

void fclaw2d_clawpack_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------- */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    /* ------------------------------------------------------- */
    // Classic call to b4step2(..)
    b4step2_(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux);
}

void fclaw2d_clawpack_src2(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------- */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    /* ------------------------------------------------------- */
    // Classic call to src2(..)
    src2_(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt);
}


/* Use this to return only the right hand side of the clawpack algorithm */
double fclaw2d_clawpack_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs)
{
    /* This should evaluate the right hand side, but not actually do the update.
       This will be useful in cases where we want to use something other than
       a single step method.  For example, in a RK scheme, one might want to
       call the right hand side to evaluate stages. */
    return 0;
}


void fclaw2d_clawpack_bc2(fclaw2d_domain *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      double t,
                      double dt,
                      fclaw_bool intersects_phys_bdry[],
                      fclaw_bool time_interp)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------ */
    // Set up boundary conditions
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    fclaw2d_block_data_t *bdata = get_block_data(this_block);
    int *block_mthbc = bdata->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[NumFaces];
    for(int i = 0; i < NumFaces; i++)
    {
        if (intersects_phys_bdry[i])
        {
            mthbc[i] = block_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }

    /* ------------------------------------------------------- */
    // Pointers needed to pass to Fortran


    /*
      We may be imposing boundary conditions on time-interpolated data;
      and is being done just to help with fine grid interpolation.
      In this case, this boundary condition won't be used to update
      anything
    */
    double *q = cp->q_time_sync(time_interp);

    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    /* ------------------------------------------------------- */
    // Classic call to bc2(..)
    bc2_(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc);
}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double fclaw2d_clawpack_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    const amr_options_t* gparms              = get_domain_parms(domain);
    ClawPatch *cp                            = get_clawpatch(this_patch);
    // fclaw2d_clawpack_patch_data_t *clawpack_patch_data = get_clawpack_patch_data(cp);
    fclaw2d_clawpack_parms_t * clawpack_parms   = get_clawpack_parms(domain);

    set_block_(&this_block_idx);

    int level = this_patch->level;

    double* qold = cp->q();

    double *aux;
    int maux;
    fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

    cp->save_current_step();  // Save for time interpolation

    // Global to all patches
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    // Specific to the patch
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    // Specific to solver
    int mwaves = clawpack_parms->mwaves;

    int maxm = max(mx,my);

    double cflgrid = 0;

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    // Replace this with a call to "step2" at some point...
    printf("here ...1\n");
    clawpatch2_(maxm, meqn, maux, mbc, clawpack_parms->method,
                clawpack_parms->mthlim, clawpack_parms->mcapa, mwaves,
                mx, my, qold,
                aux, dx, dy, dt, cflgrid, xlower, ylower,
                level,t, fp, fm, gp, gm);
    printf("here ...2\n");
    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}

double fclaw2d_clawpack_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt)
{
    const fclaw2d_clawpack_parms_t* clawpack_parms = get_clawpack_parms(domain);

    fclaw2d_clawpack_b4step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    double maxcfl = fclaw2d_clawpack_step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    if (clawpack_parms->src_term > 0)
    {
        fclaw2d_clawpack_src2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);
    }
    return maxcfl;
}



static
void cb_dump_auxarray(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    int dump_patchno = *((int *) user);

    int numb4 = domain->blocks[this_block_idx].num_patches_before;
    if (this_patch_idx == dump_patchno + numb4)
    {
        const amr_options_t* gparms              = get_domain_parms(domain);
        ClawPatch *cp                            = get_clawpatch(this_patch);
        // fclaw2d_clawpack_patch_data_t *clawpack_patch_data = get_clawpack_patch_data(cp);

        double *aux;
        int maux;
        fclaw2d_clawpack_get_auxarray(domain,cp,&aux,&maux);

        int mx = gparms->mx;
        int my = gparms->my;
        int mbc = gparms->mbc;

        int k = 0;
        for (int m = 0; m < maux; m++)
        {
            for(int j = 1-mbc; j <= my+mbc; j++)
            {
                for(int i = 1-mbc; i <= mx+mbc; i++)
                {
                    printf("q[%2d,%2d,%2d] = %24.16e\n",i,j,m,aux[k]);
                    k++;
                }
                printf("\n");
            }
            printf("\n");
            printf("\n");
        }
    }
}


void dump_auxarray(fclaw2d_domain_t *domain, int dump_patchno)
{
    printf("Dumping patch (time_interp) %d\n",dump_patchno);
    fclaw2d_domain_iterate_patches(domain,
                                   cb_dump_auxarray,
                                   &dump_patchno);

}


fclaw2d_clawpack_parms_t*  fclaw2d_clawpack_parms_new(sc_options_t *opt)
{
    /* This is a good example of a place where a [Section] would be nice.
       Something like [clawpack].  See fclaw_defaults.ini */

    fclaw2d_clawpack_parms_t *clawpack_parms;

    clawpack_parms = FCLAW2D_ALLOC_ZERO(fclaw2d_clawpack_parms_t, 1);

    /* Array of SpaceDim many values, with no defaults is set to all 0's */
    amr_options_add_int_array (opt, 0, "order", &clawpack_parms->order_string, NULL,
                               &clawpack_parms->order, SpaceDim,
                               "Normal and transverse orders");

    sc_options_add_int (opt, 'M', "mcapa", &clawpack_parms->mcapa, -1,
                        "Location of capacity function in aux array [-1]");

    sc_options_add_int (opt, 0, "maux", &clawpack_parms->maux, 0,
                        "Number of auxiliary variables [0]");

    sc_options_add_int (opt, 0, "src_term", &clawpack_parms->src_term, 0,
                        "Source term option [0]");

    sc_options_add_int (opt, 0, "mwaves", &clawpack_parms->mwaves, 1,
                        "Number of waves [1]");

    /* Array of mwaves many values */
    amr_options_add_int_array (opt, 0, "mthlim", &clawpack_parms->mthlim_string, NULL,
                               &clawpack_parms->mthlim, clawpack_parms->mwaves,
                               "Waves limiters (one for each wave)");
    /*
     * At this point amropt->mthlim is allocated, but for precisely 1 wave,
     * since that is the default value provided to clawpack_parms->waves above.
     * If the number of waves (clawpack_parms->mwaves) is changed by option
     * parsing, then this array needs to be converted anew. This can be done
     * by calling fclaw2d_clawpack_postprocess_parms.
     */


    /* -----------------------------------------------------------------------
       Read in options from file
       ----------------------------------------------------------------------- */

    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-W' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.

       WARNING: The long option name must be unique within the whole program.
                Just 'inifile' is already used in amr_options.c.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'W', "clawpack_inifile",
                            "Read clawpack options from this file [default : fclaw2d_clawpack.ini]");

    /* -----------------------------------------------------------------------
       This is the default file that will be read if no command line options are
       given.  This file is tracked by Git.
       ----------------------------------------------------------------------- */
   sc_options_load (sc_package_id, SC_LP_ALWAYS, opt, "fclaw2d_clawpack.ini");

   fclaw2d_clawpack_postprocess_parms(clawpack_parms);

   return clawpack_parms;

}

void fclaw2d_clawpack_checkparms(fclaw2d_clawpack_parms_t* clawpack_parms, amr_options_t* gparms)
{
    /* -----------------------------------------------------------------------
       Set up 'method' vector used by Clawpack.
       ------------------------------------------------------------------------ */
    clawpack_parms->method[0] = gparms->use_fixed_dt;

    clawpack_parms->method[1] = clawpack_parms->order[0];
    if (SpaceDim == 2)
    {
        clawpack_parms->method[2] = clawpack_parms->order[1];
    }
    else
    {
        clawpack_parms->method[2] = 10*clawpack_parms->order[1] + clawpack_parms->order[2];
    }
    clawpack_parms->method[3] = gparms->verbosity;
    clawpack_parms->method[4] = clawpack_parms->src_term;
    clawpack_parms->method[5] = clawpack_parms->mcapa;
    clawpack_parms->method[6] = clawpack_parms->maux;

    /* Should also check mthbc, mthlim, etc. */

}

void fclaw2d_clawpack_postprocess_parms(fclaw2d_clawpack_parms_t* clawpack_parms)
{
    /* -----------------------------------------------------------------------
       Some post-processing of arrays
       ------------------------------------------------------------------------ */

    amr_options_convert_int_array (clawpack_parms->mthlim_string, &clawpack_parms->mthlim,
                                   clawpack_parms->mwaves);

    amr_options_convert_int_array (clawpack_parms->order_string, &clawpack_parms->order,
                                   SpaceDim);
}


void fclaw2d_clawpack_parms_delete(fclaw2d_clawpack_parms_t* clawpack_parms)
{
    SC_FREE(clawpack_parms->order);
    SC_FREE(clawpack_parms->mthlim);
    FCLAW2D_FREE(clawpack_parms);
    clawpack_parms = NULL;
}

static
void fclaw2d_clawpack_patch_data_new(void** wp)
{
    fclaw2d_clawpack_patch_data_t* clawpack_patch_data;
    clawpack_patch_data = new fclaw2d_clawpack_patch_data_t;
    *wp = (void*) clawpack_patch_data;

    // or?
    // *wp = (void*) FCLAW2D_ALLOC_ZERO (fclaw2d_clawpack_patch_data_t, 1);
}

static
void fclaw2d_clawpack_patch_data_delete(void **wp)
{
    fclaw2d_clawpack_patch_data_t *clawpack_patch_data = (fclaw2d_clawpack_patch_data_t*) *wp;
    delete clawpack_patch_data;

    // or?
    // FCLAW2D_FREE(clawpack_patch_data);

    *wp = (void*) NULL;
}


void fclaw2d_clawpack_link_to_clawpatch()
{
    /* These are called whenever a new ClawPatch is created. */
    ClawPatch::f_clawpack_patch_data_new = &fclaw2d_clawpack_patch_data_new;
    ClawPatch::f_clawpack_patch_data_delete = &fclaw2d_clawpack_patch_data_delete;
}

void  fclaw2d_clawpack_link_solvers(fclaw2d_domain_t* domain)
{
    const fclaw2d_clawpack_parms_t* clawpack_parms = get_clawpack_parms(domain);

    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    if (clawpack_parms->maux > 0)
    {
        sf->f_patch_setup          = &fclaw2d_clawpack_setaux;
    }
    else
    {
        sf->f_patch_setup          = &amr_dummy_patch_setup;
    }

    sf->f_patch_initialize         = &fclaw2d_clawpack_qinit;
    sf->f_patch_physical_bc        = &fclaw2d_clawpack_bc2;

    /* Calls b4step2, step2 and src2 */
    sf->f_patch_single_step_update = &fclaw2d_clawpack_update;

    /* This is needed so that a ClawPatch knows how to create data for a clawpack solver */
    fclaw2d_clawpack_link_to_clawpatch();
}

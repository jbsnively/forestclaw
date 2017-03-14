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

#include "swirl_user.h"

#include <fclaw2d_forestclaw.h>



static void *
options_register_user (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t* user = (user_options_t*) package;

    /* [user] User options */
    sc_options_add_double (opt, 0, "period", &user->period, 4,
                           "Period of the flow field [4]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                           "Clawpack_version (4 or 5) [5]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_check_user (fclaw_app_t * app, void *package, void *registered)
{
    return FCLAW_NOEXIT;
}

static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register_user,
    NULL,
    options_check_user,
    NULL
};


static
void register_user_options (fclaw_app_t * app,
                            const char *configfile,
                            user_options_t* user)
{
    FCLAW_ASSERT (app != NULL);

    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);
}

const user_options_t* swirl_user_get_options(fclaw2d_global_t* glob)
{
    return swirl_user_get_options_old(glob->domain);
}

const user_options_t* swirl_user_get_options_old(fclaw2d_domain_t* domain)
{
    fclaw_app_t* app;
    // app = fclaw2d_global_get_app(glob);
    app = fclaw2d_domain_get_app(domain);

    const user_options_t* user = (user_options_t*) fclaw_app_get_user(app);

    return (user_options_t*) user;
}

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, amr_options_t* gparms)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL;

    /* Map unit square to disk using mapc2m_disk.f */

    gparms->manifold = 0;
    conn = p4est_connectivity_new_unitsquare();
    cont = fclaw2d_map_new_nomap();

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
    return domain;
}

static
void run_program(fclaw2d_global_t* glob, fclaw_app_t* app)
{
    user_options_t           *user;     
    fclaw2d_clawpatch_options_t *clawpatchopt;


    clawpatchopt = fclaw2d_clawpatch_get_options(glob);
    printf("mx = %d\n", clawpatchopt->mx);
    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);
    fclaw2d_domain_set_app (glob->domain,app);

    user = (user_options_t*) fclaw_app_get_user(app);
    if (user->claw_version == 4)
    {
      fc2d_clawpack46_set_vtable_defaults();
    }
    else if (user->claw_version == 5)
    {
      fc2d_clawpack5_set_vtable_defaults();
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */

    fclaw2d_initialize(glob);
    fclaw2d_run(glob);
    
    fclaw2d_map_context_t* cont = fclaw2d_domain_get_map_context(glob->domain);
    fclaw2d_map_destroy(cont);
    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;
    user_options_t              suser, *user = &suser;
    amr_options_t               *gparms;

    fclaw2d_global_t         *glob;
    fclaw2d_domain_t         *domain;
    sc_MPI_Comm mpicomm;

    int retval;

    glob = fclaw2d_global_new();

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, user);
    fclaw_forestclaw_register(app,"fclaw_options.ini");  /* Register gparms */

    /* All libraries that might be needed should be registered here */
    fclaw2d_clawpatch_register(glob,app,"fclaw_options.ini");
    fc2d_clawpack46_register(app,"fclaw_options.ini");    /* [clawpack46] */
    fc2d_clawpack5_register (app,"fclaw_options.ini");     /* [clawpack5] */
    register_user_options   (app,"fclaw_options.ini",user);  /* [user] */

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    fclaw2d_clawpatch_link_app(app);
    
    /* at this point gparms is valid */
    gparms = fclaw_forestclaw_get_options(app);
    mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
    domain = create_domain(mpicomm, gparms);
    
    fclaw2d_global_set_domain(glob, domain);
    fclaw2d_global_set_gparms (glob, gparms);

    /* Run the program */
    if (!retval & !vexit)
    {
        run_program(glob, app);
    }
    
    fclaw2d_global_destroy(glob);
    fclaw_forestclaw_destroy(app);
    fclaw_app_destroy (app);

    return 0;
}

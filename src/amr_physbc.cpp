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

#include "amr_forestclaw.H"
#include "amr_solver_typedefs.H"

/* -----------------------------------------------------------------------------
   Physical boundary conditions
   ----------------------------------------------------------------------------- */

static
void cb_set_phys_bc(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    // int numfaces = get_faces_per_patch(domain);
    fclaw_bool intersects_bc[NumFaces];
    double curr_time = *((double*) user);
    double dt = 1e20;
    get_phys_boundary(domain,this_block_idx,this_patch_idx,intersects_bc);

    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    (ddata->f_patch_physbc_ptr)(domain,
                                this_patch,
                                this_block_idx,
                                this_patch_idx,
                                curr_time,dt,
                                intersects_bc);
}


void set_phys_bc(fclaw2d_domain_t *domain, int a_level, double a_level_time)
{
    fclaw2d_domain_iterate_level(domain, a_level,cb_set_phys_bc,(void *) &a_level_time);
}

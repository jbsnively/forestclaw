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

#include <fclaw2d_options.h>
#include <fclaw2d_global.h>
#include <fclaw_base.h>
#include <fclaw_package.h>

static int s_fclaw2d_options_package_id = -1;

/* ---------------------------------------------------------
   Public interface to ForestClaw options
   --------------------------------------------------------- */

void fclaw2d_options_store (fclaw2d_global_t *glob, fclaw_options_t* gparms)
{
    int id;

    //FCLAW_ASSERT(s_fclaw2d_options_package_id == -1);
    id = fclaw_package_container_add_pkg(glob,gparms);
    s_fclaw2d_options_package_id = id;
}

fclaw_options_t* fclaw2d_get_options(fclaw2d_global_t* glob)
{
    fclaw_options_t *gp = (fclaw_options_t*) 
            fclaw_package_get_options(glob, s_fclaw2d_options_package_id);
    return gp;
}


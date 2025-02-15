/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "slosh_user.h"

#include <fclaw_pointer_map.h>

static void *
slosh_register (slosh_user_options_t *user, sc_options_t * opt)
{
    user->is_registered = 1;
    return NULL;
}

static void
slosh_destroy(slosh_user_options_t *user)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    slosh_user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (slosh_user_options_t*) package;

    return slosh_register(user,opt);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    slosh_user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (slosh_user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    slosh_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    NULL,
    options_destroy
};

slosh_user_options_t* slosh_options_register (fclaw_app_t * app,
                                        const char *section,
                                        const char *configfile)
{
    slosh_user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (slosh_user_options_t, 1);
    fclaw_app_options_register (app, section, configfile, &options_vtable_user,
                                user);

    return user;
}

void slosh_options_store (fclaw2d_global_t* glob, slosh_user_options_t* user)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"slosh-user") == NULL);
    fclaw_pointer_map_insert(glob->options, "slosh-user", user, NULL);
}

slosh_user_options_t* slosh_get_options(fclaw2d_global_t* glob)
{
    slosh_user_options_t* user = (slosh_user_options_t*) 
                              fclaw_pointer_map_get(glob->options, "slosh-user");
    FCLAW_ASSERT(user != NULL);
    return user;   
}


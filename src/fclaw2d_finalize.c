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

#include <sc_statistics.h>
#include <fclaw2d_global.h>

#include <fclaw2d_diagnostics.h>
#include <fclaw2d_options.h>
#include <fclaw2d_map.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_forestclaw.h>

/**
 *  This macro defines two utility functions for reading and writing CSV data of a specified type,
 *  typename, to a file. The generated functions are named read_input_and_check_<extension>
 *  and write_output_<extension>, where <extension> is provided as a macro argument.
 *
 *  Usage:
 *
 *  DEFINE_READ_AND_WRITE(int, int32, d) will generate the following functions:
 *  static int read_input_and_check_int32(FILE* file, const char* name, int actual)
 *  static void write_output_int32(FILE* file, const char* name, int actual)
 *
 *  Arguments:
 *
 *  typename: The data type to read/write (e.g., int, float, etc.)
 *  extension: A unique identifier for the generated functions (e.g., int32, float32, etc.)
 *  format: The format specifier used by fscanf and fprintf for the given typename (e.g., d for int, f for float, etc.)
 *
 *  Functions generated by the macro:
 *
 *  read_input_and_check_<extension>():
 *  - Reads CSV data of type 'typename' from the provided file using the specified format.
 *  - Compares the read value with the provided 'actual' value.
 *  - If reading fails prints an error message, returns 1 (failure).
 *  - If the read value does not match the actual value, prints an error message and returns 1 (failure).
 *  - If the read is successful, returns 0 (success).
 *
 *  write_output_<extension>():
 *  - Writes CSV data of type 'typename' to the provided file using the specified format.
 *  - If writing fails, prints an error message and exits the program with a non-zero status.
 */
#define DEFINE_READ_AND_WRITE(typename, extension, format) \
static int read_input_and_check_##extension(FILE* file, const char* name, typename actual) \
{ \
    char actual_name[50]; \
    typename expected; \
    int fprintf_value = fscanf(file, "%49[^,],%"#format"\n", actual_name, &expected); \
 \
    if(fprintf_value != 2) \
    { \
        fclaw_global_essentialf("Failed to read form input file\n"); \
        return 1; \
    } \
 \
    if(strcmp(name, actual_name) == 0) \
    { \
        if(expected != actual) \
        { \
            fclaw_global_essentialf("Expected %s to be %"#format", but was %"#format"\n", name, expected, actual); \
            return 1; \
        } \
    } \
    else \
    { \
        fclaw_global_essentialf("Expected to be reading %s, but got %s\n", name, actual_name); \
        return 1; \
    } \
 \
    return 0; \
} \
static void write_output_##extension(FILE* file, const char* name, typename actual) \
{ \
 \
    int fprintf_value = fprintf(file, "%s,%"#format"\n",name, actual); \
 \
    if(fprintf_value < 0) \
    { \
        fclaw_global_essentialf("Failed to write to output file\n"); \
        exit(1); \
    } \
}

DEFINE_READ_AND_WRITE(int,int,d)
DEFINE_READ_AND_WRITE(long int,long_int,ld)

/**
 * @brief Outputs the expected values for regression tests.
 *
 * This function is used to output the expected values for regression tests. If there are any
 * discrepancies, it will create a new file with the actual values and terminate the program.
 *
 * @param[in] glob Pointer to the fclaw2d_global_t object containing the global state.
 * @param[in] filename Name of the file containing the expected values for the regression tests.
 */
static void 
output_expected_values(fclaw2d_global_t* glob, const char* filename)
{
    if(glob->mpirank == 0)
    {
        // open file to get expected results
        FILE* file = fopen(filename, "r");
        if (file == NULL)
        {
            fclaw_global_essentialf("Failed to open regressions file %s\n", filename);
            exit(FCLAW_EXIT_ERROR);
        }

        int num_failures = 0;
        num_failures += read_input_and_check_int(file, "count_amr_advance", glob->count_amr_advance);
        num_failures += read_input_and_check_long_int(file, "global_num_patches", glob->domain->global_num_patches);
        num_failures += read_input_and_check_int(file, "count_amr_new_domain", glob->count_amr_new_domain);

        // close file
        fclose(file);

        if(num_failures > 0){
            // append .actual to the filename
            char* actual_filename = FCLAW_ALLOC(char, strlen(filename) + 8);
            strcpy(actual_filename, filename);
            strcat(actual_filename, ".actual");

            // open the file
            file = fopen(actual_filename, "w");
            if (file == NULL)
            {
                fclaw_global_essentialf("Failed to open actual values file %s\n", actual_filename);
                exit(FCLAW_EXIT_ERROR);
            }

            //write output file
            write_output_int(file, "count_amr_advance", glob->count_amr_advance);
            write_output_long_int(file, "global_num_patches", glob->domain->global_num_patches);
            write_output_int(file, "count_amr_new_domain", glob->count_amr_new_domain);

            //close file
            fclose(file);

            FCLAW_FREE(actual_filename);
            exit(1);
        }
    }
}
/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void fclaw2d_finalize(fclaw2d_global_t* glob)
{
    const fclaw_options_t *gparms = fclaw2d_get_options(glob);

    fclaw_global_essentialf("Finalizing run\n");
    fclaw2d_diagnostics_finalize(glob);
    if (glob->cont != NULL) {
        fclaw2d_map_destroy(glob->cont);
    }
    fclaw2d_domain_barrier (glob->domain);

    if (gparms->report_timing)
    {
        if (gparms->outstyle > 0)
        {
            /* Only call this if we have taken time steps.  For time-independent problems, we
               probably need a different report (no "amr_advance_steps") */
            fclaw2d_timer_report(glob);
        }
        else
        {
            fclaw_global_essentialf("Timing reports not generated for outstyle=0\n");
        }
    }
    if (gparms->regression_check)
    {
        fclaw_global_essentialf("Checking regression values.\n");
        output_expected_values(glob, gparms->regression_check);
        fclaw_global_essentialf("Regression check passed.\n");
    }
    fclaw2d_domain_reset(glob);
}

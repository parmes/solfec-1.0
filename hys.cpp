/* HYBRID_SOLVER implementation --> PARMEC C++ to C glue code */

/*
The MIT License (MIT)

Copyright (c) 2016 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contributors: Tomasz Koziara */

#include <stdio.h>
#include <stdlib.h>
#include "parmec8.h"

extern "C" { /* C */

/* initialize parmec */
void parmec_init (char *path, char **argv, int argc)
{
  parmec::init();

  parmec::input(path, argv, argc);
}

/* allocate force and torque exchange buffers; must be
 * called prior to parmec_{set,get}_force_and_forque routines */
void parmec_force_and_torque_alloc (int num)
{
  parmec::prescribed_body_force *p = (parmec::prescribed_body_force*)
          MAP_Find (parmec::prescribed_body_forces, (void*)(long)num, NULL);

  if (!p)
  {
    p = new parmec::prescribed_body_force;
    p->particle = num;
    p->outer_force[0] = 0.0;
    p->outer_force[1] = 0.0;
    p->outer_force[2] = 0.0;
    p->outer_torque[0] = 0.0;
    p->outer_torque[1] = 0.0;
    p->outer_torque[2] = 0.0;
    p->inner_force[0] = 0.0;
    p->inner_force[1] = 0.0;
    p->inner_force[2] = 0.0;
    p->inner_torque[0] = 0.0;
    p->inner_torque[1] = 0.0;
    p->inner_torque[2] = 0.0;
    MAP_Insert (NULL, &parmec::prescribed_body_forces, (void*)(long)num, p, NULL);
  }
}

/* prescribe body force and torque --> forces set this way are kept constant
 *                                  and used during all following time steps
 */
void parmec_set_force_and_torque (int num, double *force, double *torque)
{
  parmec::prescribed_body_force *p = (parmec::prescribed_body_force*)
          MAP_Find (parmec::prescribed_body_forces, (void*)(long)num, NULL);

  if (!p)
  {
    fprintf (stderr, "ERROR: Solfec-Parmec boundary force with Paremc id = %d has not been found", num);
    exit (1);
  }

  p->outer_force[0] = force[0]; /* prescribe outer forces */
  p->outer_force[1] = force[1];
  p->outer_force[2] = force[2];
  p->outer_torque[0] = torque[0];
  p->outer_torque[1] = torque[1];
  p->outer_torque[2] = torque[2];
}

/* read body force and torque --> these are inner parmec forces
 *                                accumulated at the latest step
 */
void parmec_get_force_and_torque (int num, double *force, double *torque)
{
  parmec::prescribed_body_force *p = (parmec::prescribed_body_force*)
          MAP_Find (parmec::prescribed_body_forces, (void*)(long)num, NULL);

  if (!p)
  {
    fprintf (stderr, "ERROR: Solfec-Parmec boundary force with Paremc id = %d has not been found", num);
    exit (1);
  }

  force[0] = p->inner_force[0];
  force[1] = p->inner_force[1];
  force[2] = p->inner_force[2];
  torque[0] = p->inner_torque[0];
  torque[1] = p->inner_torque[1];
  torque[2] = p->inner_torque[2];
}

/* perform single time integration step */
void parmec_one_step (double step, double* interval, void** interval_func, int* interval_tms, char* prefix)
{
  parmec::dem (step, step, interval, interval_func, interval_tms, prefix, 0, 0);
}

/* read angular and linear velocities of a particle */
void parmec_get_angular_and_linear (int num, double *angular, double *linear)
{
  angular[0] = parmec::angular[0][num];
  angular[1] = parmec::angular[1][num];
  angular[2] = parmec::angular[2][num];

  linear[0] = parmec::linear[0][num];
  linear[1] = parmec::linear[1][num];
  linear[2] = parmec::linear[2][num];
}

/* write angular and linear velocities of a particle */
void parmec_set_angular_and_linear (int num, double *angular, double *linear)
{
  parmec::angular[0][num] = angular[0];
  parmec::angular[1][num] = angular[1];
  parmec::angular[2][num] = angular[2];

  parmec::linear[0][num] = linear[0];
  parmec::linear[1][num] = linear[1];
  parmec::linear[2][num] = linear[2];
}

/* read rotation and position of a particle */
void parmec_get_rotation_and_position (int num, double *rotation, double *position)
{
  rotation[0] = parmec::rotation[0][num];
  rotation[1] = parmec::rotation[1][num];
  rotation[2] = parmec::rotation[2][num];
  rotation[3] = parmec::rotation[3][num];
  rotation[4] = parmec::rotation[4][num];
  rotation[5] = parmec::rotation[5][num];
  rotation[6] = parmec::rotation[6][num];
  rotation[7] = parmec::rotation[7][num];
  rotation[8] = parmec::rotation[8][num];

  position[0] = parmec::position[0][num];
  position[1] = parmec::position[1][num];
  position[2] = parmec::position[2][num];
}

/* write rotation and position of a particle */
void parmec_set_rotation_and_position (int num, double *rotation, double *position)
{
  parmec::rotation[0][num] = rotation[0];
  parmec::rotation[1][num] = rotation[1];
  parmec::rotation[2][num] = rotation[2];
  parmec::rotation[3][num] = rotation[3];
  parmec::rotation[4][num] = rotation[4];
  parmec::rotation[5][num] = rotation[5];
  parmec::rotation[6][num] = rotation[6];
  parmec::rotation[7][num] = rotation[7];
  parmec::rotation[8][num] = rotation[8];

  parmec::position[0][num] = position[0];
  parmec::position[1][num] = position[1];
  parmec::position[2][num] = position[2];
}

/* disable particle dynamics */
void parmec_disable_dynamics (int num)
{
  parmec::flags[num] |= parmec::SKIP;
}

/* get number of defined time series objects */
int parmec_get_tmsnum ()
{
  return parmec::tmsnum;
}

} /* extern C */

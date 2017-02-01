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

/* this file lists the C functions exported in hys.cpp facilitaing C-level access to PARMEC's C++ functions and global data */

#ifndef __hyshpp__
#define __hyshpp__

/* initialize parmec */
void parmec_init(char *path);

/* prescribe body force and torque --> forces set this way are kept constant
 *                                  and used during all following time steps
 */
void parmec_set_force_and_torque (int num, double *force, double *torque);

/* read body force and torque --> these are inner parmec forces
 *                                accumulated over a number of steps
 */
void parmec_get_force_and_torque (int num, int nstep, double *force, double *torque);

/* perform single time integration step */
void parmec_one_step (double step, double interval[2], char *prefix);

/* read angular and linear velocities of a body */
void parmec_get_angular_and_linear (int num, double *angular, double *linear);

/* read rotation and position of a body */
void parmec_get_rotation_and_position (int num, double *rotation, double *position);

#endif

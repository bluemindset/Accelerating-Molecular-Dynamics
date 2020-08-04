/**
 * @Author: B159973
 * @Date:	10/4/2019
 * @Course: Performance Programming - 2020
 * @University of Edinburgh
 * 
 * Providing only the most optimal solution
*/
#include <stdio.h>
#include <math.h>
#include "coord.h"
#include "omp.h"
#include "immintrin.h" // for AVX
#include <xmmintrin.h> /* Streaming SIMD Extensions Intrinsics include file */

void vis_force(int N, double *f, double *vis, double *vel);
void add_norm(int N, double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N, double *f, double *vis, double vel);

void evolve(int count, double dt)
{
  int step;
  int i, j, k, l;
  int collided;
  double Size;

  for (step = 1; step <= count; step++)
  {
    printf("timestep %d\n", step);
    printf("collisions %d\n", collisions);

/***************************************************************************//**
 * First Part - As described in report page 6
 ******************************************************************************/
    for (k = 0; k < Nbody; k++)
    {
      double pos0 = *(*(pos + 0) + k);
      double pos1 = *(*(pos + 1) + k);
      double pos2 = *(*(pos + 2) + k);
      double m = mass[k];
      double r = pos0 * pos0;
      r += pos1 * pos1;
      r += pos2 * pos2;
      r = sqrt(r);

      f[0][k] = -(vis[k] * velo[0][k]) - (vis[k] * wind[0]) - force(G * m * M_central, pos0, r);
      f[1][k] = -(vis[k] * velo[1][k]) - (vis[k] * wind[1]) - force(G * m * M_central, pos1, r);
      f[2][k] = -(vis[k] * velo[2][k]) - (vis[k] * wind[2]) - force(G * m * M_central, pos2, r);
    }
    k = 0;
/***************************************************************************//**
 * Second Part - As described in report page 8
 ******************************************************************************/
    for (i = 0; i < Nbody; i++)
    {
      double f0_i = f[0][i];
      double f1_i = f[1][i];
      double f2_i = f[2][i];
      double pos0_i = pos[0][i];
      double pos1_i = pos[1][i];
      double pos2_i = pos[2][i];
      double m, size, delta_r =0 ;
      int re;
#pragma omp simd reduction(+ \
                           : f0_i, f1_i, f2_i,k,delta_r) private(m, size, re) aligned(pos : 32,f:32)

#pragma prefetch pos,f    /* No use */
      for (j = i + 1; j < Nbody; j++)
      {

        size = radius[i] + radius[j];
        m = mass[i] * mass[j] * G;

        double delta_pos0 = pos0_i - pos[0][j];
        double delta_pos1 = pos1_i - pos[1][j];
        double delta_pos2 = pos2_i - pos[2][j];
        delta_r = delta_pos0 * delta_pos0;
        delta_r += delta_pos1 * delta_pos1;
        delta_r += delta_pos2 * delta_pos2;
        delta_r = pow(delta_r,1.5);

        union {
          double _f;
          unsigned long long _x;
        } _u;
        double r = delta_r - size;  /*size can be replaced with 1*/
        _u._f = r;
        re = (int)(_u._x >> 63);
        re = 1 - 2 * re;

        double s = m/delta_r; 
        double a =  delta_pos0 * s;
        double b = delta_pos1  * s;
        double c = delta_pos2  * s;
        /*1*/ f0_i += (-re) * a;
        /*2*/ f1_i += (-re) * b;
        /*3*/ f2_i += (-re) * c;
        /*1*/ f[0][j] += (+re) * a;
        /*2*/ f[1][j] += (+re) * b;
        /*3*/ f[2][j] += (+re) * c;

        collisions -=  (re>>15);
         
        k += 1;
      }
      pos[0][i] += dt * velo[0][i];
      velo[0][i] += dt * (f0_i / mass[i]);
      pos[1][i] += dt * velo[1][i];
      velo[1][i] += dt * (f1_i / mass[i]);
      pos[2][i] += dt * velo[2][i];
      velo[2][i] += dt * (f2_i / mass[i]);
    }
  }
}
/**
 * @Author: B159973
 * @Date:	10/4/2019
 * @Course: Performance Programming - 2020
 * @University of Edinburgh
*/
/*
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
 *
 *  Nbody	  Number of particles
 *  Npair	  Number of particle pairs
 *  pos		  Position of the particles
 *  r         distance of partice from central mass 
 *  vel		  velocity of the particles
 *  f		  Forces acting on each particle
 *  vis       viscosity coefficient for each particle
 *  mass	  mass of each particle
 *  delta_pos separation vector for each particle pair
 *  delta_r	  separation for each particle pair
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };
    
DEF double *__attribute__((aligned(32))) pos[Ndim];
DEF double *__attribute__((aligned(32))) velo[Ndim];
DEF double *__attribute__((aligned(32)))f[Ndim], *vis, *mass, *radius;
DEF double *__attribute__((aligned(32)))delta_pos[3];
DEF double  __attribute__((aligned(32))) wind[Ndim];
DEF double *r;
DEF double *delta_r;
DEF int collisions;


#define G 2.0
#define M_central 1000.0

void evolve(int Nstep, double dt);
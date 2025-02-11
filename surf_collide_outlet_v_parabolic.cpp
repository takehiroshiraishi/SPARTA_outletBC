/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "surf_collide_outlet_v_parabolic.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_extra.h"
#include "error.h"

#include <iostream>

using namespace SPARTA_NS;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideOutletVParabolic::SurfCollideOutletVParabolic(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 5 ) error->all(FLERR,"Illegal surf_collide outletv command");

  // macroscopic velocity
  meanv = atof(arg[2]);
  directionv = atoi(arg[3]);
  maxx = atof(arg[4]);

  // initialize RNG
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVParabolic::~SurfCollideOutletVParabolic()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   general outlet boundary condition with prescribed velocity
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVParabolic::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  double tangent1[3];
  double *v = ip->v; // particle velocity
  double *x = ip->x; // particle position
  double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  double u_p = -dot;
  double v_p = MathExtra::len3(tangent1);

  MathExtra::norm3(tangent1); // normalize vector
  double u0dotn = 4 * meanv / (maxx * maxx) * x[directionv] * (maxx - x[directionv]); 
  double u0dottan = 0;

  v[0] = -(2*u0dotn - u_p)*norm[0] + (2*u0dottan - v_p)*tangent1[0];
  v[1] = -(2*u0dotn - u_p)*norm[1] + (2*u0dottan - v_p)*tangent1[1];
  v[2] = -(2*u0dotn - u_p)*norm[2] + (2*u0dottan - v_p)*tangent1[2];

  if((2*u0dotn - u_p) >= 0){
    //remove particle
    ip = NULL;
  }else if( random->uniform() >= abs((2*u0dotn - u_p)/u_p)){
    //remove particle
    ip = NULL;
  }
  
  return NULL;
}
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
#include "surf_collide_outlet_vgrad_yfree.h"
#include "surf.h"
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

SurfCollideOutletVgradYfree::SurfCollideOutletVgradYfree(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 8 ) error->all(FLERR,"Illegal surf_collide outlet_vgrayfree command");

  // macroscopic velocity
  directionv = atoi(arg[2]);
  directiongrad = atoi(arg[3]);
  grad = atof(arg[4]);
  refvel = atof(arg[5]);
  width = atof(arg[6]);
  maxwidth = atof(arg[7]);
  margin = (maxwidth - width) * 0.5;

  // initialize RNG
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradYfree::~SurfCollideOutletVgradYfree()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVgradYfree::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  double tangent1[3];
  double u0[3];

  double *v = ip->v; // particle velocity
  double *x = ip->x; // particle position
  double dot = MathExtra::dot3(v,norm);

  std::fill(u0, u0+3, 0.0);
  if (x[directiongrad] < margin) {
    u0[directionv] = refvel;
  } else if (x[directiongrad] > maxwidth - margin){
    u0[directionv] = width * grad + refvel;
  } else {
    u0[directionv] = (x[directiongrad]-margin) * grad + refvel;
  }

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  double u_p = -dot;
  double u0dotn = -(u0[0]*norm[0] + u0[1]*norm[1] + u0[2]*norm[2]);

  v[0] = -(2*u0dotn - u_p)*norm[0] + tangent1[0];
  v[1] = -(2*u0dotn - u_p)*norm[1] + tangent1[1];
  v[2] = -(2*u0dotn - u_p)*norm[2] + tangent1[2];

  if((2*u0dotn - u_p) >= 0){
    //remove particle
    ip = NULL;
  }else if( random->uniform() >= abs((2*u0dotn - u_p)/u_p)){
    //remove particle
    ip = NULL;
  }
  
  return NULL;
}
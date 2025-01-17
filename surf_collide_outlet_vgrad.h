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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(outletvgrad,SurfCollideOutletVgrad)

#else

#ifndef SPARTA_SURF_COLLIDE_OUTLET_VGRAD_H
#define SPARTA_SURF_COLLIDE_OUTLET_VGRAD_H

#include "surf_collide.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfCollideOutletVgrad : public SurfCollide {
 public:
  SurfCollideOutletVgrad(class SPARTA *, int, char **);
  SurfCollideOutletVgrad(class SPARTA *sparta) : SurfCollide(sparta) {} // needed for Kokkos
  virtual ~SurfCollideOutletVgrad();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);

 protected:
  int directionv;         // velocity direction x:0 y:1 z:2
  int directiongrad;      // direction in which velocity has a gradient x:0 y:1 z:2
  double grad;            // gradient
  double refvel;          // reference velocity at the zeropoint
  double width;           // width of region where gradient exists
  double maxwidth;        // width of the system
  double margin;          // width of the region without gradient
  class RanKnuth *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

/**
 * @file die-swell_ViscoElastic.c
 * @brief This file contains the simulation code for the die-swell of a viscoelastic liquid being extruded out of a die. 
 * @author Vatsal Sanjay
 * @version 2.0
 * These are viscoelastic simulations!
 * @date Oct 22, 2024
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "log-conform-viscoelastic-scalar-2D.h"

#define FILTERED // Smear density and viscosity jumps
#include "two-phaseVE.h"

#define logFile "log-die-swell_ViscoElastic.dat"

#include "navier-stokes/conserving.h"
#include "tension.h"

#define tsnap (1e-2)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define AErr (1e-3)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)

#define epsilon (1e-1)
#define R2(x,y,z) (sq(y) + sq(z))


int MAXlevel;
// We -> Weber number
// Re_s -> Solvent Reynolds number
// Re_a -> air Reynolds number
// Wi -> Weissberger number
// El -> elasticity number (this is \tilde{G}*sq(Wi)), El = G\lambda^2/(\rho R^2).
// muR -> ratio of air viscosity to solvent viscosity

double We, Re_s, Re_a, muR, Wi, El, tmax;
char nameOut[80], dumpFile[80];

// boundary conditions
u.n[top] = neumann(0.0);
p[top] = dirichlet(0.0);
u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);

// f[left] = dirichlet(y > 1. + epsilon ? 0.0 : y < 1. - epsilon ? 1.0 : 0.5 * (1.0 + tanh((1e0 - R2(x,y,z)) / epsilon)));
f[left] = neumann(0.);
u.n[left] = dirichlet(clamp(f[],0.,1.)*2e0*(1-R2(x,y,z)));
u.t[left] = dirichlet(0.0);
p[left] = dirichlet(0.0);

int main(int argc, char const *argv[]) {


  // Values taken from the terminal
  L0 = 1e1; //atof(argv[1]);
  MAXlevel = 11; //atoi(argv[2]);

  We = 1e2; //atof(argv[3]);
  Re_s = 1e0; //atof(argv[4]);
  muR = 1e-2;

  tmax = 0.1; //atof(argv[5]);

  Wi = 1e0; //atof(argv[6]);
  El = 1e0; //atof(argv[7]);

  TOLelastic = (1.0-1e-3);

  A12[left] = dirichlet(clamp(f[],0.,1.)*El*(-2e0*y)/Wi);
  A11[left] = dirichlet(clamp(f[],0.,1.)*(1 + 2*El*sq(-2e0*y)) + (1.0-clamp(f[],0.,1.)));

  init_grid (1 << 6);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");


  rho1 = 1., rho2 = 1e-3;
  mu1 = 1e0/Re_s, mu2 = muR/Re_s;
  lambda1 = Wi, lambda2 = 0.;
  G1 = El/sq(Wi), G2 = 0.;
  f.sigma = 1.0/We;

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine(x < 2*epsilon && R2(x,y,z) < sq(1+2*epsilon) && level < MAXlevel);
    // fraction (f, 1-R2(x,y,z)-x/epsilon);
    fraction (f, intersection(1-R2(x,y,z), epsilon-x));
    foreach(){
      u.x[] = clamp(f[],0.,1.)*2e0*(1-R2(x,y,z));
      A12[] = clamp(f[],0.,1.)*El*(-2e0*y)/Wi;
      A11[] = clamp(f[],0.,1.)*(1 + 2*El*sq(-2e0*y))+(1-clamp(f[],0.,1.));
      u.y[] = 0.0;
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA[];
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, A22, A11, A12, AThTh, KAPPA},
      (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr},
      MAXlevel, 4);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Writing Log
*/
event logWriting (i++) {

  scalar ke_field[];
  foreach(){ 
    ke_field[] = (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  double ke = statsf(ke_field).sum;

  scalar pos[];
  position (f, pos, {0,1,0});
  double ymax = statsf(pos).max;
  scalar posX[];
  position (f, posX, {1,0,0});
  double xmax = statsf(posX).max;

  static FILE * fp = NULL;
  if (!fp) {
    char name[80];
    sprintf(name, logFile);
    fp = fopen(name, t == 0 ? "w" : "a");
    if (fp == NULL) {
      fprintf(stderr, "Error opening log file\n");
      return 1;
    }
    if (t == 0) {
      fprintf(ferr, "Level %d, We %2.1e, Re_s %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re_s, muR, Wi, El);
      fprintf(ferr, "t dt ke xmax ymax\n");
      fprintf(fp, "Level %d, We %2.1e, Re_s %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re_s, muR, Wi, El);
      fprintf(fp, "t dt ke xmax ymax\n");
    }
  }

  fprintf(fp, "%g %g %g %g %g\n", t, dt, ke, xmax, ymax);
  fprintf(ferr, "%g %g %g %g %g\n", t, dt, ke, xmax, ymax);
  fflush(fp);

  if(ke < -1e-10) return 1;
  if (xmax > 0.9*L0){ 
    fprintf(stderr, "Liquid is close to the right boundary. Stopping simulation!\n");      
    return 1;
  }
}

/**
## Ending Simulation
*/
event end (t = end) {
  fprintf(stderr, "Level %d, We %2.1e, Re_s %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re_s, muR, Wi, El);
}

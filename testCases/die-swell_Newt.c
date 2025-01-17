/**
 * @file die-swell_Newt.c
 * @brief This file contains the simulation code for the die-swell of a viscoelastic liquid being extruded out of a die. 
 * @author Vatsal Sanjay
 * @version 1.0
 * This is the Newtonian limit!
 * @date Oct 22, 2024
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
#include "two-phase.h"
// #include "src-local/two-phaseVE.h"
// #include "src-local/log-conform-viscoelastic-scalar-2D.h"
#define logFile "log-die-swell_Newt.dat"

#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

#define tsnap (1e-2)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J

#define epsilon (4e-2)
#define R2(x,y,z) (sq(y) + sq(z))


int MAXlevel;
// We -> Weber number
// Oh -> Solvent Ohnesorge number
// Oha -> air Ohnesorge number
// Bo -> Bond number
// De -> Deborah number
// Ec -> Elasto-capillary number
// for now there is no viscoelasticity

double We, Oh, Oha, Bo, De, Ec, tmax;
char nameOut[80], dumpFile[80];

// boundary conditions
u.n[top] = neumann(0.0);
p[top] = dirichlet(0.0);
u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);

f[left] = y > 1. + epsilon ? 0.0 : y < 1. - epsilon ? 1.0 : 0.5 * (1.0 + tanh((1e0 - R2(x,y,z)) / epsilon));
u.n[left] = dirichlet(f[]*2e0*(1-R2(x,y,z)));
u.t[left] = dirichlet(0.0);

int main(int argc, char const *argv[]) {


  // Values taken from the terminal
  L0 = 10.0;
  MAXlevel = 9;
  tmax = 100;
  We = 1e0;
  Oh = 1e-1;
  Oha = 1e-4;
  Bo = 0.0;

  De = 0.0;
  Ec = 0.0;

  if (De != 0.0 || Ec != 0.0) {
    fprintf(ferr, "De and Ec must be zero! This is a Newtonian simulation!\n");
    return 1;
  }

  init_grid (1 << 6);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");


  rho1 = 1., rho2 = 1e-3;
  mu1 = Oh/sqrt(We), mu2 = Oha/sqrt(We);
  // lambda1 = De*sqrt(We), lambda2 = 0.;
  // G1 = Ec/We, G2 = 0.;
  G.x = Bo/We;
  f.sigma = 1.0/We;

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine(x < 2*epsilon && R2(x,y,z) < sq(1+epsilon) && level < MAXlevel);
    fraction (f, 1-R2(x,y,z)-x/epsilon);
    foreach(){
      u.x[] = f[]*2e0*(1-R2(x,y,z));
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
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
      (double[]){fErr, VelErr, VelErr, KErr},
      MAXlevel, 4);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, We %2.1e, Ohs %2.1e, Oha %2.1e, De %2.1e, Ec %2.1e\n", MAXlevel, We, Oh, Oha, De, Ec);
}

/**
## Log writing
*/
event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  scalar pos[];
  position (f, pos, {0,1,0});
  double ymax = statsf(pos).max;
  scalar posX[];
  position (f, posX, {1,0,0});
  double xmax = statsf(posX).max;

  static FILE * fp;
  if (pid() == 0) {
    const char* mode = (i == 0) ? "w" : "a";
    fp = fopen(logFile, mode);
    if (fp == NULL) {
      fprintf(ferr, "Error opening log file\n");
      return 1;
    }

    if (i == 0) {
      fprintf(ferr, "Level %d, We %2.1e, Ohs %2.1e, Oha %2.1e, De %2.1e, Ec %2.1e\n", MAXlevel, We, Oh, Oha, De, Ec);
      fprintf(ferr, "i dt t ke xmax ymax\n");
      fprintf(fp, "Level %d, We %2.1e, Ohs %2.1e, Oha %2.1e, De %2.1e, Ec %2.1e\n", MAXlevel, We, Oh, Oha, De, Ec);
      fprintf(fp, "i dt t ke xmax ymax\n");
    }

    fprintf(fp, "%d %g %g %g %g %g\n", i, dt, t, ke, xmax, ymax);
    fprintf(ferr, "%d %g %g %g %g %g\n", i, dt, t, ke, xmax, ymax);

    fflush(fp);
    fclose(fp);
  }

  if(ke < -1e-10) return 1;
  if (xmax > 0.9*L0){ 
    fprintf(ferr, "The drop has come very close to the end of domain! Stopping simulation!\n");
    return 1;
  }

  if (i > 1e1 && pid() == 0) {
    if (ke > 1e2 || ke < 1e-8) {
      const char* message = (ke > 1e2) ? 
        "The kinetic energy blew up. Stopping simulation\n" : 
        "kinetic energy too small now! Stopping!\n";
      
      fprintf(ferr, "%s", message);
      
      fp = fopen(logFile, "a");
      fprintf(fp, "%s", message);
      fflush(fp);
      fclose(fp);
      
      dump(file=dumpFile);
      return 1;
    }
  }

}

/**
 * @file testing mac mpi
 * @brief This file contains the simulation code for the die-swell of a viscoelastic liquid being extruded out of a die. 
 * @author Vatsal Sanjay
 * @version 2.0
 * These are viscoelastic simulations!
 * @date Oct 22, 2024
*/

#include "navier-stokes/centered.h"
#define logFile "log-die-swell_ViscoElastic.dat"

/**
The default boundary conditions are symmetry (i.e. slip walls). We
need no-slip on three boundaries and $u=1$ on the top
boundary i.e. */

u.t[top] = dirichlet(1);

/**
For the other no-slip boundaries this gives */

u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

char nameOut[80], dumpFile[80];
int main(int argc, char const *argv[]) {

  origin (-0.5, -0.5);
  init_grid (1 << 6);
  mu[] = {1e-3,1e-3};
  DT = 0.1;
  CFL = 0.5;

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    foreach(){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
  }
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += 0.01; t <= 1e1) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Writing Log
*/
event logWriting (i++) {

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
      fprintf(ferr, "i t dt\n");
      fprintf(fp, "i t dt\n");
    }
  }

  fprintf(fp, "%d %g %g\n", i, t, dt);
  fprintf(ferr, "%d %g %g\n", i, t, dt);
  fflush(fp);

}
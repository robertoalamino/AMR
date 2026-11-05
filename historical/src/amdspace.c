// Calculates asymptotic death probabilities
// for the whole range of parameters defining the AMDs


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "genotype.h"
#include "random.h"
#include "response.h"
#include "lattice.h"
#include "config.h"

#ifndef AMR_STEPS
#define AMR_STEPS 3000
#endif

#ifndef AMR_REALISATIONS
#define AMR_REALISATIONS 100
#endif

// ====
// Main
// ====
int main(int argc, char *argv[])
{

  RunConfig config = {
  .seed = -1,
  .steps = AMR_STEPS,
  .realisations = AMR_REALISATIONS,
  .treatment_start = 0,
  .treatment_end = 0
  };

  parse_config(argc, argv, &config, 0);


  int
    cull,               // Cull population? 0 - No | 1 - Yes
    cn,                 // Chosen neighbour
    iaux,iaux2,         // Auxiliary integer variable
    ia,ib,              // AMD loops
    i,j,t,z,            // Loop variables
    L = 50,             // Linear Grid Size
    N = L*L,            // Number of sites in the grid
    T = config.steps,   // Total number of time steps
    D = 50,             // Plasmid dimension
    M,                  // Number of realisations with the same AMD
    no,nno,             // Number of occupied cells
    cc,                 // Chosen cell index
    sc,                 // Cell to spread to
    s[N+1],             // Grid occupation {0,1}
    oc[N+1],            // List of occupied cells (contains their positions)
    noc[N+1],           // New list of occupied cells
    P[N+1][D+1];        // Plasmids {+-1}
  double
    wtp,                // Average pi of wild type
    wtpp,               // Wild type +1 probability
    aa,ab,              // Parameters of AMDs
    r,                  // Reproduction rate
    m,                  // Plasmid mutation probability
    iq,                 // Initial death probability
    c,                  // Threshold concentration
    phi,                // c parameter
    h,                  // Step rigidity
    omega,              // p parameter
    delta,              // Concentration above threshold
    tp,                 // Total death probability
    rc,                 // Antibiotic reference concentration
    cThres,             // Sampling threshold
    cRate,              // Sampling rate
    fapt,               // Final average q
    apt[T+1];           // Average death probability over cells
  long
    seed;
  FILE
    *qav;

  seed = config.seed;

  qav = fopen("asymdrug.txt","w");
  if(qav == NULL)
  {
      perror("Cannot open asymdrug.txt");
      return EXIT_FAILURE;
  }

  // Parameter Initialization
  M = config.realisations;

  iq = 0.3;

  wtp = 0;
  wtpp = (wtp+1.0)/2.0;

  m   = 0.0;
  r   = 1.0;

  cull   = 1;
  cThres = 1.0;
  cRate  = 0.5;

  fprintf(stderr,
        "Program: amdspace\n"
        "Lattice: %d x %d\n"
        "Genotype dimension: %d\n"
        "Steps: %d\n"
        "Realisations per parameter combination: %d\n"
        "Parameter combinations: 121\n"
        "Seed: %ld\n",
        L, L, D, T, M, -seed);

  aa=-1.0;
  for(ia=0;ia<=10;ia++)
  {
    ab=-1.0;
    for(ib=0;ib<=10;ib++)
    {
      // Calculates the reference concentration for the current AMD
      rc = conc(aa,ab,wtp,iq);

      // Realisations Loop
      fapt=0;
      for(z=1;z<=M;z++)
      {
        // Measurements are unavailable until first calculated.
        for(t=0;t<=T;t++) apt[t] = NAN;
        
        // Clearing the grid
        for(i=1;i<=N;i++)
        {
          s[i] = 0;
          oc[i] = 0;

        }

        // Populating the grid with DNA
        // according to the given probability wtpp
        for(i=1;i<=N;i++)
        {
          for(j=1;j<=D;j++)
          {
            if (ran2(&seed)<=wtpp) P[i][j] = 1;
            else P[i][j]=-1;
          }
        }

        // Distributing living cells on the grid uniformly randomly
        no = 0;
        for(i=1;i<=N;i++)
        {
          if(ran2(&seed)>0.5)
          {
            s[i]=1;
            no++;
            oc[no] = i;
          }
        }

        // Cell natural life cycle
        for(t=1;t<=T;t++)
        {

          nno = no;

          if(no!=0) // Go through the whole process if there is a cell alive
          {

            // CELL REPRODUCTION
            nno = no;
            for(i=1;i<=no;i++)
            {
              iaux = (int)(ran2(&seed)*(no-i+1)+1);
              cc = oc[iaux];

              if(ran2(&seed)<=r)
              {
                oc[iaux] = oc[no-i+1];     // Rotates occupation vector
                oc[no-i+1] = cc;

                // Chooses one site to spread
                iaux2 = (int)(4.0*ran2(&seed)+1.0);
                cn = lattice_neighbour(cc, iaux2, L);

                // Checks if neighbour is free and spreads
                if (s[cn]==0)
                {
                  // Creates new cell
                  sc = cn;
                  s[sc]=1;
                  nno++;
                  oc[nno]=sc;

                  // Creates new plasmid with possible mutation
                  for(j=1;j<=D;j++)
                  {
                    if (ran2(&seed)<m) P[sc][j] = -P[cc][j];
                    else P[sc][j] = P[cc][j];
                  }
                }
              }
            }
            no = nno;       // Updating the number of occupied sites

            // AMD DEATH
            if(no!=0)
            {

              for(i=1;i<=no;i++)
              {
                cc = oc[i];         // Chosen cell

                // Threshold Concentration
                phi = aa*total(P[cc],D)/D;
                c = (1.0+phi)/(1.0-phi);

                delta = rc-c;


                // Threshold hardness
                omega = ab*total(P[cc],D)/D;
                h = (1.0+omega)/(1.0-omega);

                // Cell death probability
                tp = 0.5*(1.0+tanh(h*delta));

                if(ran2(&seed)<=tp) s[cc]=0;
              }

              nno=0;
              for(i=1;i<=no;i++)
              {
                if(s[oc[i]])
                {
                  nno++;
                  noc[nno]=oc[i];
                }
              }

              no=nno;
              for(i=1;i<=no;i++) oc[i]=noc[i];

            } // If statement for AMD death

            // AVERAGE VALUES OF PARAMETERS
            if (no!=0)
            {
              // Calculating average values
              apt[t]=0;

              for(i=1;i<=no;i++)
              {
                phi = aa*total(P[oc[i]],D)/D;
                c = (1.0+phi)/(1.0-phi);

                omega = ab*total(P[oc[i]],D)/D;
                h = (1.0+omega)/(1.0-omega);

                delta  = rc-c;

                apt[t]+=0.5*(1.0+tanh(h*delta));

              }

              apt[t]/=(float)no;

            }
            else  // If there are no cells alive...
            {
              apt[t]= apt[t-1];
            }
          } // End of main if-else checking whether some cell is alive
          else  // If there are no cells alive...
          {
            apt[t]= apt[t-1];
          }

          // Culling cell population if necessary
          while(cull&&(((float)nno/(float)N)>=cThres))
          {

            for(i=1;i<=nno;i++)
            {
              if (ran2(&seed)<cRate)
              {
                s[oc[i]]=0;      // Kills the cell
                oc[i]=oc[nno];
                oc[nno]=0;
                nno--;    // Decreases the number of living cells

              }
            }
          } // while loop for sampling

          no=nno;

        } // Cell Life Cycle Loop - t


        fapt+=apt[T];
      } // Realisations Loop - z

      fprintf(qav,"%f\t%f\t%f\n",aa,ab,fapt/(float)M);

      ab+=0.2;
    }  // End of ib loop
    aa+=0.2;
  } // End of ia loop

  if(fclose(qav) == EOF)
  {
    perror("Cannot finish writing asymdrug.txt");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}

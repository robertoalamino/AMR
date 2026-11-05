// Single AMD treatment
// D -> oo

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "genotype.h"
#include "random.h"
#include "response.h"
#include "lattice.h"
#include "config.h"
#include "population.h"

// Default experiment settings -- override with compiler -D options.
#ifndef AMR_STEPS
#define AMR_STEPS 15000
#endif

#ifndef AMR_TREATMENT_START
#define AMR_TREATMENT_START 5000
#endif

#ifndef AMR_TREATMENT_END
#define AMR_TREATMENT_END 5500
#endif

#ifndef AMR_REALISATIONS
#define AMR_REALISATIONS 1000
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
  .treatment_start = AMR_TREATMENT_START,
  .treatment_end = AMR_TREATMENT_END
  };

  parse_config(argc, argv, &config, 1);

  int
    iaux,                         // Auxiliary integer variable
    i,j,t,z,                      // Loop variables
    useamd,                       // Use AMD? 1 - Yes | 0 - No
    cull,                         // Cull population? 1 - Yes | 0 - No
    L = 50,                       // Linear Grid Size
    N = L*L,                      // Number of sites in the grid
    D = 50,                       // DNA dimension
    T = config.steps,             // Total number of time steps
    T1 = config.treatment_start,  // Beginning of treatment
    T2 = config.treatment_end,    // End of treatment
    M,                            // Number of realisations with the same AMD
    no,nno,                       // Number of occupied cells
    extinction_time,              // 0 = initially empty; -1 = no extinction observed
    cc,                           // Chosen cell index
    cn,                           // Chosen neighbour
    sc,                           // Cell to spread to
    s[N+1],                       // Grid occupation: 1 - Occupied | 0 - Empty
    oc[N+1],                      //  List of occupied cells (contains their positions)
    noc[N+1],                     // New list of occupied cells
    P[N+1][D+1];                  // Cell DNA
  double
    wtp,                          // Average phenotype of wild type
    aai,abi,                      // Environment
    aaf,abf,                      // AMD
    d,                            // Natural death rate
    r,                            // Reproduction rate
    m,                            // Mutation probability per base
    iq,                           // Initial death probability
    c,                            // Threshold concentration
    phi,                          // c parameter
    h,                            // Step rigidity
    omega,                        // p parameter
    delta,                        // Concentration above threshold
    tp,                           // Total death probability
    rci,                          // Environment reference concentration
    rcf,                          // AMD reference concentration
    cThres,                       // Sampling threshold
    cRate,                        // Sampling rate
    apt[T+1];                     // Average death probability over cells
  long
    seed;
  FILE
    *population,
    *qav,
    *extinction;

  seed = config.seed;

  population = fopen("population.txt", "w");
  if(population == NULL)
  {
      perror("Cannot open population.txt");
      return EXIT_FAILURE;
  }

  qav = fopen("qav.txt", "w");
  if(qav == NULL)
  {
      perror("Cannot open qav.txt");
      fclose(population);
      return EXIT_FAILURE;
  }

  extinction = fopen("extinction.txt", "w");
  if(extinction == NULL)
  {
    perror("Cannot open extinction.txt");
    fclose(population);
    fclose(qav);
    return EXIT_FAILURE;
  }

  // Parameter Initialization
  M = config.realisations;

  iq = 0.3;

  m   = 0.001;
  d   = 0.0;
  r   = 1.0;

  cull   = 1;
  cThres = 1.0;
  cRate  = 0.5;

  fprintf(stderr,
        "Program: reversaltimes\n"
        "Lattice: %d x %d\n"
        "Genotype dimension: %d\n"
        "Steps: %d\n"
        "Realisations: %d\n"
        "Treatment interval: [%d, %d)\n"
        "Seed: %ld\n",
        L, L, D, T, M, T1, T2, -seed);

  // Realisations Loop
  for(z=1;z<=M;z++)
  {
    // Measurements are unavailable until first calculated.
    for(t=0;t<=T;t++) apt[t] = NAN;

    useamd = 0;

    // Treatment reference concentration is unavailable before treatment.
    rcf = NAN;

    // Preparing the grid
    for(i=1;i<=N;i++)
    {
      s[i] = 0;                           // Clearing the grid
      oc[i] = 0;                          // No occupied cells
      for(j=1;j<=D;j++)
        if (ran2(&seed)<=.5) P[i][j] = 1;
        else P[i][j]=-1;
    }

    // Distributing living cells on the grid randomly
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

    extinction_time = no == 0 ? 0 : -1;

    // Choosing the initial environment
    aai = ran2(&seed)*2.0-1.0;
    abi = ran2(&seed)*2.0-1.0;
    rci = conc(aai,abi,0,iq);

    // Choosing the AMD
    aaf = ran2(&seed)*2.0-1.0;
    abf = ran2(&seed)*2.0-1.0;

    // Cell natural life cycle
    fprintf(population,"%d\t",no);

    for(t=1;t<=T;t++)
    {
      nno = no;
      if(t==T1)
      {
        useamd = 1;
        wtp=0;
        for(i=1;i<=N;i++)
          for(j=1;j<=D;j++)
            wtp+=P[i][j];
        wtp=wtp/(float)(N*D);
        rcf = conc(aaf,abf,wtp,iq);
      }
      if(t==T2) useamd = 0;

      if(no!=0) // Go through the whole process if there is a cell alive, otherwise just calculate the probability
      {

        // CELL REPRODUCTION
        nno = no;                   // New number of occupied cells
        for(i=1;i<=no;i++)
        {
          iaux = (int)(ran2(&seed)*(no-i+1)+1);
          cc = oc[iaux];

          // Remove this parent from the remaining selection pool.
          oc[iaux] = oc[no-i+1];
          oc[no-i+1] = cc;

          if(ran2(&seed)<=r)
          {
            // Chooses one site to spread
            int empty[4];
            int empty_count = lattice_empty_neighbours(cc, L, s, empty);

            // Reproduce only when space is available.
            if(empty_count > 0)
            {
              int choice = (int)(ran2(&seed) * empty_count);
              cn = empty[choice];
              // Create new cell
              sc = cn;
              s[sc]=1;
              nno++;
              oc[nno]=sc;

              // Create new plasmid with possible mutation
              for(j=1;j<=D;j++)
                if (ran2(&seed)<m) P[sc][j] = -P[cc][j];
                else P[sc][j] = P[cc][j];
            }
          }
        }
        no = nno;       // Updating the number of occupied sites

        // CELL NATURAL DEATH
        if (no!=0)  // (Only if there's a living cell)
        {
          for(i=1;i<=no;i++)
            // Cell's death rate

            if (ran2(&seed)<d) s[oc[i]]=0;       // Kills the cell

          nno=0;
          for(i=1;i<=no;i++)
          {
            if(s[oc[i]])
            {
              nno++;
              noc[nno]=oc[i];
            }
          }
        }
        // AMD DEATH - standard environment
        if(no!=0)
        {

          for(i=1;i<=no;i++)
          {
            cc = oc[i];         // Chosen cell

            // Threshold Concentration
            phi = aai*total(P[cc],D)/D;
            c = (1.0+phi)/(1.0-phi);

            delta = rci-c;

            // Maximum death rate
            omega = abi*total(P[cc],D)/D;
            h = (1.0+omega)/(1.0-omega);

            // Cell death probability
            tp = .5*(1.0+tanh(h*delta));

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

        } // If statement for AMD death - standard

        // AMD DEATH - antibiotic
        if(no!=0&&useamd)
        {

          for(i=1;i<=no;i++)
          {
            cc = oc[i];         // Chosen cell

            // Threshold Concentration
            phi = aaf*total(P[cc],D)/D;
            c = (1.0+phi)/(1.0-phi);

            delta = rcf-c;

            // Maximum death rate
            omega = abf*total(P[cc],D)/D;
            h = (1.0+omega)/(1.0-omega);

            // Cell death probability
            tp = .5*(1.0+tanh(h*delta));

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

        } // If statement for AMD death - antibiotic

        // AVERAGE VALUES OF PARAMETERS
        if(t>=T1)
        {

          if (no!=0)
          {
            // Calculating average values
            apt[t]=0;

            for(i=1;i<=no;i++)
            {
              phi = aaf*total(P[oc[i]],D)/D;
              c = (1.0+phi)/(1.0-phi);

              h = (1.0+abf*total(P[oc[i]],D)/D)/(1.0-abf*total(P[oc[i]],D)/D);
              delta  = rcf-c;
              apt[t]+=.5*(1.0+tanh(h*delta));
            }

            apt[t]/=(float)no;
            fprintf(qav,"%f\t",apt[t]);
          }
          else  // If there are no cells alive...
          {
            apt[t]= apt[t-1];
            fprintf(qav,"%f\t",apt[t]);
          }
        }
      } // End of main if-else checking whether some cell is alive
      else  // Population was already extinct at the start of this step.
      {
        apt[t]= apt[t-1];

        if(t>=T1)
          fprintf(qav,"%f\t",apt[t]);
      }

      fprintf(population,"%d\t",no);

      // Culling cell population if necessary
      // Repeat passes if the population remains at the threshold.
      while(cull && (double)nno / N >= cThres)
        nno = population_cull(s, oc, nno, cRate, &seed);

      no=nno;

      if(no == 0 && extinction_time == -1)
        extinction_time = t;

    } // Cell Life Cycle Loop - t

    fprintf(population,"\n");
    fprintf(qav,"\n");
    fprintf(extinction, "%d\t%d\n", z, extinction_time);

  } // Realisations Loop - z

  int status = EXIT_SUCCESS;

  if(fclose(population) == EOF)
  {
    perror("Cannot finish writing population.txt");
    status = EXIT_FAILURE;
  }

  if(fclose(qav) == EOF)
  {
    perror("Cannot finish writing qav.txt");
    status = EXIT_FAILURE;
  }

  if(fclose(extinction) == EOF)
  {
    perror("Cannot finish writing extinction.txt");
    status = EXIT_FAILURE;
  }

  return status;
}

// Calculates asymptotic death probabilities
// for the whole range of parameters defining the AMDs


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ran2 definitions
// ----------------
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// ==============================
// Random number generator
// Uniform on (0,1)
// Source: Numerical Recipes in C
// ==============================
double ran2(long *idum)
{
  int
    j;
  long
    k;
  static long
    idum2=123456789;
  static long
    iy=0;
  static long
    iv[NTAB];
  double
    temp;

  if (*idum<=0)
  {
    if(-(*idum)<1) *idum=1;
    else *idum=-(*idum);
    idum2=(*idum);
    for(j=NTAB+7;j>=0;j--)
    {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum<0) *idum+=IM1;
      if (j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }

  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum<0) *idum+=IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2<0) idum2+=IM2;

  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j]=*idum;
  if (iy<1) iy+=IMM1;

  if ((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
}

// ========================
// Random Integer Generator
// ========================
int randi(int max, long *idum)
{
  int
    r;

  r = (int)(max*ran2(idum))+1;

  return r;
}

// ===============================
// Function: Hyperbolic Arctangent
// ===============================
double atanh(double x)
{
  double r;

  r=.5*log((1.0+x)/(1.0-x));

  return r;
}

// ===============
// Function: total
// ===============
double total(int *v, int n)
{
  int
    i;
  double
    r;

  r=0;
  for(i=1;i<=n;i++) r+=v[i];

  return r;
}

// =======================
// Function: concentration
// =======================
double conc(double a, double b, double p,double q)
{
  double
    r;

  r = (1.0+a*p)/(1.0-a*p)+(1.0-b*p)*atanh(2*q-1.0)/(1.0+b*p);

  return r;
}

// ====
// Main
// ====
int main(void)
{
  int
    cull,               // Cull population? 0 - No | 1 - Yes
    cn,                 // Chosen neighbour
    iaux,iaux2,         // Auxiliary integer variable
    ia,ib,              // AMD loops
    i,j,t,z,            // Loop variables
    L = 50,             // Linear Grid Size
    N = L*L,            // Number of sites in the grid
    T = 3000,           // Total number of time steps
    D = 50,             // Plasmid dimension
    M,                  // Number of realisations with the same AMD
    no,nno,             // Number of occupied cells
    cc,                 // Chosen cell index
    sc,                 // Cell to spread to
    fb,                 // Bit to be flipped in mutations
    s[N+1],             // Grid occupation {0,1}
    oc[N+1],            // List of occupied cells (contains their positions)
    noc[N+1],           // New list of occupied cells
    P[N+1][D+1];        // Plasmids {+-1}
  double
    wtp,                // Average pi of wild type
    wtpp,               // Wild type +1 probability
    aa,ab,              // Parameters of AMDs
    aux,                // Auxiliary real variable
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

  seed=-1;

  qav = fopen("C:\\Users\\alaminrc\\Dropbox\\Work\\AMR\\data\\asymdrug.txt","w");

  // Parameter Initialization
  M  = 100;

  iq = 0.3;

  wtp = 0;
  wtpp = (wtp+1.0)/2.0;

  m   = 0.0;
  r   = 1.0;

  cull   = 1;
  cThres = 1.0;
  cRate  = 0.5;

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
                switch(iaux2)
                {
                  case 1:
                    if (cc%L==0) cn = cc-L+1;           // Right neighbour
                    else cn = cc+1;
                    break;
                  case 2:
                    if ((cc-1)%L==0) cn = cc+L-1;       // Left neighbour
                    else cn = cc-1;
                    break;
                  case 3:
                    if (cc-L<1) cn = cc-L+N;            // Upper neighbour
                    else cn=cc-L;
                    break;
                  case 4:
                    if (cc+L>N) cn = cc+L-N;            // Down neighbour
                    else cn=cc+L;
                    break;
                }

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

}

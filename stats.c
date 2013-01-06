/** 
  @file stats.c
  Accumulate statistical information during the run of the program.
  They will be printed onto the screen after every block (that is a certain number of
  computational steps) and a summary will be presented after every run.
  Those statistics are not saved to any file.
 */

#include "Gillespie.h"
#include "stats.h"
#include "queue.h"

/*------------------------Globally defined variables------------------------*/  

Stats *Xblk, ///< Array for all components with statistical information in current block
      *Xrun; ///< Array for all components with statistical information in total run
double tau_step, tau_out;


/**
  Initialize statistical data of the block with zero
*/
void block_init()
{
  int i;

  tau_out = tau_step = sys.tau;
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    Xblk[i].sum   = 0.;
    Xblk[i].sumsq = 0.;
    Xblk[i].acc   = 0;
    Xblk[i].err   = 0.;
    Xblk[i].noise = 0.;
  }
}


/**
  Accumulates statistical information
*/
void block_acc()
{
  int i;
  double dt;

  // calculate the time difference to the last accumulation
  dt = sys.tau - tau_step;
  if( sys.output_conc )
    dt /= sys.volume;
  tau_step = sys.tau;
  
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    // accumulate values for each component
    Xblk[i].sum   += X[i] * dt;
    Xblk[i].sumsq += X[i] * X[i] * dt;
    Xblk[i].acc++;
  }
}


/**
  Finalizes a block and prints output onto the screen.
  @param[in]  block  ID of the block
*/
void block_finish(int block)
{
  int i;
  double Xtot = 0.0, dt;

  // calculate the time difference to the last write out
  dt = sys.tau - tau_out;
  tau_out = sys.tau;
  
  for( i=0; i<sys.Ncomp; i++ ) 
  {

    // average the count (dividing the integral by the total time)
    Xblk[i].sum   /= dt;
    // average the count squared (dividing the integral by the total time)
    Xblk[i].sumsq /= dt;

    // calculate the error sigma^2 using the above values
    Xblk[i].err    = Xblk[i].sumsq - Xblk[i].sum * Xblk[i].sum; /*sigma^2*/
    // calculate the standard deviation, if possible
    if (Xblk[i].err > 0.) {
      Xblk[i].err = sqrt(Xblk[i].err); /*sigma*/
      // calculate the noise
      if (Xblk[i].sum != 0.) 
        Xblk[i].noise = Xblk[i].err / Xblk[i].sum;
      else
        Xblk[i].noise = 0.;
    }
    else { // error calculation not possible
      Xblk[i].noise = 0.;
    }
  }

  // transfer statistics of the block to that of the total run
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    Xrun[i].sum   += Xblk[i].sum;
    Xrun[i].sumsq += Xblk[i].sum * Xblk[i].sum;
    Xrun[i].noise += Xblk[i].noise;
    Xrun[i].acc   ++;
  }

  // print statistics about the run of the current block
  printf("\n\nThe results of block %d\n\n",block);
  printf("Elapsed time         %8.4f\n\n",dt);
  
  printf("Component\tblock average\tdeviation\t   noise\n");
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    printf("%s\t\t%8.4f\t%8.4f\t%8.4f\n",
           Xname[i],
           Xblk[i].sum,
           Xblk[i].err,
           Xblk[i].noise
          );
    Xtot += Xblk[i].sum;
  }
  printf( "Total:\t\t%8.4f\n", Xtot );
  printf( "INFO:\tQueuelength: %d  \tTotal_tau: %8.4f\tVolume: %8.4f\n", 
          queue_length(), sys.tau, sys.volume
        );
}


/**
  Initialize statistical variables of the run with zero.
*/
void run_init()
{
  int i;

  sys.tau = sys.tau_init; // initialize run-time with previous times
  printf("INIT: %f - %f", sys.tau, sys.tau_init );
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    Xrun[i].sum   = 0.;
    Xrun[i].sumsq = 0.;
    Xrun[i].acc   = 0;
    Xrun[i].err   = 0.;
    Xrun[i].noise = 0.;
  }
  
  block_init();
}


/**
  Output the statistical information for the complete run
  @param[in]  run   ID of the current run
 */
void run_finish(int run)
{
  int  i;
 
  printf("\n==============================================================================\n\n");
  if (run==EQUIL) 
    printf("The averages of the equilibration run.\n\n");
  else
    printf("The average of the production run.\n\n");

  printf("Elapsed time         %8.4f\n\n",sys.tau);
  
  printf("Run averages.\n\n");
  printf("Component\tRun average\t   error\t   noise\t\t#counts\n");
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    Xrun[i].sum   /= Xrun[i].acc;
    Xrun[i].sumsq /= Xrun[i].acc;
    Xrun[i].err    = (Xrun[i].sumsq - Xrun[i].sum * Xrun[i].sum) / 
        (double) Xrun[i].acc;
    if (Xrun[i].err > 0.)  Xrun[i].err = sqrt(Xrun[i].err);
    Xrun[i].noise /= Xrun[i].acc;
    printf("%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n",
           Xname[i],Xrun[i].sum,Xrun[i].err,Xrun[i].noise,Xrun[i].acc);
  }
}


/**
  Write statistics of the run to an output file
*/
void run_acc()
{
  int   i;
  Stats stats_tp; //< struct containing the statistical information 
  FILE  *fp;
  char  filename[40];

  sprintf(filename,"%s.so",sys.name);
  fp = fopen(filename,"a");
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    // average the count (dividing the integral by the number of blocks)
    stats_tp.sum   = Xrun[i].sum / Xrun[i].acc;
    // average the count squared (dividing the integral by the number of blocks)
    stats_tp.sumsq = Xrun[i].sumsq / Xrun[i].acc;
    stats_tp.err   = (stats_tp.sumsq - stats_tp.sum * stats_tp.sum) / 
      (double) Xrun[i].acc;
    if (stats_tp.err > 0.)  stats_tp.err = sqrt(stats_tp.err);
    stats_tp.noise = Xrun[i].noise / Xrun[i].acc;

    fprintf( fp, "%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n",
             Xname[i], stats_tp.sum, stats_tp.err, stats_tp.noise, stats_tp.acc
           );

  }
  fclose(fp);
}


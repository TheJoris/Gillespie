/** 
  @file analyse.c
  Gathers results during the calculation.
  At certain points during the calculation the state of the system is saved. It is written
  to the hard disk after the run has finished. The data is not shown on the screen.
 */

#include "Gillespie.h"
#include "analyse.h"
#include "queue.h"

/**
  Structure for gathering data in each step
*/
typedef struct analyseData_type {
  double tau,      ///< the elapsed time
         *X,       ///< array of the concentrations of the components
         volume;   ///< volume of the simulation
  long long int ns,          ///< number of total steps 
      queuelength; ///< length of the queue
} AnalyseData;

/*------------------------Globally defined functions-------------------------*/

AnalyseData *adata=NULL; ///< global array storing the statistical information
unsigned int adata_id=0,         ///< ID which saves the current position, where data is going to be written to
     adata_length=0;     ///< current length of the adata array
double adata_time=0.;    ///< the last time a data set has been written

/**
  Initialize analyzation by allocating memory.
*/
void analyse_init()
{
  int i;
  
  // reset ID to ensure array is filled from the beginning onwards
  adata_id=0;
  adata_time = sys.tau_init - sys.dt;
  
  // initial container for data
  adata_length = ANALYSIS_EXTEND;
  adata = malloc( adata_length * sizeof(AnalyseData) );
  for( i=adata_length-ANALYSIS_EXTEND; i<adata_length; i++ ) 
  {
    adata[i].tau = 0.;
    adata[i].X = malloc( sys.Ncomp * sizeof(double) );
  }
}

/**
  Keep data of current step for later analysis.
  @param[in]  steps     Number of current step
*/
void analyse( long long int steps )
{
  int  i;
  
  // check, if actual writeout has been earlier
  while( adata_time <= sys.tau )
  {
    // check, if enough memory has been allocated
    if( adata_id >= adata_length )
    {
      // allocate memory for the data structure
      adata_length += ANALYSIS_EXTEND;
      adata = realloc( adata, adata_length * sizeof(AnalyseData) );
      // initialize values
      for( i=adata_length-ANALYSIS_EXTEND; i<adata_length; i++ )
      {
        adata[i].tau = 0.;
        adata[i].X = malloc( sys.Ncomp * sizeof(double) );
      }
    }
    
    // step to next time
    adata_time += sys.dt;
    
    // save important information
    adata[adata_id].tau = adata_time;
    adata[adata_id].ns = steps;
    adata[adata_id].queuelength = queue_length();
    adata[adata_id].volume = sys.volume;
    
    // copy concentrations
    for( i=0; i<sys.Ncomp; i++ )
    {
      if( sys.output_conc == OUTPUT_COPYNUMBER )
        adata[adata_id].X[i] = (double) X[i]; // convert to double to be compatible with concentration output
      else if( sys.output_conc == OUTPUT_CONCENTRATION )
        adata[adata_id].X[i] = ( (double) X[i] ) / sys.volume;
    }

//     memcpy( &(adata[adata_id].X), &X, sys.Ncomp * sizeof(int) );

    // step counter
    adata_id++;
  }
}


#define OUTPUT_PHOS filename_phos != NULL
#define OUTPUT_TOTAL filename_total != NULL
/**
  Does special output used for the Kai system.
  If any of the given filenames is set to NULL, the respective output is not made.
  @param[in]  filename_total    Filename where the total amount of KaiC is written to
  @param[in]  filename_phos     Filename where the phosphorylation ratio of KaiC is written to
  @warning This function relies on a special scheme for the naming of the components:
    All KaiC complexes which should go into the statistics have to start with C#
    where the # stands for the amount of phosphorylation of KaiC.
*/
void output_special( char *filename_total, char *filename_phos )
{
  FILE *fp_phos, *fp_total;
  int *w;
  unsigned int i, j;
  double total=0, phos=0, norm=0;

  w = malloc( sys.Ncomp * sizeof(unsigned int) );
  
  // get weigths for all components
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    if( Xname[i][0] == 'C' ) // condition for a component to be KaiC
    {
      if( OUTPUT_PHOS )
        sscanf( &Xname[i][1], "%d", &w[i] ); // phosphorylation number has to be the second character
      else
        w[i] = 1;
      
      if( w[i] > norm ) // find maximal phosphorylation count
        norm = w[i];
    }
    else
    {
      w[i] = -1; // ignore this component
    }
  }
  norm = 1./norm;
  
  if( OUTPUT_TOTAL )
  {
    fp_total = fopen( filename_total, "w" );
    fprintf( fp_total, "#time\t\t" );
    fprintf( fp_total, "Total KaiC\t" );
    fprintf( fp_total, "step\n" );
  }
  if( OUTPUT_PHOS )
  {
    fp_phos = fopen( filename_phos, "w" );
    fprintf( fp_phos, "#time\t\t" );
    fprintf( fp_phos, "Phos. Ratio\t" );
    fprintf( fp_phos, "step\n" );
  }
  
  // iterate through all steps
  for( j=0; j<adata_id; j++ ) 
  {
    total = 0.;
    if( OUTPUT_TOTAL ) 
      fprintf( fp_total, "%e\t", adata[j].tau );
    if( OUTPUT_PHOS ) 
    {
      phos = 0.;
      fprintf( fp_phos, "%e\t", adata[j].tau );
    }
    
    // iterate through components
    for( i=0; i<sys.Ncomp; i++ )
    {
      if( w[i] >= 0 )
      {
        total += adata[j].X[i];
        if( OUTPUT_PHOS )
          phos += w[i]*adata[j].X[i];
      }
    }
    if( OUTPUT_TOTAL ) 
    {
//      if( sys.output_conc == OUTPUT_CONCENTRATION )
//        total /= sys.volume;
      fprintf( fp_total, "%e\t", total );
      fprintf( fp_total, "%lld\n", adata[j].ns );
    }
    if( OUTPUT_PHOS )
    {
      fprintf( fp_phos, "%e\t", phos/total*norm );
      fprintf( fp_phos, "%lld\n", adata[j].ns );
    }
  }
  
  if( OUTPUT_TOTAL )
    fclose(fp_total);
  if( OUTPUT_PHOS )
    fclose(fp_phos);
  
  free(w);
}
#undef OUTPUT_PHOS
#undef OUTPUT_TOTAL


/**
  Write out the gathered data into individual files for each component.
  @param[in]  run   ID number of the run
*/
void analyse_finish( int run )
{
  int i,j;
  char filename[128], filenameCp[128];
  FILE *fp;
  
  // check, if we have to do the special output
  if( sys.output_kaic )
  {
    sprintf( filename, "%s.%d.Ct", sys.name, run );
    if( sys.output_phos )
    {  
      sprintf( filenameCp, "%s.%d.Cp", sys.name, run );
      output_special( filename, filenameCp );
    }
    else
    {
      output_special( filename, NULL );
    }
  }
  else
  {
    if( sys.output_phos )
    {
      sprintf( filenameCp, "%s.%d.Cp", sys.name, run );
      output_special( NULL, filenameCp );
    }
  }
    // ordinary output is requested
  
  if( ! sys.output_quite )
  {
    // write file for every component
    for( i=0; i<sys.Ncomp; i++ ) 
    {
      if( sys.species == NULL || 0 == strcmp( Xname[i], sys.species ) )
      {
        // open file for writing
        sprintf( filename, "%s.%d.%s", sys.name, run, Xname[i] );
        fp = fopen( filename, "w" );
    
        // write header
        fprintf( fp, "#time\t\t" );
        fprintf( fp, "%s\t", Xname[i] );
        fprintf( fp, "step\n" );
    
        // iterate through all steps
        for( j=0; j<adata_id; j++ ) 
        {
          fprintf( fp, "%e\t", adata[j].tau );
          fprintf( fp, "%e\t", adata[j].X[i] );
          fprintf( fp, "%lld\n", adata[j].ns );
        }
    
        // close file
        fclose(fp);
      }
    }
  }

  // output queue if it has been used
  if( sys.needs_queue )
  {
    // open file for writing the queue
    sprintf(filename,"%s.%d.queue",sys.name,run);
    fp = fopen(filename,"w");
  
    // write header
    fprintf(fp,"#time\t\tqueue\tstep\n");
  
    // iterate through all steps
    for( j=0; j<adata_id; j++ ) 
    {
      fprintf( fp, "%e\t", adata[j].tau );
      fprintf( fp, "%lld\t", adata[j].queuelength );
      fprintf( fp, "%lld \n", adata[j].ns );
    }
    
    // close file
    fclose(fp);
  }

  // output volume if it has been changed
  if( HAS_GROWTH )
  {
    // open file for writing the queue
    sprintf( filename, "%s.%d.volume", sys.name, run );
    fp = fopen( filename, "w" );
  
    // write header
    fprintf(fp,"#time\tvolume\tstep\n");
  
    // iterate through all steps
    for( j=0; j<adata_id; j++ ) 
    {
      fprintf( fp, "%e\t", adata[j].tau );
      fprintf( fp, "%e\t", adata[j].volume );
      fprintf( fp, "%lld\n", adata[j].ns );
    }
    
    // close file
    fclose(fp);
  }
}

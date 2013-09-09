/** 
  @file propagate.c
  central calculations and definition of the algorithm
 */

#include "Gillespie.h"
#include "analyse.h"
#include "random.h"
#include "stats.h"
#include "propagate.h"
#include "queue.h"

/*------------------------Locally defined functions-------------------------*/

void calculate_components();
void determine_propensity_functions();
double ran_gaussian( const double sigma );
void select_reaction( int *j, double sum_a );
void update_concentrations( int j );
void update_concentrations_queuereactants( int j );
void update_concentrations_queueproducts( int j );
void update_propensity_functions( int *ids, int len );
void calculate_volume();

void growth_init();
void growth_step();
void growth_dependent_reactions();

/*------------------------Locally defined variables-------------------------*/

double sum_a,         ///< total of the propensity function over all reaction channels
       volume_min,    ///< minimal volume, which should be reached right after division
       growth_const,  ///< constant which is precalculated for the growth
       growth_dt,     ///< time between two successive volume increases
       ccycle_duptime;///< Time of next gene duplication in current cell cycle.

int gdrl;			  //length of the gdr array == number of volume dependent reactions

    
/**
  Do central calculations.
  This is the core function of the program having most of the logic for each time step.
  A descriptive explanation of the algorithm can be found in the documentation.
  @param[in]  run           ID number of the run
  @param[in]  total_time    The total time, the run should last
  @param[in]  total_steps   The total number of steps, the run should last
*/
void run( int run, double total_time, long total_steps )
{
  int    react, //< next reaction
         i, j;
  long long int   steps; //< number of reaction steps
  double sum_a, //< sum of all propensity functions
         last_volume_update = 0., //< save the last time the volume has updated
         time; //> temporary variable

  // initialize variables
  sys.tau = sys.tau_init;
  sys.current_gene_copynbr = sys.init_gene_copynbr;
  steps = 0;
  analyse_init();
  if( sys.output_stats )
    block_init(); // in stats.c
  
  // initialize growth, if neccessary
  if( HAS_GROWTH )
  {  
    growth_dependent_reactions();
	  growth_init();
 	  last_volume_update = sys.tau - 2.*growth_dt; // make sure, volume is set before first calculation
  }
  else
  {
    last_volume_update = INF_POS;
    growth_dt = INF_POS;
  }
  
  // calculate the propensity function for each reaction channels (stored in a[])
  determine_propensity_functions();

  // run until final time or final number of steps is reached
  while( sys.tau < sys.tau_init + total_time && steps < total_steps )
  {   
    // check if volume or gene copy number should be updated
    if( HAS_GROWTH && sys.tau >= last_volume_update + growth_dt )
    {
      // Store time of this volume update.
      last_volume_update = sys.tau;
    
      // Update volume or #gene and divide if necessary,
	    // influenced propensities are updated.
      growth_step();
    }
    
    // advance the step counter and save current status
    steps++;
    analyse( steps );
    if( sys.output_stats )
      block_acc();
    
    // do statistics and output
    if( sys.stats_steps > 0 && 0 == steps % sys.stats_steps )
    {
      if( sys.output_stats )
      {
        block_finish( steps / sys.stats_steps );
        block_init();
      }
      else
      {
        printf( "%8.4f time elapsed after %'lld steps. (queuelength: %d, volume: %f)\n", sys.tau, steps, queue_length(), sys.volume );
      }
    }
    
    // calculate sum of the propensity functions
    sum_a = 0.;
    for( i=0; i<sys.Nreact; ++i )
      sum_a += a[i];

    // determine the time tau at which the next reaction takes place
    if( sum_a > 0. )
    {
      sys.tau -= log( 1. - ran_get() ) / sum_a;
      // 1.-ran_get() creates a random variable out of (0,1] instead of [0,1)
      
      // check, if a queue has ended in the last step and return the reaction
      react = queue_check( &sys.tau );
    }
    else // no reactions might occur
    {
      // try to step to next item in queue
      react = queue_advance( &sys.tau );
      
      if( react < 0 ) // no item found ==> print error, since algorithm cannot advance
      {
        printf( "Not a single reaction can occur.\n" );
        printf( "The run will be terminated.\n" );
        printf( "sum_a is %f\n", sum_a );
        if( HAS_GROWTH )
          printf( "The volume is %f\n", sys.volume );
        if( sys.needs_queue )
          printf( "The length of the queue is %d\n", queue_length() );
        for( i=0; i<sys.Ncomp; ++i  ) 
          printf( "%d: [%s] is %d\n", i, Xname[i], X[i] );
        abort();
      }
    }
    
    // check if the last time step was too large to account for growth
    if( sys.tau >= last_volume_update + growth_dt )
    {
      // if we got an item from the queue push it back up
      if( react >= 0 ) 
        queue_push( react, sys.tau );

      // make a smaller timestep to account for growth
      sys.tau = last_volume_update + growth_dt + EPSILON;
    }
    else // otherwise: do reaction
    {
      // check if an reaction has already be determined (came off queue)
      if( react >= 0 ) 
      {
        update_concentrations_queueproducts( react );
      }
      else // no reaction determined => draw Poissonian one
      {

        // select next reaction on a random basis
        select_reaction( &react, sum_a );
      
        if( REAC_DELAYED == R[react].type ) // delayed reaction
        {
          // draw gaussian time for this reaction
          do
            time = R[react].time + ran_gaussian( R[react].sigma );
          while( time <= EPSILON ); // make sure that time does not lie in the past

          // push reaction onto the queue
          queue_push( react, time + sys.tau );
          update_concentrations_queuereactants( react );
        }
        else // ordinary reaction
        {
          update_concentrations( react );
        }

      }
    }
  }

  // output statistics about the whole run
  analyse_finish(run);
  //run_finish(run);
  
  // shift origin of time for the next run to current end
  sys.tau_init = sys.tau;
}


/**
  Determine propensity functions for all reaction channels and store their value in a[i].
  Makes sure that all propensity functions are up to date.
*/
void determine_propensity_functions()
{
  int i, *ids;
  
  // setup array with all ids of all reactions
  ids = malloc( sys.Nreact * sizeof( int ) );
  for( i=0; i<sys.Nreact; ++i )
    ids[i] = i;
  
  //calculate_components();						//knocked out
    
  // calculate all propensity functions
  update_propensity_functions( ids, sys.Nreact );
  
  free( ids );
}


/**
  Performs the calculation of components, if they are defined in such a way
*/
void calculate_components()
{
  int id, j;
  
  for( id=0; id<sys.Ncomp; ++id )
  {
    if ( Xconst[id] == 2 )
    {
      X[id] = 0;
      for( j=0; j<Xcalc_count[id]; ++j )
        X[id] += X[Xcalc[id][j]];
    }
  }
}


/**
  Determines selected propensity functions.
  Updates propensity functions for all reactions which IDs are given in an array.
  Special care has to be taken w.r.t. to the volume dependence of the reaction constants.
  @param[in]  ids   A pointer to an array of the reaction ids, which should be updated
  @param[in]  len   The length of the array
*/
void update_propensity_functions( int *ids, int len )				
{
  int i;
  double x;
  React *r;

  // calculate all propensity function which change
  for( i=0; i<len; ++i ) 
  {

    // store pointer to reaction struct
    r = &R[ids[i]];
    
    // check if the reaction constant is definied by a Hill-function
    if( r->HillComp >= 0 )
    {    
      if ( r->HillCoeff > 0 )
      {
        // calculate reaction constant for activation
        x = pow( (double) X[r->HillComp]/sys.volume, r->HillCoeff );
        r->k = r->Hillk * x / ( x + pow( r->HillConst, r->HillCoeff ) );
      }
      else if( r->HillCoeff < 0 ) 
      {
        // calculate reaction constant for repression 
        x = pow( r->HillConst, -r->HillCoeff );
        r->k = r->Hillk * x / ( x + pow( (double) X[r->HillComp] / sys.volume, -r->HillCoeff ) );
      }
    }

    // remember: maximum number of reactants is two by definition!
    if( 0 == r->Nreact ) 
    {
      // production without reactants (Volume dependent!)
      a[ids[i]] = r->k * sys.volume;

    }
    else if( 1 == r->Nreact ) 
    {
      // exactly one reactant
      a[ids[i]] = r->k * X[r->react[0].index];

    }
    else if( r->react[0].index == r->react[1].index ) 
    {
      // both reactants are identical
      a[ids[i]] = r->k / sys.volume * 
                  X[r->react[0].index] * (X[r->react[1].index] - 1);
	
    }
    else 
    {
      // both reactants are different
      a[ids[i]] = r->k / sys.volume * 
                  X[r->react[0].index] * X[r->react[1].index];

    }						
  }
}


/**
  Selects next reaction randomly.
  This is done based on the propensity functions, which have to be calculated beforehand.
  @param[out]  j      ID of the selected reaction
  @param[in]   sum_a  Total of the propensity functions
*/
void select_reaction( int *j, double sum_a )
{
  double rs, cumu_a;

  // pick some random real number out of [0,sum_a)
  rs = ran_get() * sum_a;
  
  // find related reaction
  *j = 0;
  cumu_a = a[*j];
  while( cumu_a <= rs ) 
  {
    ++(*j);
    cumu_a += a[*j];
  }

}


/**
  Updates the concentrations of the components after reaction j has happend.
  The function handles the effects of a firing of reaction `j`
  @param[in]  j   ID of the reaction, that takes place
*/
void update_concentrations( int j )				//should be: update_copy_numbers. this function is NOT called after growth step
{
  int i;

  // lower concentration of reactants
  for( i=0; i<R[j].Nreact; ++i )
    --X[R[j].react[i].index];

  // increase concentration of products
  for( i=0; i<R[j].Nprod; ++i )
    X[R[j].prod[i].index] += R[j].prod[i].change;

  // increase reaction counter
  ++react_count[j];

  // update the propensity functions, which have changed
  update_propensity_functions( react_network[j].val, react_network[j].len );		
}

void update_concentrations_queuereactants( int j )
{
  int i;

  // lower concentration of reactants
  for( i=0; i<R[j].Nreact; ++i )
    --X[R[j].react[i].index];

  update_propensity_functions( react_network[j].val, react_network[j].len );
}

void update_concentrations_queueproducts( int j )
{
  int i;

  // increase concentration of products
  for( i=0; i<R[j].Nprod; ++i )
    X[R[j].prod[i].index] += R[j].prod[i].change;

  // increase reaction counter
  ++react_count[j];

  // update the propensity functions which have changed
  update_propensity_functions( react_network[j].val, react_network[j].len );
}


/**
  Calculates the current volume accoring to the time.
  The volume is defined such that it doubles within the doubling time sys.doubling_time.
  Additionally, it is normalized such that the mean volume over time is one.
*/
void calculate_volume()
{
  if( GROWTH_LINEAR == sys.growth_type )
  {
    sys.volume = volume_min + growth_const * (sys.tau - sys.last_division);
  }
  else if( GROWTH_EXPONENTIAL == sys.growth_type ) 
  {
    sys.volume = volume_min * exp( growth_const*(sys.tau - sys.last_division) );
  }
}


/**
  Initializes a new cell cycle.
*/
void growth_init()
{
  double doubling_time; //< tentative doubling time for the current growth process.
  
  // draw a doubling time from a Gaussian distribution
  do
    doubling_time = sys.doubling_time + ran_gaussian( sys.doubling_time_std );
  while ( sys.tau >= sys.last_division + doubling_time ); 
  
  // draw duplication time from Gaussian distribution.
  // Make sure duplication time lies inside the cell cycle time.
  do {
    ccycle_duptime = doubling_time * sys.duplication_phase + 
      ran_gaussian( doubling_time * sys.duplication_phase_std ); }
  while ( sys.tau > sys.last_division + ccycle_duptime || ccycle_duptime >= doubling_time ); 
  ccycle_duptime += sys.last_division;
  
  // precalculation for the volume
  if( GROWTH_LINEAR == sys.growth_type )
  {
    volume_min = 2./3.;
    growth_const = volume_min / doubling_time;
    
    growth_dt = .001 * doubling_time;
  }
  else if( GROWTH_EXPONENTIAL == sys.growth_type )
  {
    /*
      The inition volume is set to Ln(2), such that the average volume over a period
      <V(t)> is equal to 1.    
    */
    growth_const = LOG2 / doubling_time;
    volume_min = LOG2;
    
    growth_dt = .001 * doubling_time;
  }
  else
  {
    printf( "Unknown growth type `%d`", sys.growth_type );
    abort();
  }
  
  // calculate the current volume
  calculate_volume();
}


/**
  Adjustes the volume and divides the bacterium if neccessary.
*/
void growth_step()
{
  unsigned int i;
  
  calculate_volume();

  if( sys.volume >= 2. * volume_min )
  {
    // store time of cell division.
    sys.last_division = sys.tau;

    //initialize the next growth period
    growth_init();
  
    // distribute chemical species between cells.
    for( i=0; i<sys.Ncomp; ++i )
    {
      if( !Xconst[i] )
      {
        if( DIVIDE_BINOMIALLY == sys.division_type )
          X[i] = ran_binomial( 0.5, X[i] );
        else if( DIVIDE_HALF == sys.division_type )
          X[i] /= 2;
        else
        {
          printf( "Unknown division type `%d`", sys.division_type );
          abort();
        }
      }
    }
  
    // disregard half the items on the queue
    if( sys.needs_queue )
      queue_disregard_randomly( 0.5 );

    // Half all prefactors of Hill functions with 'CanDuplicate' on
    for( i=0; i<sys.Nreact; ++i )
    {
      if( R[i].HillComp >= 0 ) 
      {
        if ( R[i].CanDuplicate )
        {
          //Only half gene copynbr, when gene duplication is not at cell division.
          if( ccycle_duptime > sys.last_division ) 
          {   
            R[i].Hillk = R[i].Hillk0 * sys.init_gene_copynbr;
            sys.current_gene_copynbr = sys.init_gene_copynbr;
          
            //Register duplication as data.      
            duplications.time[duplications.cntr] = sys.tau;
            duplications.number[duplications.cntr] = sys.init_gene_copynbr;
            duplications.cntr++;
          }
          else
          {
            duplications.time[duplications.cntr] = sys.tau;
            duplications.number[duplications.cntr] = sys.init_gene_copynbr;
            duplications.time[duplications.cntr+1] = sys.tau;
            duplications.number[duplications.cntr+1] = 2*sys.init_gene_copynbr;
            duplications.cntr += 2;        
          }          
        }
      }
    }

	  //Update all propensities (All numbers have changed).
	  determine_propensity_functions();
  }
  else{ 
    // check if genes duplicate.
	  if( sys.tau >= ccycle_duptime && sys.current_gene_copynbr < 2 * sys.init_gene_copynbr )
	  {
	    for( i=0; i<sys.Nreact; ++i )
      {
        if( R[i].HillComp >= 0 )
        {
          if ( R[i].CanDuplicate )
          {
            R[i].Hillk = R[i].Hillk0 * (++sys.current_gene_copynbr);
            duplications.time[duplications.cntr] = sys.tau;
            duplications.number[duplications.cntr] = sys.current_gene_copynbr;
            duplications.cntr++;
          }  
        }
      }
	  }

	  //update propensities which change under growth.
    update_propensity_functions( gdr, gdrl );
  }
}

void growth_dependent_reactions()          
{
	int i, j;
	
	j=0 ;
	
	for( i=0; i<sys.Nreact; ++i )
	{
		if( 2 == R[i].Nreact || R[i].HillComp >= 0 || 0 == R[i].Nreact )
		{
			gdr[j] = i ;
			j++ ;
//			printf("reaction %d dependent on volume \n", i);
		}	
		
	}
	gdrl = j ;
//a	printf("length of the array of growth dependent reactions: %d\n", gdrl);
	
	return ;
}




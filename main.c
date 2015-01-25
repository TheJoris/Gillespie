/**  
  @file main.c
  Main code of the program.
*/

#include <getopt.h>
#include <limits.h>

#include "Gillespie.h"
#include "propagate.h"
#include "queue.h"
#include "stats.h"
#include "random.h"
#include "read_comp_react.h"




/*------------------------Globally defined variables------------------------*/  

Sys      sys;            ///< general information on the system
int      *X;             ///< array for the current number of the components
long long int *react_count;   ///< array, which counts how often each reaction has been fired
char     **Xname;        ///< array for the names of the components
boolean  *Xconst;        ///< array stating if the respective component is constant over time 
int      **Xcalc;        ///< array for the components from which the value is calculated
int      *Xcalc_count;   ///< size of the previous array
double   *a;             ///< array for the propensity function of each reaction channel
int	 *gdr ;		           ///< array of reactions whose propensity functions chance when the volume chances
React    *R;             ///< array of structs used to store the reaction channels
IntArray *react_network; ///< array containting structural information about the reaction network
Eventpair duplications;  ///< Store all duplication events.

/*------------------------Locally defined functions--------------------------*/ 

int  start();
void allocate_memory();
void print_reactions( boolean show_count );
void print_initial_conditions();
int  finish();
void get_reaction_network();
void show_help();


/**
  Main program.
  @return error code or EXIT_SUCCESS if successful
*/
int main(int argc, char *argv[])
{
  int i;

  sys.input = NULL;
  sys.output = NULL;
  sys.needs_queue = FALSE;
  sys.volume = 1.;              // initial volume factor
  sys.doubling_time = -1;       // no growth at all
  sys.doubling_time_std = 0;    // the doubling time is exact
  sys.duplication_period_std = 0.; // Std deviation in gene doubling.
  sys.init_gene_copynbr = 1; // Number of gene doublings before doubling whole genome.
  sys.last_division = NAN;      // to distinguish from given division times
  sys.growth_type = GROWTH_EXPONENTIAL;
  sys.division_type = DIVIDE_HALF;
  sys.output_stats = FALSE;
  sys.output_conc = OUTPUT_AUTOMATIC;
  sys.output_kaic = FALSE;
  sys.output_phos = FALSE;
  sys.output_quite = FALSE;
  sys.species = NULL;
  sys.tau_init = 0.;
  
  // check command line parameters
  while( TRUE )
  {
    static struct option long_options[] =
    {
      // These options set a flag.
      {"a_bionmial",    no_argument, 0,  'a'},
      {"conc",          no_argument, 0,  'c'},
      {"exponential",   no_argument, 0,  'e'},
      {"linear",        no_argument, 0,  'l'},
      {"help",          no_argument, 0,  'h'},
      {"input",   required_argument, 0,  'i'},
      {"number",        no_argument, 0,  'n'},
      {"output",  required_argument, 0,  'o'},
      {"phos",          no_argument, 0,  'p'},
      {"quite",         no_argument, 0,  'q'},
      {"seed",    optional_argument, 0,  's'},
      {"total",         no_argument, 0,  't'},
      {"verbose",       no_argument, 0,  'v'},
      {"species", required_argument, 0, 1000},
    };

    // getopt_long stores the option index here.
    int option_index = 0;
    
    // parse the next option
    i = getopt_long( argc, argv, "acehlnpqstv", long_options, &option_index );
    
    // Detect the end of the options.
    if( i == -1 )
      break;
    
    switch(i) 
    {
      case 'a': // also --a_binomial
        sys.division_type = DIVIDE_HALF;
      break;
      case 'c': // also --conc
        sys.output_conc = OUTPUT_CONCENTRATION;
      break;
      case 'e': // also --exponential
        sys.growth_type = GROWTH_EXPONENTIAL;
      break;
      case 'h': // also --help
        show_help();
        return EXIT_SUCCESS;
      break;
      case 'l': // also --linear
        sys.growth_type = GROWTH_LINEAR;
      break;
      case 'i': // only --input
        if( NULL == optarg )
        {
          fprintf( stderr, "If the option '--input' is given, a filename has to be added" );
          return EXIT_FAILURE;
        }
        sys.input = malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
        strcpy( sys.input, optarg );
      break;
      case 'n': // also --number
        sys.output_conc = OUTPUT_COPYNUMBER;
      break;
      case 'o': // only --output
        if( NULL == optarg )
        {
          fprintf( stderr, "If the option '--output' is given, a filename has to be added" );
          return EXIT_FAILURE;
        }
        sys.output = malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
        strcpy( sys.output, optarg );
        printf( "files: %s - %s\n", sys.output, optarg );
      break;
      case 'p':
        sys.output_phos = TRUE;
      break;
      case 'q':
        sys.output_quite = TRUE;
      break;
      case 's':
        if( NULL == optarg )
        {
          // init random number generator with seed based on current time
          ran_init(0);
        }
        else
        {
          long seed;
          if( 1 != sscanf( optarg, "%ld%*s", &seed ) )
          {
            fprintf( stderr, "If a value for '--seed' is given, it has to be a long integer number.\n" );
            return EXIT_FAILURE;
          }
          ran_init( seed );
        }
      break;
      case 't':
        sys.output_kaic = TRUE;
      break;
      case 'v':
        sys.output_stats = TRUE;
      break;
      case 1000: // only --species
        if( NULL == optarg )
        {
          fprintf( stderr, "If the option '--species' is given, an identifier has to be added" );
          return EXIT_FAILURE;
        }
        sys.species = malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
        strcpy( sys.species, optarg );
      break;
      default:
        abort();
    }
  }

	// load data from files and initialize program
  start();
	
	// cleanup
  return finish ();
}


/**
  Load data from files and initialize program.
*/
int start()
{
  int  i,Ncomp,Nqueue;
  double tau_equi, tau_prod;
  long steps_equi, steps_prod;
  char *dummy;
  dummy = (char*) malloc(40);
  FILE *fp;

  // load general information from Gillespie.inp
  if( (fp = fopen("Gillespie.inp","r")) == NULL ) 
  {
    printf("Cannot open Gillespie.inp.\n");
    abort();
  }
  
  // read the information
  sys.name = calloc( MAXNAMELENGTH, sizeof(char) );
  fscanf( fp, "%s%*s", sys.name ); 
  fscanf( fp, "%d%*s", &sys.Ncomp );
  fscanf( fp, "%d%*s", &sys.Nreact );
  fscanf( fp, "%ld\t%ld\t%ld%*s", &steps_equi, &steps_prod, &sys.stats_steps );
  fscanf( fp, "%lg%*s", &sys.dt );
  fscanf( fp, "%lg\t%lg%*s", &tau_equi, &tau_prod );
  fscanf( fp, "%lg\t%lg%*s", &sys.doubling_time, &sys.doubling_time_std );
  fscanf( fp, "%lg\t%d%*s", &sys.duplication_period_std, &sys.init_gene_copynbr );
  
  // convert the number of blocks, which the calculation should run, into steps
  // or set it to the maximum number, when it should be neglected
  if( steps_equi < 0 )
    steps_equi = LONG_MAX;
  else
    steps_equi *= sys.stats_steps;

  if( steps_prod < 0 )
    steps_prod = LONG_MAX;
  else
    steps_prod *= sys.stats_steps;
  
  // ignore negative times
  if( tau_equi < 0 )
    tau_equi = INF_POS;
  if( tau_prod < 0 )
    tau_prod = INF_POS;
    
  // close file
  fclose(fp);

  printf("===============================================================================\n");
  printf("This program propagates the chemical master equation according\n");
  printf("to the Gillespie-algorithm.\n");
/*  printf("Log book information.\n\n");
  if (!(log_init(sys.name,TRUE)==0))
    printf("No log book information could be obtained.\n");*/
  printf("-------------------------------------------------------------------------------\n");
  printf("System parameters.\n\n");
  printf("Name of the run                   %8s\n",sys.name);
  printf("Number of components              %8d\n",sys.Ncomp);
  printf("Number of reaction channels       %8d\n",sys.Nreact);
/*  printf("Number of equilibrium blocks      %8d\n",*n_blk_eq);
  printf("Number of production  blocks      %8d\n",*n_blk_run);
  printf("Number of steps per block         %8d\n",*n_steps);*/
  printf("Timesteps of writeout             %8f\n",sys.dt);
  printf("Total time(steps) of equilibrium run     %8f\n",tau_equi);
  printf("Total time(steps) of production run      %8f\n",tau_prod);
  if( HAS_GROWTH )
  {
    printf( "Doubling time                     %8f +- %8f\n", sys.doubling_time, sys.doubling_time_std );
    printf( "Duplication phase                  %8f +- %8f, Number of gene duplications: %d\n",  sys.doubling_time, sys.duplication_period_std, sys.init_gene_copynbr );   
  }

  // use this information to allocate memory and read in components and reactions
  allocate_memory();
  read_components();
  read_reactions();
  get_reaction_network();

  //Set duplication events array.
  duplications.len = (int) 2 * (tau_prod + tau_equi)  / sys.doubling_time;
  duplications.cntr = 0;
  duplications.time = malloc( (sys.init_gene_copynbr + 1) * duplications.len * sizeof( double ) );
  duplications.number = malloc( (sys.init_gene_copynbr + 1) * duplications.len * sizeof( double ) );
  
  // check, if the state of the system has to be loaded from a file
  if( sys.input != NULL )
  {
    printf( "Load initial state of the system from '%s'...\n", sys.input );
    
    if( (fp = fopen(sys.input, "r")) == NULL) 
    {
      printf("The file '%s' could not be opened\n", sys.input);
      return EXIT_FAILURE;
    }
    
    fscanf( fp, "%lg\t%d\t%d\t%lg\t%*s", 
            &sys.tau_init, &Ncomp, &Nqueue, &sys.last_division 
          );
    if( Ncomp != sys.Ncomp )
    {
      printf("The number of components do not match (%d in the definition and %d in the state)\n",sys.Ncomp,Ncomp);
      return EXIT_FAILURE;
    }
    
    for( i=0; i<Ncomp; i++ )
    {
      fscanf( fp, "%d\t%s", &X[i], dummy );
    }
    
    queue_read( fp, Nqueue, sys.tau_init );
    
    fclose( fp );
  }
  else
  {
    sys.tau_init = 0; // initialize values, which are loaded from file otherwise
  }
  
  // print the chemical system that will be investigated
  if( sys.output_stats )
    print_reactions( FALSE );
  print_initial_conditions();
  
  if( HAS_GROWTH )
  {
	printf("the system has GROWTH \n");
    if( sys.growth_type == GROWTH_LINEAR )
    {
      fprintf( stderr, "Growth is linear\n" );
      if( isnan( sys.last_division ) )
      {
        // set last_division to current time.
        sys.last_division = sys.tau_init;
      }
    }
    else if ( sys.growth_type == GROWTH_EXPONENTIAL )
    {
      fprintf( stderr, "Growth is exponential\n" );
      if( isnan( sys.last_division ) )
      {
        // set last_division to current time.
        sys.last_division = sys.tau_init;// - 0.528766*sys.doubling_time;
        // 0.528766 == -log(log(2))/log(2)
      }
    }
    else
    {
      printf( "Unknown growth type `%d`", sys.growth_type );
      abort();
    }
  }
  
  if( sys.output_conc == OUTPUT_AUTOMATIC )
  {
    if( HAS_GROWTH )
      sys.output_conc = OUTPUT_CONCENTRATION;
    else 
      sys.output_conc = OUTPUT_COPYNUMBER;
  }
  
  if( sys.output_conc == OUTPUT_CONCENTRATION )
  {
    fprintf( stderr, "Output data are concentrations (mean volume = 1).\n" );
  }
  else
  {
    fprintf( stderr, "Output data is copynumbers.\n" );
  }
  
  printf( "===============================================================================\n" );
  
  // Fire up the Gillespie algorithm.
  printf( "Start run...\n" );
  run( tau_prod, tau_equi, steps_prod );
 
  return EXIT_SUCCESS;
}


/**
  Allocates memory for all variables
*/
void allocate_memory ()
{
  int i;

  // variables used to describe the reactions
  R = calloc( sys.Nreact, sizeof(React) );
  react_count = calloc( sys.Nreact, sizeof(long long int) );
  a = malloc( sys.Nreact*sizeof(double) );
  gdr = malloc( sys.Nreact*sizeof(int) );
  for( i=0; i<sys.Nreact; i++ )
  {
    a[i] = 0.;
	gdr[i] = 0.;
  }
  
  // variables used to describe the components
  X           = calloc( sys.Ncomp, sizeof(int) );
  Xname       = calloc( sys.Ncomp, sizeof(char *) );
  Xconst      = calloc( sys.Ncomp, sizeof(boolean) );
  Xcalc       = calloc( sys.Ncomp, sizeof(int *) );
  Xcalc_count = calloc( sys.Ncomp, sizeof(int) );

  for( i=0; i<sys.Ncomp; ++i ) 
    Xname[i] = calloc( MAXNAMELENGTH, sizeof(char) );
   
  // output stats system
  if( sys.output_stats )
  {
    Xblk  = calloc( sys.Ncomp, sizeof(Stats) );
    Xrun  = calloc( sys.Ncomp, sizeof(Stats) );
  }
  
}


/**
  Print all reaction channels on screen.
  @param[in]  show_count  Flag determining, if the number of times a reaction has fired should be plotted
*/
void print_reactions( boolean show_count )
{
  int i,j;
  
  printf( "\nThe following reactions are simulated:\n\n" );
  for( i=0; i<sys.Nreact; i++ ) 
  {
    if( show_count )
      printf( "%3d (Count: %11lld): ", i, react_count[i] );
    else
      printf( "%3d: ", i );
    
    // print reactants
    if( R[i].Nreact == 0 ) 
      printf(" Null ");
    else 
    {
      for( j=0; j<R[i].Nreact; j++ ) 
      {
        if( j > 0 )
          printf( "+" );
        if( R[i].react[j].change == 1 )
          printf( " %s ", Xname[R[i].react[j].index] );
        else
          printf( " %2d %s ", R[i].react[j].change, Xname[R[i].react[j].index] );
      }
    }
    
    // print type of reaction
    if( REAC_DELAYED == R[i].type )
      printf( "\t==>" );
    else
      printf( "\t-->" );
    
    // print products
    if( R[i].Nprod == 0 ) 
      printf( " Null " );
    else
    {
      for( j=0; j<R[i].Nprod; j++ )
      {
        if( j > 0 )
          printf( "+" );
        if( R[i].prod[j].change == 1 )
          printf( " %s ", Xname[R[i].prod[j].index] );
        else
          printf( " %2d %s ", R[i].prod[j].change, Xname[R[i].prod[j].index] );
      }
    }
    
    // print conditions
    if( 0 == R[i].HillCoeff )
    {
      printf( "\tk = %4.3f", R[i].k );
    }
    else
    {
      printf( "\tk = %4.3f * [%s]^n / ( K^n+[%s]^n ), K=%4.3f, n=%.1f, CanDuplicate=%1d",
              R[i].Hillk, Xname[R[i].HillComp], Xname[R[i].HillComp], R[i].HillConst, 
                R[i].HillCoeff, R[i].CanDuplicate );
    }
    
    if( R[i].time > 0 || R[i].sigma > 0 )
    {
      printf( "\tT = %4.3f\ts = %4.3f ", R[i].time, R[i].sigma );
    }
/*    printf("\n\tinfluencing ");
    for( j=0; j<react_network[i].len; j++ )
      printf( " %2d", react_network[i].val[j] );*/
    printf("\n");
  }
}


/**
  Prints the initial conditions of the raection system.
*/
void print_initial_conditions()
{
  int i;
  
  printf( "\nThe initial conditions are:\n\n" );
  for( i=0; i<sys.Ncomp; i++ ) 
  {
    if( X[i] != 0 )
      printf( "[%s] \t= %d\n", Xname[i], X[i] );
  }
  if( queue_length() > 0 )
  {
    printf( "The length of the queue is %d.\n", queue_length() );
  }
  printf( "The origin of the timeline is %f.\n", sys.tau_init );
  printf( "\n" );
}


/**
  Finish program and output final statistics to screen.
*/
int finish ()
{
  FILE *fp;
  int i;
  
  // write reaction count to stdout
  print_reactions( TRUE );
  
  // output, that program has finished
  printf( "\nFinal time: %f", sys.tau );
  printf( "\nFinal volume: %f", sys.volume );
  printf("\n\n===============================================================================\n");
  printf("Run completed.\n");
//   log_exit();
  
  // check, if the state of the system has to be written to a file
  if( sys.output != NULL )
  {
    printf( "Write final state of the system to '%s'...\n", sys.output );
    
    if( (fp = fopen(sys.output, "w")) == NULL) 
    {
      printf("The file '%s' could not be created\n", sys.output);
      return EXIT_FAILURE;
    }
    
    fprintf( fp, "%f\t%d\t%d\t%f\tTime_Ncomp_QueueLength_LastDivision\n", 
             sys.tau, sys.Ncomp, queue_length(), sys.last_division 
           );
    
    for( i=0; i<sys.Ncomp; i++ )
      fprintf( fp, "%d\t%s\n", X[i], Xname[i] );
    
    if( queue_length() > 0 )
    {
      fprintf( fp, "Qindex_Qtime\n" );
      //queue_shift(sys.tau-sys.tau_init); //correct for previous shift
      queue_write( fp );
    }
    
    fclose( fp );
  }
  return EXIT_SUCCESS;
}

/**
  Determines the reaction network of the system.
  The reaction network describes which reaction directly
  influence other reaction and especially their propensity functions.
  For delayed reactions, it's important, that both after decreasing 
  the reactants and also after increasing the products, the propensity
  functions may change, whereas in reactions of type A --> A + B without
  any delay the component A may be disregarded, since its concentration
  does not change.
*/
void get_reaction_network()
{
  int i, j, m, n_r, *c, *r;
    
  // allocate memory
  react_network = malloc( sys.Nreact * sizeof(IntArray) );
  c = malloc( sys.Ncomp * sizeof(int) );       // affected components
  r = malloc( (sys.Nreact+1) * sizeof(int) );  // affected reactions
  
  // run through all reactions and check which reactions they influence
  for( i=0; i<sys.Nreact; i++ )
  {
    // printf( "\nReaction %d influences components: ", i );
    
    // reset buffer
    for( m=0; m<sys.Ncomp; m++ ) // buffer for affected components
    {
      c[m] = 0;
    }
    for( j=0; j<sys.Nreact+1; j++ )  // buffer for affected reactions
      r[j] = -1;
    n_r = 0; // number of affected reactions
    
    // find all affected components of reaction i

    // look for reactants
    for( m=0; m<R[i].Nreact; m++ )
//       if( R[i].type == REAC_DELAYED )
//         c[R[i].react[m].index] = 1;  // ensure that c is positive!
//       else
        c[R[i].react[m].index]--;  // count right
    
    
    // look for products
    for( m=0; m<R[i].Nprod; m++ )
      c[R[i].prod[m].index] += R[i].prod[m].change;
    
    // print output
	  if(i == 1)
	  {
      for( m=0; m<sys.Ncomp; m++ )
        if( c[m] != 0 )
          printf( "%s (%d), ", Xname[m], c[m] );
      printf( "\nand reactions: " );
	  }
    
    // find all reaction which are influenced by those components
    for( j=0; j<sys.Nreact; j++ )
    {
      // check, if reactants are influenced
      for( m=0; m<R[j].Nreact; m++ )
      {
        // check whether reactant is calculated
        if( Xconst[m] == 2 )
        {
          r[n_r++] = j;
          goto next; // leave two loops at once
        }
        
        // check whether reaction is influenced by the current reaction i
        if( 0 != c[R[j].react[m].index] )
        {
          r[n_r++] = j;
          goto next; // leave two loops at once
        }
      }
      
      // check, if Hillfunction is influenced
      if( R[j].HillComp >= 0 && R[j].HillConst != 0 &&
          ( c[R[j].HillComp] != 0 || Xconst[R[j].HillComp] == 2 ) )
	    {
        r[n_r++] = j;
	    }
      
      next: continue; // check the next reaction
    }
    
    // save those reactions
    react_network[i].val = (int*) malloc( n_r * sizeof(int) );				
    react_network[i].len = n_r;								
    
    for( j=0; j<n_r; j++ )
    {
    react_network[i].val[j] = r[j];

    }
  }
  
  free(c);
  free(r);
}


/**
  Shows the help message for the program.
 */
void show_help() 
{
  printf(
      "Usage: Gillespie.exe [OPTION]\n"\
      "Program running the Gillespie-Algorithm. The model is defined by the compulsory file 'Gillespie.inp',\n"\
      "which must reside in the current working directory.\n"\
      "\n"\
      "Possible options:\n"\
      "\t-a  --a_bionmial   do cell division by exactly dividing in half\n"\
      "\t-c  --conc         output concentrations\n"\
      "\t-n  --number       output copynumbers\n"\
      "\t-e  --exponential  set the growth mechanism to exponential growth\n"\
      "\t-l  --linear       set the growth mechanism to linear growth\n"\
      "\t-s, --seed[=nr]    initializes the used random number generator with the integer `nr`\n"\
      "\t                   If `nr` is not given or is zero a random seed is used\n"\
      "\t-q, --quite        Prevent the output of single files for each species\n"\
      "\t-p, --phos         Outputs the phosphorylation ratio for the Kai system\n"\
      "\t-t, --total        Outputs the total amount of KaiC for the Kai system\n"\
      "\t    --species=id   Outputs only the chemical species which matches `id`\n"\
      "\n"\
      "\t    --input=file   define input file, where the initial state of the system is given\n"\
      "\t    --output=file  define output file, where the final state of the system is written to\n"\
      "\t-h, --help         display this help and exit\n"\
      "\t-v, --verbose      outputs intermediate stats while calculating\n"
  );
}

/** 
  @file Gillespie.h
  Header file for Gillespie simulation.
  Definition of usefull constants, macros and most structures.
 */

#ifndef _GIL_H_
#define _GIL_H_

/*---------------------------INCLUDE'S---------------------------------------*/
 
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------Some uesful definitions------------------------------*/


#define EQUIL       0
#define RUN         0

#define MAXPROD         10 ///< maximum number of products in one reaction
#define MAXNAMELENGTH   30

// types of reactions
#define REAC_NORMAL    0
#define REAC_DELAYED   1

// defines the growth-type
#define GROWTH_LINEAR      0
#define GROWTH_EXPONENTIAL 1

// output types
#define OUTPUT_AUTOMATIC      0
#define OUTPUT_COPYNUMBER     1
#define OUTPUT_CONCENTRATION  2

// division types
#define DIVIDE_HALF         0
#define DIVIDE_BINOMIALLY  1

/**
  Fake boolean type
 */
typedef int boolean;
#define TRUE        1
#define FALSE       0
#ifndef NAN
  #define NAN         0.0/0.0  ///< definition for not-a-number
#endif

#define HAS_GROWTH  ( sys.doubling_time > 0 )

#define CANCELLED   -1
#define F_COMP(a,b) (fabs ( a-b ) < 1e-10)
#define MIN(a,b)    ( a < b ? a : b )
#define MAX(a,b)    ( a > b ? a : b )
#define EPSILON     1e-9
#define PI          3.141592653589793
#define INF_POS     1.0/0.0

/*-----------------------Structure definitions-------------------------------*/

/**
  Structure containing statistical information for a component.
*/
typedef struct stats_type 
{
  int     acc;    ///< Number of steps 
  double  sum,    ///< integral of the component count over time
          sumsq,  ///< integral of the component count squared over time
          err,    ///< err := &lt;X^2&gt; - &lt;X&gt;^2
          noise;  ///< noise := err / &lt;X&gt;
} Stats; 

/**
  Structure containing general information on the system.
*/
typedef struct sys_type 
{
  int     Ncomp,         ///< Number of components
          Nreact,        ///< Number of reactions
          needs_queue,   ///< Switch determining, if the simulation needs a queue
          growth_type,   ///< Equation to use for determining the growth type
          division_type, ///< Determines the method used for division
          init_gene_copynbr, ///< number of gene copies in the cell at the start of the cell cycle.
          current_gene_copynbr; ///> number of gene copies at the current phase of the cell cycle.
  long    stats_steps,   ///< Number of steps calculation inbetween each output
          equi_steps,    ///< Number of steps for the first run
          prod_steps;    ///< Number of steps for the second run
  double  volume,        ///< Volume factor
          doubling_time, ///< Doubling time of the bacterium
          doubling_time_std, ///< Standard deviation of the doubling time
          last_division, ///< Time of the last division
          duplication_period_std, ///< Standard deviation of the duplication time.
          tau,           ///< Current time
          tau_init,      ///< Time of previous runs, read from input-file
          dt;            ///< Time step of data write out
  char    *name,         ///< Name of the run
          *input,        ///< possible file for input
          *output,       ///< possible file for output
          *species;      ///< possible restriction onto the species, which should be output
  boolean output_stats,  ///< Switch determining, if stats should be output
          output_conc,   ///< Switch determining, if the output should be concentrations instead of copynumbers
          output_kaic,   ///< output total amount of KaiC instead of all components
          output_phos,   ///< output phosphorylation ratio instead of all components
          output_quite;  ///< switch determining whether all species should be written to files
} Sys;

/**
  Information on one component for a reaction.
*/
typedef struct stoch_type 
{
  int  index,  ///< Index of the component in the component list
       change; ///< Number of the components, which are created by one reaction
} Stoch;

/**
  Structure describing a reaction channel.
*/
typedef struct reaction_type 
{

  // general information
  int    Nreact,    ///< Number of reactants 
         Nprod,     ///< Number of products
         type;      ///< type of reaction
  double k,         ///< reaction constant
  // delayed reactions
         time,      ///< average delay of the reaction 
         sigma;     ///< sigma of the Gaussian distribution used for specifing the delay
  // reactions of Hill type
  int    HillComp,  ///< Index of the component used in the Hill function
         CanDuplicate;///< Bollean if the gene is duplicated during cell cycle.
  double Hillk,     ///< Prefactor of the Hill function; k = k0 * active_genes.
         Hillk0,    ///< Prefactor for Hillfunction for a single gene.
         HillCoeff, ///< The Hill coefficient, n
         HillConst; ///< The concentration of the HillComp for half maximum effect, K

        
  Stoch  react[2],      ///< Array for describing the reactants
         prod[MAXPROD]; ///< Array for describing the products

} React;

/**
  Structure describing an array of integers
*/
typedef struct int_array 
{
  int len,  ///< The length of the array
      *val; ///< A pointer to the first element
} IntArray;

typedef struct eventpair_type 
{
  double *time; ///< Store event times
  int *number,  ///< of when the number changes
  len,          ///< so many times.
  cntr;         /// Counter of current position in array.
} Eventpair;

/*-----------------------Externally declared variables----------------------*/

extern Sys      sys;          ///< general information on the system
extern int      *X;           ///< array for the current number of the components
extern long long int *react_count; ///< array which counts how often each reaction has been fired
extern char     **Xname;      ///< array for the names of the components
extern boolean  *Xconst;      ///< array stating if the respective component is constant over time 
extern int      **Xcalc;      ///< array for the components from which the value is calculated
extern int      *Xcalc_count; ///< size of the previous array
extern double   *a;           ///< array for the propensity function of each reaction channel
extern int	*gdr;	      ///array of reactions whose propensity functions chance when the volume chances
extern Eventpair duplications; ///storing all the gene duplication events.
extern React    *R;           ///< array of structs used to store the reaction channels
extern IntArray *react_network; ///< array containting structural information about the reaction network

#endif

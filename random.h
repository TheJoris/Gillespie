// Header file for the random number generator

// default seed, if random number generator is not initialized
#define MSEED 161803398

// constant for a fake boolean value
#define FALSE 0
#define TRUE  1

// functions
unsigned long int random_seed();
void ran_init( long seed );
double ran_get();
double ran_gaussian( const double sigma );
int ran_binomial( double pp, int n );
double gammln( double xx );

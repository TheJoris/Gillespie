#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#define Pi 3.141592654

#define BLOCK 100
#define BREAK 10000000

#define PHOSP	0
#define TOTAL	1

#define NO	0
#define YES	1

void mem_alloc();
void mem_realloc();
void read_input();
void mem_alloc2();	
void determine_peak();
void determine_first_peak();	
void determine_next_peak();
void print_output();
double get_kappa(double t_peak, int &i);
			

int input_type;
int calc_force;
char *inputfile;
char *duplicationfile;
double Td;					//doubling time
double threshold;   //threshold for proteine production in hill function.
double Tk;					//estimate of the clock time  --always around 25h, right?
double first_div;	//time of first division --has to be entered
double num_bins;		//number of bins of the histograms --set equal to T of the 1st period for the time being 

double upper;			  //bounds for the search window for the next peak
double lower;
double dup;					//bounds for reporting warnings about possibly false peaks
double dlow;


double *x;					//array of data
double *p;					//array of peaks times
double *ph;				//array of peak heights
double *l;					//array of period lengths 
int   *T;					//histogram of period lengths
double *patp;			//array of phases at peaks   
int   *H;					//histogram of peak phases
						//array of length num_bins
						//at output this has to be converted to phases

double *duplication_times;	//array of duplication_times.
		
double mean_ph ;					//mean peak height
double mean_patp ;				//mean phase at peak --which I assume to equal the perfectly locked state if the run is long enough
						

double time_init ;
double time_end ;
double time_step ;
int array_mem ;
int array_id ;	

int prev_peak ;					//id of the previous peak
int d_count ;
int p_count ; 
int delta ;					//number of steps between peaks -- 				
						//this number is used to increase the efficiency of the algorithm by 
						//'guessing' the position of the next peak once one is found


//Warning: beware of first peaks right at the start of the run, it's hard for the algorithm to determine these

int main(int argc, char *argv[])
{

	//int j ;

	d_count = 0 ;
	p_count = 0 ;
	mean_ph = 0. ;
	mean_patp = 0. ;
	array_id = 0 ;
	Td = 0 ;
	Tk = 27. ;
	threshold = 1.;
	upper = 1.25 ;
	lower = 0.75 ;
	num_bins = Tk ;	
	time_step = 0 ;
	time_init = 0 ;
	time_end = 0 ;
	first_div = 0 ;
	inputfile = NULL ;
	input_type = PHOSP ;
	calc_force = NO ;


	int gt;

	int option_index = 0;

	static const struct option long_options[] = {
   	{ "duplicationfile", required_argument, NULL, 'b' },
  	{ "threshold", required_argument, NULL, 'a' },
		{ "doubling time", required_argument, NULL, 'd' },
		{ "upper boundary", required_argument, NULL, 'u' },
		{ "lower boundary", required_argument, NULL, 'l' },
		{ "first division", required_argument, NULL, 'f' },
		{ "clock period", required_argument, NULL, 'k' },
		{ "time step", required_argument, NULL, 's' },
		{ "inputfile", required_argument, NULL, 'i' },
		{ "total", no_argument, NULL, 't' },
		{ "phosph. ratio", no_argument, NULL, 'p' },
		{NULL, no_argument,NULL, 0}
	};

	while(gt != -1)
	{
	
	  gt = getopt_long(argc, argv, "a:b:d:f:i:k:tpu:l:s:r", long_options, &option_index);

		switch(gt)
		{
		  case 'a' :	if( NULL == optarg )
						{
							printf("An argument is required for the threshold\n" );
							abort() ;
						}
						
						if( 1 != sscanf( optarg, "%lf%*s", &threshold ) )		
						{
							printf("warning: no threshold given\n");
							return EXIT_FAILURE;
						}
						printf("threshold is %f \n", threshold);
						break ;
		  
		  case 'b' :	if( NULL == optarg )
						{
							printf("A name is required for the duplication file\n" );
							abort() ;
						}
						duplicationfile = (char*) malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
						strcpy( duplicationfile, optarg );
						printf( "duplication file: %s\n", duplicationfile );
						break ;
		
		
			case 'i' :	if( NULL == optarg )
						{
							printf("A name is required for the inputfile\n" );
							abort() ;
						}
						inputfile = (char*) malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
						strcpy( inputfile, optarg );
						printf( "proteine files %s\n", inputfile );
						break ;
						
			case 'd':	if( NULL == optarg )				
						{
							printf("An argument is required for the doubling time\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &Td ) )		
						{
							printf("error, doubling time not a valid number: %f\n",Td);
							return EXIT_FAILURE;
						}
						printf("doubling time is %f \n", Td);
						break ;

			case 'f':	if( NULL == optarg )				
						{
							printf("An argument is required for the first division time\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &first_div ) )		
						{
							printf("error, first division time not a valid number\n");
							return EXIT_FAILURE;
						}
						printf("first division time is %f \n", first_div );
						break ;

			case 'k':	if( NULL == optarg )				
						{
							printf("An argument is required for the clock period\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &Tk ) )		
						{
							printf("error, clock period not a valid number: %f\n",Tk);
							return EXIT_FAILURE;
						}
						printf("estimated clock period is %f \n", Tk);
						break ;

			case 'u':	if( NULL == optarg )				
						{
							printf("An argument is required for the upper boundary\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &upper ) )		
						{
							printf("error, upper boundary not a valid number \n");
							return EXIT_FAILURE;
						}
						printf("upper boundary for the search window for the next peak is %f \n", upper);
						break ;

			case 'l':	if( NULL == optarg )				
						{
							printf("An argument is required for the lower boundary \n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &lower ) )		
						{
							printf("error, lower boundary not a valid number \n");
							return EXIT_FAILURE;
						}
						printf("lower boundary for the search window for the next peak is %f \n", lower);
						break ;

			case 's':	if( NULL == optarg )				
						{
							printf("An argument is required for the time step\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%lf%*s", &time_step ) )		
						{
							printf("error, time_step is not a valid number\n");
							return EXIT_FAILURE;
						}
						printf("time_step is %f \n", time_step);
						break ;

      			case 't':	input_type = TOTAL;
					break;

      			case 'p':	input_type = PHOSP;
					break;

      			case 'r':	calc_force = YES;
					break;
			
			case -1 : break;
			    
			default  :	printf("no valid choice for switch \n") ;
						abort() ;
		}
	}
	
	mem_alloc();

	read_input();

  //(Too) Very important parameter for peak searching!
	delta = (int) (Tk / time_step) ;
	printf("delta is %d \n", delta);	

	mem_alloc2();

  prev_peak = 0;
	determine_peak();
	do
	{
		p_count++;
		determine_peak();
	}
	while( p[p_count] < (time_end-Tk) );
			
  printf("Number of peaks found: %d\n",p_count);

	mean_patp = mean_patp / (p_count + 1 ) ;

	print_output() ;

	free(x) ;
	free(p) ;
	free(ph) ;
	free(T) ;
	free(H) ;
	free(l) ;
	free(patp) ;

	return 0;
}


void mem_alloc()
{
	int i ;
	x = (double*) malloc( BLOCK * sizeof(double) );

	array_mem = BLOCK;

	for(i=0; i<array_mem; i++)     
	{
		x[i] = 0.;
  }
  
	return;
}


void mem_alloc2()
{
	int i,
      num_cycl, 
      num_dupl;   

	num_cycl = (int) (2 * (time_end - time_init) / Tk) ;

	printf("estimated number of proteine cycles: %d \n", num_cycl/2) ;

	p = (double*) malloc( num_cycl * sizeof(double) );
	ph = (double*) malloc( num_cycl * sizeof(double) );
	l = (double*) malloc( num_cycl * sizeof(double) );
	patp = (double*) malloc( num_cycl * sizeof(double) );
	
	for( i=0; i<num_cycl; i++ )    
	{
   	p[i] = 0.;
		ph[i] = 0.;
		l[i] = 0.;
		patp[i] = 0.;
  }

	return;
}


void mem_realloc()
{
	int i ;

  // allocate memory for the data structure
	array_mem += BLOCK;
	x = (double*) realloc( x, array_mem * sizeof(double) );

  // initialize values
	for( i=array_mem-BLOCK; i<array_mem; i++)
	{
    x[i] = 0.;
 	}
 	
	return ;
}


void read_input()										
{
	double dummy;
	int dummy2,cellcycle_nbr;

	FILE *fp;

  //Read proteine concentration timetrace.

	if( (fp = fopen(inputfile,"r")) == NULL )
	{
		printf("cannot open inputfile \n");
		abort();
	}

	if(input_type == TOTAL)
	{
		fscanf( fp, "#time		Total KaiC	step \n");
	}
	else if(input_type == PHOSP)
	{
		fscanf( fp, "#time		Phos. Ratio	step \n");
	}
	
	char *mufassa;
  mufassa = (char*) calloc( 30, sizeof(char) );
	fscanf( fp, "\n%[^\n]", mufassa );

	array_id = 0;

	fscanf( fp, "%le\t%le\t%d\n", &dummy, &x[array_id] , &dummy2 );
	array_id++;	
	
	time_init = dummy;

	fscanf( fp, "%le\t%le\t%d\n", &dummy, &x[array_id], &dummy2 );
	array_id++;
	
	//printf("2### %f - %d - %f\n",dummy,array_id,x[array_id-1]);

	time_step = dummy - time_init;

	printf("time_init is %f \n", time_init);
	printf("time_step is %f \n", time_step);

	while (!feof(fp) && array_id <BREAK)   
	{
		  fscanf( fp, "%le\t%le\t%d \n", &dummy, &x[array_id], &dummy2);
		  //printf("%f - %f\n",dummy,x[array_id]);
		  array_id++;							

		  if( array_id >= array_mem ) 
		    mem_realloc();
	}
	time_end = dummy; 

	fclose(fp);

  
  if(Td > 0)
  {
	  cellcycle_nbr = (time_end - time_init) / Td + 1;
	  printf("Time_end is %f, with %d data points. Est. cell cycle#: %d\n", time_end, array_id,     cellcycle_nbr-1);  

    //Process duplication times.
 	  duplication_times = (double*) malloc( 2 * cellcycle_nbr * sizeof(double) );
  
    if( (fp = fopen(duplicationfile,"r")) == NULL )
	  {
		  printf("Cannot open duplication file: %s \n",duplicationfile);
		  abort();
	  }

    //ignore first line.
	  fscanf( fp, "%s%*s\n", mufassa );

    array_id = 0;
    while (!feof(fp) && array_id <BREAK)   
	  {
	    fscanf( fp, "%le\t%d\n", &duplication_times[array_id], &dummy2);
	 	  //printf("Id %d, dup. time: %f\n", array_id, duplication_times[array_id]);
		  array_id++;
	  }
	
	  printf("Read %d gene number modifications events.\n", array_id);

	  fclose(fp); 
	}
	else {
	  printf("No Cell division. Time_end is %f, with %d data points.\n", time_end, array_id);
	}

	return;
}


void determine_peak()
{
	double dx1,dx2;
	int i, iend, istart;

  if(p_count == 0)
    istart = 0;
  else
  	istart = (int) prev_peak + lower*delta;
  	
	iend   = (int) prev_peak + upper*delta;
	
	for( i = istart; i<= iend; i++ )
	{
    dx1 = x[i]-x[i-1];
    dx2 = x[i+1]-x[i];
    
    //When sign change in derivatives: peak or valley.
    if(dx1 * dx2 < 0 && dx1 > dx2)
    {
       break;
    }
	}
	
	//extra check if this is indeed a peak and not a annomaly due to noise.
	dx1 = x[i-1]-x[i-2];
  dx2 = x[i+2]-x[i+1];
	if( (dx1 < 0 || dx2 > 0) && x[i] > threshold )
	{
	  printf("WARNING: noisy peak. Around i: %f - %f - %f \n",x[i-1],x[i],x[i+1]);
	  printf("WARNING: i - 1, i + 1. %f - %f // %f - %f \n",x[i-2],x[i-1],x[i+1],x[i+2]);
	}
	else
	{
  	p[p_count] = time_init + i * time_step;				
    ph[p_count] = x[i];
    mean_ph += x[i];
  
    if(p_count > 0)
      l[(p_count-1)] = p[p_count] - p[(p_count - 1)];	
	
	  prev_peak = i;
	
	  //printf("%d - %d - %d --- %dth peak at t = %f \n",istart,iend,i, p_count, p[p_count]) ;		

	  if(i == istart || i == iend)  
	  {
	  	printf("WARNING: peak at search boundary. %d at %f\n", p_count, p[p_count]);
	  }
  
  }

	return;
}


double get_kappa(double t_peak, int &i){

  int dupl_nbr = (int) 2 *((time_end - time_init) * Td + 10);
  
  while(t_peak > duplication_times[i] && i < dupl_nbr) i += 2;
  
  if(i < 2)
    return -1;
  else   
    return t_peak - duplication_times[i-2];
}


void print_output()
{
	int i,j;
	double phase;
	double time;

	FILE** fp = (FILE**) malloc(sizeof(FILE*) * 2);

	if( (fp[0] = fopen("oscillations.dat","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	mean_ph = mean_ph / (p_count + 1 );						
	fprintf(fp[0], "#number of peaks is %d \n", p_count);
	fprintf(fp[0], "#mean height of peaks is %f \n", mean_ph); 

  j=0;
  if(Td > 0){
  	fprintf(fp[0], "#time_peak - period from this peak - height peak - kappa \n"); 	
  	for(i=0; i<p_count; i++)    
     	fprintf( fp[0], "%f\t%f\t%f\t%f \n", p[i], l[i], ph[i], get_kappa(p[i],j) ); 
     	
  }
  else
  {
  	fprintf(fp[0], "#time_peak - period from this peak - height peak \n"); 	
    for(i=0; i<p_count; i++)    
     	fprintf( fp[0], "%f\t%f\t%f\n", p[i], l[i], ph[i] ); 
  }

	fclose(fp[0]);

	return;
}

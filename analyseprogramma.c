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
void filter_peaks(double delta);
void moving_average(int avr_len);
double get_kappa(double t_peak, int &i);
			

int input_type;
int avr_len;

char *inputfile;
char *duplicationfile;

double Td;					//doubling time
double threshold;   //threshold for proteine production in hill function.
double Tk;					//estimate of the clock time  --always around 25h, right?
double first_div;	  //time of first division --has to be entered

double upper;			  //bounds for the search window for the next peak
double lower;
double dup;					//bounds for reporting warnings about possibly false peaks
double dlow;
double filter_interval;

double *x;					//array of data
double *movavr_x;		//array of data, using a moving average of avr_len
double *p;					//array of tentative peak times
double *ph;				  //array of tentative peak heights
double *l;					//array of period lengths 
double *duplication_times;	//array of duplication_times.
		
double mean_ph ;					//mean peak height

double time_init;
double time_end;
double time_step;
int array_mem;
int array_id;	

int prev_peak;			//id of the previous peak
int maxima_cntr;		// Counter for the number of period maxima.
int p_count; 				// peak counter.
int p_count_max;		//maximum number of peaks.	

int main(int argc, char *argv[])
{
	p_count = 0 ;
	mean_ph = 0. ;
	array_id = 0 ;
	avr_len = 10;
	Td = 0 ;
	Tk = 27. ;
	threshold = 1.;
	upper = 1.25 ;
	lower = 0.75 ;
	time_step = 0 ;
	time_init = 0 ;
	time_end = 0 ;
	first_div = 0 ;
	inputfile = NULL ;
	input_type = PHOSP ;

	int gt;

	int option_index = 0;

	static const struct option long_options[] = {
   	{ "duplicationfile", required_argument, NULL, 'b' },
  	{ "threshold", required_argument, NULL, 'a' },
  	{ "average length", required_argument, NULL, 'c' },
  	{ "filter interval", required_argument, NULL, 'e' },
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
	
	  gt = getopt_long(argc, argv, "a:b:c:d:e:f:i:k:tpu:l:s:r", long_options, &option_index);

		switch(gt)
		{
			case 'c' :	if( NULL == optarg )
						{
							printf("An argument is required for the moving average +/-length\n" );
							abort() ;
						}
						
						if( 1 != sscanf( optarg, "%d%*s", &avr_len ) )		
						{
							printf("warning: no moving average length given\n");
							return EXIT_FAILURE;
						}
						printf("moving average length is %d \n", avr_len);
						break ;

			case 'e' :	if( NULL == optarg )
						{
							printf("An argument is required for the filter interval\n" );
							abort() ;
						}
						
						if( 1 != sscanf( optarg, "%lf%*s", &filter_interval ) )		
						{
							printf("warning: no filter interval given\n");
							return EXIT_FAILURE;
						}
						printf("Filter interval is %f \n", filter_interval);
						break ;

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
						printf( "proteine number timetrace file %s\n", inputfile );
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
		
			case -1 : break;
			    
			default  :	printf("no valid choice for switch \n") ;
						abort() ;
		}
	}
	
	mem_alloc();

	printf("Reading from input file: %s...\n",inputfile);
	read_input();

  //Setting parameters.
	p_count_max = (int) (10 * (time_end - time_init) / Tk);
	printf("Max. expected # of peaks: %d \n", p_count_max) ;

	mem_alloc2();

	//Calculate the moving average of signal, and use this to find peaks.
	moving_average(avr_len);

	p_count=-1;
  prev_peak = 2;
	printf("Starting peak search within (%f,%f).\n",time_init,time_end-Tk);
	do
	{
		p_count++;
		determine_peak();

		if(p_count >= p_count_max)
		{
			printf("Max. # of allowed peaks reached.\n");
			break;
		}
	} while( p[p_count] < (time_end-Tk) );
  printf("Number of tentative peaks found: %d\n", p_count + 1 );

	//Only keep the highest peaks within the interval [Delta/2,Delta/2].
	filter_peaks( filter_interval * Tk );
  printf("Number of period maxima found: %d\n", maxima_cntr );

	printf("Printing output to file...");
	print_output();
	printf("done.\n");

	free(x);
	free(p);
	free(ph);
	free(l);
	free(duplication_times);
	//free(movavr_x); Freeing this causes havoc.
	return 0;
}


void mem_alloc()
{
	int i;
	x = (double*) malloc( BLOCK * sizeof(double) );
	movavr_x = (double*) malloc( BLOCK * sizeof(double) );

	array_mem = BLOCK;

	for(i=0; i<array_mem; i++)     
	{
		x[i] = 0.;
		movavr_x[i] = 0.;
  }
  
	return;
}


void mem_realloc()
{
	int i;

  // allocate memory for the data structure
	array_mem += BLOCK;
	x = (double*) realloc( x, array_mem * sizeof(double) );
	movavr_x = (double*) realloc( x, array_mem * sizeof(double) );

  // initialize values
	for( i=array_mem-BLOCK; i<array_mem; i++)
	{
    x[i] = 0.;
    movavr_x[i] = 0.;
 	}
 	
	return ;
}


void mem_alloc2()
{
	int i;

	p  = (double*) malloc( p_count_max * sizeof(double) );
	ph = (double*) malloc( p_count_max * sizeof(double) );
	l  = (double*) malloc( p_count_max * sizeof(double) );
	
	for( i=0; i<p_count_max; i++ )    
	{
   	p[i] = 0.;
		ph[i] = 0.;
		l[i] = 0.;
  }

	return;
}


void read_input()										
{
	double dummy;
	int dummy2,cellcycle_nbr,i,j;

	FILE *fp;

  //Read proteine concentration timetrace.

	if( (fp = fopen(inputfile,"r")) == NULL )
	{
		printf("cannot open inputfile: %s \n",inputfile);
		abort();
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
	  cellcycle_nbr = (int) (time_end - time_init) / Td + 1;
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

    j = 0;
    while (!feof(fp) && j <BREAK)   
	  {
	    fscanf( fp, "%le\t%d\n", &duplication_times[j], &dummy2);
	 	  //printf("Id %d, dup. time: %f\n", array_id, duplication_times[array_id]);
		  j++;
	  }
	
	  printf("Read %d gene number modifications events.\n", j);

	  fclose(fp); 
	}
	else {
	  printf("No Cell division. Time_end is %f, with %d data points.\n", time_end, array_id);
	}

	return;
}


void moving_average(int avr_len)
{
	double avr_x;
	int i, j, iend, istart;

	printf("Moving average; %d data points, averaging length %d...",array_id,2*avr_len);
	
	for(i = avr_len; i < array_id - avr_len; i++)
	{
		istart = i - avr_len;
  	iend   = i + avr_len;
		avr_x = 0;
  	for(j = istart; j < iend; j++)
		{
			avr_x += x[j];
		}

		movavr_x[i - avr_len] = avr_x/(2*avr_len);
	}

	printf(" done\n");

	return;
}


void determine_peak()
{
	double dx1,dx2;
	int i;

	for( i = prev_peak + 1; i < array_id - 2 * avr_len; i++ )
	{
    dx1 = movavr_x[i]-movavr_x[i-1];
    dx2 = movavr_x[i+1]-movavr_x[i];
    
    //When sign change in derivatives, change is at top and above threshold.
    if( dx1 * dx2 <= 0 && dx1 > dx2 && movavr_x[i] > threshold)
	  	break;
	}
	
	p[p_count] = time_init + (i + avr_len) * time_step;				
  ph[p_count] = movavr_x[i];

  if(p_count > 0)
    l[(p_count-1)] = p[p_count] - p[(p_count - 1)];	
		
	//printf("%d - %d - %d --- %dth peak at t = %f \n",istart,iend,i, p_count, p[p_count]) ;		
  
  prev_peak = i;

	return;
}


void filter_peaks(double delta)
{
	int i,j,local_max_i;
	double local_max;

	//printf("Filtering out non-maxima peaks...");
	
	maxima_cntr=0;
	for(i = 0; i < p_count-1; i++)
	{
		//Find te first peak that is seperated within
		// a distance delta from the next.
		//printf("%f - %f\n",l[i],delta);

		local_max	= ph[i];
		local_max_i = i;
		if(l[i] < delta)
		{
			j = i;		
			do {
				j++;
				if(ph[j] > local_max)
				{
					local_max = ph[j];
					local_max_i = j;
				}
			} while(l[j] < delta && j < p_count - 1);

			i = j;
		}

		p[maxima_cntr] = p[local_max_i];
		ph[maxima_cntr] = local_max;
		if(maxima_cntr > 0)
	    l[(maxima_cntr-1)] = p[maxima_cntr] - p[(maxima_cntr - 1)];	

  	maxima_cntr++;
	}

	//For the last tent. maximum, since it has no next peak to ther right,
  //we assume it is maximum if it is more than a distance delta away 
	//from the left peak.
	if(l[p_count - 1] > delta || ph[p_count] > local_max)
	{
		local_max = ph[i];

		p[maxima_cntr] = p[i];
		ph[maxima_cntr] = local_max;
		if(maxima_cntr > 0)
	    l[(maxima_cntr-1)] = p[maxima_cntr] - p[(maxima_cntr - 1)];	

  	maxima_cntr++;
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

	FILE *fp;

	if( (fp = fopen("oscillations.dat","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	mean_ph=0;
	for(i=0;i<maxima_cntr;i++)
	{
		mean_ph += ph[i];
	}
	mean_ph /= maxima_cntr;
							
	fprintf(fp, "#number of peaks is %d \n", maxima_cntr);
	fprintf(fp, "#mean height of peaks is %f \n", mean_ph); 

  i=0,j=0;
  if(Td > 0){
  	fprintf(fp, "#time_peak - period from this peak - height peak - kappa \n"); 	
  	for(i=0; i<maxima_cntr; i++)    
     	fprintf( fp, "%f\t%f\t%f\t%f \n", p[i], l[i], ph[i], get_kappa(p[i],j) ); 
     	
  }
  else
  {
  	fprintf(fp, "#time_peak - period from this peak - height peak \n"); 	
    for(i=0; i<maxima_cntr; i++)    
     	fprintf( fp, "%f\t%f\t%f\n", p[i], l[i], ph[i] ); 
  }

	fclose(fp);

  //Write data of moving average of timetrace.
	if( (fp = fopen("moving_average.dat","w")) == NULL )
	{
		printf("cannot open file to write moving average. \n");
		abort();
	}

 	fprintf(fp, "time - proteine concentration \n"); 	
  for(i = 0; i < array_id - 2 * avr_len; i++)    
    fprintf( fp, "%f\t%f\n", (avr_len + i)*time_step, movavr_x[i] ); 
    
	fclose(fp);

	return;
}


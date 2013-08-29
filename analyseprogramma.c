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

void mem_alloc() ;
void mem_realloc() ;
void read_input() ;
void mem_alloc2() ;	
void determine_first_peak() ;	
void determine_peak() ;
void period_length()	;
void moment_of_division() ;
void restoring_force() ;
void print_output() ;
					

int input_type ;
int calc_force ;
char *inputfile ;
float Td ;					//doubling time
float Tk ;					//estimate of the clock time  --always around 25h, right?
float first_div ;		//time of first division --has to be entered
float num_bins ;		//number of bins of the histograms --set equal to T of the 1st period for the time being 

float upper ;			  //bounds for the search window for the next peak
float lower ;
float dup ;					//bounds for reporting warnings about possibly false peaks
float dlow ;


float *x ;					//array of data
float *p ;					//array of peaks times
float *ph ;					//array of peak heights
float *l ;					//array of period lengths 
int   *T ;					//histogram of period lengths
float *patp ;				//array of phases at peaks   
int   *H ;					//histogram of peak phases
						//array of length num_bins
						//at output this has to be converted to phases

float *dbeta ;					//array of phase deviations from locked phase
		
float mean_ph ;					//mean peak height
float mean_patp ;				//mean phase at peak --which I assume to equal the perfectly locked state if the run is long enough
						

float time_init ;
float time_end ;
float time_step ;
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

	int j ;

	d_count = 0 ;
	p_count = 0 ;
	mean_ph = 0. ;
	mean_patp = 0. ;
	array_id = 0 ;
	Td = 0 ;
	Tk = 27. ;
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
		{ "doubling time", required_argument, NULL, 'd' },
		{ "upper boundary", required_argument, NULL, 'u' },
		{ "lower boundary", required_argument, NULL, 'l' },
		{ "first division", required_argument, NULL, 'f' },
		{ "clock period", required_argument, NULL, 'k' },
		{ "number of bins", required_argument, NULL, 'n' },
		{ "time step", required_argument, NULL, 's' },
		{ "inputfile", required_argument, NULL, 'i' },
		{ "total", no_argument, NULL, 't' },
		{ "phosph. ratio", no_argument, NULL, 'p' },
		{ "restoring force", no_argument, NULL, 'r' },
		{NULL, no_argument,NULL, 0}
	};

	gt = getopt_long(argc, argv, "d:f:i:n:k:tpu:l:s:r", long_options, &option_index);

	while(gt != -1)
	{

		switch(gt)
		{
			case 'i' :	if( NULL == optarg )
						{
							printf("An argument is required for the inputfile\n" );
							abort() ;
						}
						inputfile = (char*) malloc( ( strlen(optarg) + 1 ) * sizeof( char ) );
						strcpy( inputfile, optarg );
						printf( "files: %s - %s\n", inputfile, optarg );
						break ;

			case 'd':	if( NULL == optarg )				
						{
							printf("An argument is required for the doubling time\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%f%*s", &Td ) )		
						{
							printf("error, doubling time not a valid number\n");
							return EXIT_FAILURE;
						}
						printf("doubling time is %f \n", Td);
						break ;

			case 'f':	if( NULL == optarg )				
						{
							printf("An argument is required for the first division time\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%f%*s", &first_div ) )		
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

						if( 1 != sscanf( optarg, "%f%*s", &Tk ) )		
						{
							printf("error, clock period not a valid number\n");
							return EXIT_FAILURE;
						}
						printf("estimated clock period is %f \n", Tk);
						break ;

			case 'u':	if( NULL == optarg )				
						{
							printf("An argument is required for the upper boundary\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%f%*s", &upper ) )		
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

						if( 1 != sscanf( optarg, "%f%*s", &lower ) )		
						{
							printf("error, lower boundary not a valid number \n");
							return EXIT_FAILURE;
						}
						printf("lower boundary for the search window for the next peak is %f \n", lower);
						break ;

			case 'n':	if( NULL == optarg )				
						{
							printf("An argument is required for the number of bins\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%f%*s", &num_bins ) )		
						{
							printf("error, number of bins is not a valid number\n");
							return EXIT_FAILURE;
						}
						printf("number of bins is %f \n", num_bins);
						break ;

			case 's':	if( NULL == optarg )				
						{
							printf("An argument is required for the time step\n" );
							abort() ;
						}

						if( 1 != sscanf( optarg, "%f%*s", &time_step ) )		
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


			default  :	printf("no valid choice for switch \n") ;
						abort() ;
		}

		gt = getopt_long(argc, argv, "d:f:i:n:k:tpu:l:s:r", long_options, &option_index);
	}

	mem_alloc();

	read_input();

	delta = (int) (Tk / time_step) ;

	dup = (1. + 2*(upper)) / 3. ;   
	dlow = (1. + 2*(lower)) / 3. ;

	printf("delta is %d \n", delta);	

	mem_alloc2() ;

	determine_first_peak();

	do
	{
		determine_peak() ;
		period_length()	;
		moment_of_division() ;
	}
	while( p[p_count] < (time_end-10*Tk) ) ;
			

	mean_patp = mean_patp / (p_count + 1 ) ;

	if(calc_force == YES) 
	  restoring_force();

	print_output() ;

	free(x) ;
	free(p) ;
	free(ph) ;
	free(T) ;
	free(H) ;
	free(l) ;
	free(patp) ;
	free(dbeta) ;

	return 0;
}


void mem_alloc()
{
	int i ;
	x = (float*) malloc( BLOCK * sizeof(float) );

	array_mem = BLOCK ;

	for(i=0; i<array_mem; i++)     
		x[i] = 0.;

	return ;
}


void mem_realloc()
{
	int i ;

  // allocate memory for the data structure
	array_mem += BLOCK;
	x = (float*) realloc( x, array_mem * sizeof(float) );
  
  // initialize values
	for( i=array_mem-BLOCK; i<array_mem; i++)
	{
    x[i] = 0.;
 	}
	return ;
}


void read_input()										
{
	float dummy;
	int dummy2;

	FILE *fp;

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
	
	array_id = 0;

	fscanf( fp, "%e\t%e\t%d \n", &dummy, &x[array_id] , &dummy2 );
	array_id++;	
	printf("1### %f - %d - %f\n",dummy,array_id,x[array_id-1]);
	
	time_init = dummy;

	fscanf( fp, "%e\t%e\t%d \n", &dummy, &x[array_id], &dummy2 );
	array_id++;
	
	printf("2### %f - %d - %f\n",dummy,array_id,x[array_id-1]);

	time_step = dummy - time_init;

	printf("time_init is %f \n", time_init) ;
	printf("time_step is %f \n", time_step) ;

	while (!feof(fp) && array_id <BREAK)   
	{
		  fscanf( fp, "%e\t%e\t%d \n", &dummy, &x[array_id], &dummy2);
		  array_id++;							

		  if( array_id >= array_mem ) 
		    mem_realloc();
	}

	fclose(fp);

	time_end = dummy ; 
	printf("time_end is %f, with %d data points\n", time_end, array_id) ;                            

	return;
}


void mem_alloc2()
{
	int i ;
	int num_cycl ;   

	num_cycl = (int) (1.5 *(time_end - time_init) / Tk) ;

	printf("estimated number of cycles: %d \n", num_cycl) ;

	p = (float*) malloc( num_cycl * sizeof(float) );
	ph = (float*) malloc( num_cycl * sizeof(float) );
	l = (float*) malloc( num_cycl * sizeof(float) );
	patp = (float*) malloc( num_cycl * sizeof(float) );
	dbeta = (float*) malloc( num_cycl * sizeof(float) );


	for( i=0; i<num_cycl; i++ )    
	{
   	p[i] = 0. ;
		ph[i] = 0. ;
		l[i] = 0. ;
		patp[i] = 0. ;
		dbeta[i] = 0. ;
  }

	T = (int*) malloc( ((int) (15*Tk)) * sizeof(int) );

	for( i=0; i<((int) (15*Tk)); i++ )    
	{
        	T[i] = 0. ;
  }

	H = (int*) malloc( ((int) (num_bins+1)) * sizeof(int) );	

	for( i=0; i<(num_bins+1); i++ )     
	{
        	H[i] = 0. ;
  }

	return ;
}


void determine_first_peak()
{
	
	int temp_id;
	float temp_h;
	int i;

	prev_peak = 0;
	
	temp_h = x[0];
	temp_id = 0.;
	
	//Find highest data point between [0,1/2 * delta)
	for(i=1; i<(0.5*delta); i++)			
	{
		if(x[i] > temp_h)
		{
			temp_h = x[i];
			temp_id = i;
		}		
	} 

  /** If highest point DOES NOT lie within [0.1 * delta, 0.4 * delta]
    Find the highest point within [0.5 * delta, delta].
  
  Rationale: If the highest point in the first half of the interval lies not in the middle of it,
  the highest point lies on the boundaries of the semi-interval. 
  Then, the peak is likely the lie in the second half.
  */
	if(temp_id > 0.4*delta || temp_id < 0.1*delta)
	{
		while(i < (delta+1))			
		{
			i++;

			if(x[i] > temp_h)
			{
				temp_h = x[i];
				temp_id = i;
			}		
		} 
	}
	
	/** If highest point lies within (.9 * delta, delta]
	      Search for highest point in (delta, 1.5 * delta );
	
	  If the highest point lies at the boundery, the peak might be 
	  outside the estimated interval.
  */
	if(temp_id > 0.9*delta)
	{
		while(i< (1.5*delta))			//maakt niet uit, hoeft niet perse een integer te zijn
		{
			i++ ;

			if(x[i] > temp_h)
			{
				temp_h = x[i] ;
				temp_id = i ;
			}		
		} 

	}
	
	p[0] = time_init + temp_id * time_step;
	prev_peak = temp_id;
	mean_ph = temp_h;
	ph[0] = temp_h;
	
	printf("first peak at t = %f \n", p[0]) ;
	
	return ;
}


void determine_peak()
{
	int temp_id;
	float temp_h;
	int i;

	i = 0;
	temp_h = 0;
	temp_id = 0;	
	p_count++;	

	for( i= (int) (prev_peak +  lower*delta); i<(prev_peak+ upper*delta); i++ )
	{
    //printf("i = %d \n", i) ;
		if(x[i] > temp_h)
		{
			temp_h = x[i] ;
			temp_id = i ;
		}
	}

	p[p_count] = time_init + temp_id*time_step ;				
	ph[p_count] = temp_h ;
	mean_ph += temp_h ;		
							
//	printf("peak at t = %f \n", p[p_count]) ;

//warning if the highest point is at the edges of the window -- we might not have found a peak in this case	

	if(temp_id > (prev_peak + dup*delta) || temp_id< (prev_peak + dlow*delta))  
	{
		printf("peak %d at %f might not be a real peak\n", p_count, p[p_count]);
	}				



//	delta = temp_id - prev_peak ;						//uitgeschakeld voor meer betrouwbaar algorithme
	prev_peak = temp_id ;
//	printf("delta = %d, prev_peak = %d \n", delta , prev_peak) ;
	

	return ;
}


void period_length()
{

	int i ;

//calculate the length of this period

	l[(p_count -1)] = p[p_count] - p[(p_count - 1)] ;				//arbitrary way of numbering

//	printf("the length of the period is = %f \n", l[(p_count-1)]) ;

//put the length in a histogram

	i=0;

	while(  (0.2*(i+0.5)) < l[(p_count-1)] )    i++ ;

//	printf("bin number for period length %d\n", i) ;
	if(i > 15*Tk)
	{
		printf("error: something goes wrong in the histogram algorithm for period length at time %f\n", p[p_count -1]) ;
		abort() ;	
	}	

	T[i]++ ;

	return ;
}


//check if there has been a division during this period	-- if so put the phase at peak into the histogram

void moment_of_division()			//naam veranderen!
{

	int i ;
	float temp_d ;
	temp_d = 0. ;

	while( (first_div + (d_count +1)*Td) < p[p_count - 1]) d_count ++ ;

	temp_d = p[p_count - 1] - (first_div + d_count*Td)  ; 

//	printf("temp_d is %f\n", temp_d) ;
	

	patp[p_count-1] = (temp_d / Td) * 2 * Pi ;			//dit moet anders!!!
	mean_patp += patp[p_count -1] ;

	i = 0 ;

	while( ((i+1) / num_bins) < (temp_d / Td) ) i++ ;		//should both be smaller than 1			
											
//	printf("bin number for moment of division is %d\n", i) ;	
	if(i > num_bins)
	{
		printf("error: something goes wrong in the histogram algorithm for moments of divison at time %f\n", p[p_count -1]) ;
		abort() ;	
	}
	H[i]++ ;
	

	return ;
}


void print_output()
{

	int i ;
	float phase ;
	float time ;

	FILE** fp = (FILE**) malloc(sizeof(FILE*) * 4);


	if( (fp[0] = fopen("phase@division.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	if( (fp[1] = fopen("Tkhistogram.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	if( (fp[2] = fopen("periodlengths.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	if( (fp[3] = fopen("peakheights.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	
	fprintf(fp[0], "#number of divisions is %d \n", d_count);
	fprintf(fp[0], "#mean phase at peak is is %f \n", mean_patp); 
	fprintf(fp[0], "#phase		number of division at this phase \n");
 
	fprintf(fp[1], "#length of period		number of periods with this length \n"); 
	fprintf(fp[2], "#time		length of period \n"); 
	

	mean_ph = mean_ph / (p_count + 1 ) ;						
	fprintf(fp[3], "#number of peaks is %d \n", p_count);
	fprintf(fp[3], "#mean height of peaks is %f \n", mean_ph); 
	fprintf(fp[3], "#time_peak	height peak \n"); 	

	for( i=0; i<num_bins; i++ )     
	{	
		phase = (i+0.5)* 2*Pi / num_bins ;
		fprintf( fp[0], "%f\t%d \n", phase, H[i] );
	}

	for( i=0; i<((int) (10*Tk)); i++ )    
	{
        	fprintf( fp[1], "%f\t%d \n", (0.2*i), T[i] );
  }

	for( i=0; i<p_count; i++ )    
	{
		time = (p[i] + p[(i)]) / 2  ;
        	fprintf( fp[2], "%f\t%f \n", time, l[i] ) ;
  }

	for( i=0; i<p_count; i++ )    
	{
        	fprintf( fp[3], "%f\t%f \n", p[i], ph[i] ) ;
  }

	fclose(fp[0]) ;
	fclose(fp[1]) ;
	fclose(fp[2]) ;
	fclose(fp[3]) ;

	return ;

}


void restoring_force()
{

	int i ;
	int j ;

	float bins ;
	float range ;
	float div ;
	
	float temp ;
	float temp2 ;
	

	float *db ;
	float *ddb_mean ; 
	float *bin_mean ;

	FILE *fq ;
	FILE *fr ;

/*
	bins = 16. ;
	range = 0.4 ;
	div = bins / (2*range) ;

	//allocate memory for *db histogram and corresponding mean ddb array
	db = malloc( ((int) (bins+1)) * sizeof(float) );
	ddb_mean = malloc( ((int) (bins+1)) * sizeof(float) );	
	bin_mean = malloc( ((int) (bins+1)) * sizeof(float) );	

	for( i=0; i<(bins+1); i++ )     
	{
        	db[i] = 0. ;
		ddb_mean[i] = 0. ;
		bin_mean[i] = 0. ;
    	}

*/
	for( i=0; i< p_count; i++ )    
	{	
		dbeta[i] = patp[i] - mean_patp ;
/*
		//neglect deviations outside of range
		if(dbeta[i-1] < range && dbeta[i-1] > -range )
		{
			j=0  ;

			while( ((j+1) / div) - range < dbeta[i-1] ) j++ ;

			if(j > bins)
			{
				printf("error: something goes wrong in the restoring force function, j =  %d , dbeta = %f \n", j , dbeta[i-1]) ;
				abort() ;	
			}
			
			db[j] ++ ;
			ddb_mean[j] += ddbeta[i-1] ;
			bin_mean[j] += dbeta[i-1] ;
		}
*/
    	}

/*
	if( (fq = fopen("p@darray.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	fprintf(fq, "#p@d - mean	d to next point		number of items in histogram \n"); 

	for( i=0; i<(bins+1); i++ )     
	{

		if(db[i] > 0)
		{
			bin_mean[i] = bin_mean[i] / db[i] ;
			ddb_mean[i] = ddb_mean[i] / db[i] ;
//			temp2 = ((i+0.5) / div) - range ;

        		fprintf( fq, "%f\t%f\t%f \n", bin_mean[i], ddb_mean[i], db[i] ) ;
		}

    	}

	fclose(fq) ;

	free(ddb_mean) ;
	free(db) ;
*/

	if( (fr = fopen("scatterplot.txt","w")) == NULL )
	{
		printf("cannot open outputfile \n");
		abort();
	}

	fprintf(fr, "#dbeta_i	dbeta_i+1 \n"); 


	for( i=1; i< p_count; i++ )    
	{
		fprintf( fr, "%f\t%f \n", dbeta[i], dbeta[i+1]) ;		
	}	
	
	fclose(fr) ;
}


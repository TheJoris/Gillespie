#include "Gillespie.h"
#include "read_comp_react.h"

/**
  Reads components from the file %name%.components
*/
void read_components()
{
  int  i,j,Ncomp;
  char filename[40];
  FILE *fp;

  // determine filename and open the file
  sprintf( filename, "%s.components", sys.name );
  if( ( fp = fopen( filename, "r" ) ) == NULL ) 
  {
    printf("could not open %s\n",filename);
    abort();
  }
  else 
  {
    // check the number of components
    fscanf( fp, "%d\n", &Ncomp );
    if( Ncomp != sys.Ncomp ) 
    {
      printf( "The number of components is %d\n", Ncomp );
      printf( "The number of components should be %d\n", sys.Ncomp );
      abort();
    }
    else
    {
      // read in all components
      for (i=0;i<sys.Ncomp;i++)	
      {
        fscanf( fp, "%d\t\t%s\t%d", &X[i], Xname[i], &Xconst[i] );
        // check whether the component must be calculated
        if ( Xconst[i] == 2 )
        {
          fscanf( fp, "\t%d", &Xcalc_count[i] );
          Xcalc[i] = calloc( Xcalc_count[i], sizeof( int ) );
          for (j=0;j<Xcalc_count[i];j++ )
            fscanf( fp, "\t%d", &Xcalc[i][j] );
        }
        fscanf( fp, "\n" );
      }
    }
    
  }
  return;
}

/**
  Reads reactions from the file %name%.reactions
*/
void read_reactions()
{
  int  i,j,Nreact;
  char filename[40],dummy[80];
  FILE *fp;

  // determine filename and open the file
  sprintf( filename, "%s.reactions", sys.name );
  if (( fp = fopen( filename, "r" ) ) == NULL ) 
  {
    printf( "could not open %s\n", filename );
    abort();
  }
  else 
  {
    // check the number of reactions
    fscanf( fp, "%d%*s\n", &Nreact );
    printf( "Nreact is %d\n", Nreact );
    if( Nreact != sys.Nreact ) 
    {
      printf( "The number of reactions is %d\n", Nreact );
      printf( "The number of reactions should be %d\n", sys.Nreact );
      abort();
    }
    else 
    {
      // read in all reactions
      for( i=0; i<sys.Nreact; i++ )
      {
        
        fscanf( fp, "%lg %d %d %lg %lg %d %lg %lg %d %s\n",
                &R[i].k, &R[i].Nreact, &R[i].Nprod, &R[i].time, &R[i].sigma,
                &R[i].HillComp, &R[i].HillConst, &R[i].HillCoeff, &R[i].CanDuplicate, &dummy
              );
        R[i].Hillk0 = R[i].k; // copy rate constant to Hill prefactor
        
        R[i].Hillk = R[i].Hillk0 * sys.init_gene_copynbr;
        
        // check, if the reaction constant is defined by an Hill function
        if( R[i].k == 0 || R[i].HillConst == 0 || R[i].HillComp < 0)
        {
          R[i].HillComp = -1;
          R[i].CanDuplicate = 0;
        }

        // scan reactants
	      if( R[i].Nreact==0 )
        {
	        fscanf( fp, "%s", &dummy );
        }
	      else
        {
	        fscanf( fp, "%s %d", &dummy, &R[i].react[0].index );
          R[i].react[0].change = 1;
        }
	      for( j=1; j<R[i].Nreact; j++ ) 
        {
	        fscanf( fp, "%s %s %d", &dummy, &dummy, &R[i].react[j].index );
          R[i].react[j].change = 1;
        }
	      
	      fscanf( fp, "%s", &dummy );
	      
        // scan products
        if( R[i].Nprod==0 )
        {
	        fscanf( fp, "%s", &dummy );
        }
	      else
        {
	        fscanf( fp, "%d %s %d\n", 
                  &R[i].prod[0].change, &dummy, &R[i].prod[0].index );
        }
	      for( j=1; j<R[i].Nprod; j++ ) 
        {
	        fscanf( fp, "%s %d %s %d\n", &dummy, &R[i].prod[j].change,
		              &dummy, &R[i].prod[j].index );
        }
        
        // check type of the equation
        if( R[i].time == 0 && R[i].sigma == 0 )
        {
          R[i].type = REAC_NORMAL;
        }
        else
        {
          R[i].type = REAC_DELAYED;
          sys.needs_queue = TRUE;
        }
      }
    }
  }
  return;
}

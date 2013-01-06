/*
*  C Implementation of a queue
*
* Description: 
*   General structure for modelling a queue.
*   The elements are kept unordered and only the information about the first item to come
*   off the queue is stored separetly.
*
* Author: David Zwicker <mail@david-zwicker.de>, (C) 2009
*
*/

#include "queue.h"
#include "random.h"

Queue *Q;            ///< array of elements on the queue
int Nqueue=0,        ///< number of items in the queue
    queuelength=0,   ///< amount of places allocated for the use as a queue
    Qnext=0;         ///< id of the next item in the queue
double Qtau=INF_POS; ///< time of the next item in the queue


/**
  Does initialization for the queue.
  This function needs only to be called, if the queue should be cleared.
  @warning  The function does not free memory
*/
void queue_init()
{
  // init state variables
  Nqueue = 0;
  Qtau = INF_POS;
  Qnext = 0;
}


/**
  Puts a new item on the queue with the index id.
  @param[in]  id     ID for the item on the queue
  @param[in]  time   the time after which the item gets of the queue
 */
void queue_push( const int id, const double time )
{
//	printf("de queue wordt gebruikt \n");
  // check, if there is enough memory, otherwise enlarge
  if( Nqueue >= queuelength )
  {
    queuelength += QUEUELENGTH;
    if( queuelength > QUEUEMAX )
    {
      fprintf( stderr, "The queue exceeded the maximum length of %d", QUEUEMAX );
      abort();
    }
    Q = realloc( Q, queuelength * sizeof( Queue ) );
  }
  
  // set new item
  Q[Nqueue].time = time;
  Q[Nqueue].index = id;
  
  // check, if it will be the first item of the queue
  if( Q[Nqueue].time < Qtau )
  {
    Qtau = Q[Nqueue].time;
    Qnext = Nqueue;
  }
  
  Nqueue++;
  
//   printf("push: %d, %d, %f\n", rid, Nqueue, Q[Nqueue-1].time );

}


/**
  Get an item of the queue
  @param[in]  pos  index in the array of items of the queue
 */
int queue_pop( const int pos )
{
  int i, id;
  
  // get index of the item
  id = Q[pos].index;
  
  // if the removed item has not been the last one in the queue, 
  // put the last item of the queue in the place of the removed item
  // to ensure a compact array of elements
  if( pos < Nqueue-1 )
  {
    Q[pos].index = Q[Nqueue-1].index;
    Q[pos].time = Q[Nqueue-1].time;
//     memcpy( &Q[qid], &Q[Nqueue-1], sizeof(Queue) );
  }

  // shorten queue (although no memory is freed)
  Nqueue--;

  // determine next item in queue
  Qtau = INF_POS;
  for( i=0; i<Nqueue; i++ )
  {
    if( Q[i].time < Qtau )
    {
      Qtau = Q[i].time;
      Qnext = i;
    }
  }

  return id;
}


/**
  Checks the queue for an item which has finished before the time `time`.
  @param[in,out]  time   The time before which the item has to come off. If an item could be found its due time is returned in this value.
  @return   The ID, which has been stored on the queue
*/
int queue_check( double *time )
{
  // check, if an event has happend, or propagate to next event, when requested
  if( Qtau <= *time ) // an item gets off the queue
  {
    // set time to this point
    *time = Qtau;
    
    // pop the item
    return queue_pop( Qnext );
  }
  return -1; // no item in the period
}


/**
  Advances the `time` to the next item on the queue and returns its ID.
  If there are no more items on the queue, -1 is returned.
  @param[out] time  The time at which the item comes off the queue. If there are no items, this value is unaltered.
  @return The ID, which has been stored on the queue.
*/
int queue_advance( double *time )
{
  if( Nqueue > 0 )
  {
    *time = Qtau;
    return queue_pop( Qnext );
  }
  else
  {
    return -1;
  }
}


/**
  Get the length of the current queue.
  @return The length of the queue, e.g. the number of elements on it.
*/
int queue_length()
{
  return Nqueue;
}


/**
  Shifts the internal time.
  @param[in]  dt  The amount by which the time is shifted.
*/
void queue_shift( const double dt )
{
  unsigned int i;
  
  for( i=0; i<Nqueue; i++ )
    Q[i].time += dt;
  
  Qtau += dt;
}


/**
  Randomly drops half the items of the queue.
  @param[in]  prob  Survival probability (= fraction that will be kept in the queue)
*/
void queue_disregard_randomly( double prob )
{
  unsigned int i, n=0;
  
  Qtau = INF_POS;
  
  for( i=0; i<Nqueue; i++ )
  {
    if( ran_get() < prob ) // keep item if random number is below `prob`
    {
      Q[i-n] = Q[i];
      // check, if it will be the first item of the queue
      if( Q[i].time < Qtau )
      {
        Qtau = Q[i].time;
        Qnext = i-n;
      }
    }
    else // drop item otherwise
    {
      n++;
    }
  }
  
  Nqueue -= n;
}


/**
  Writes the data of the queue to a file.
  @param[in]  file  The handle of the (already opened) file
*/
void queue_write(FILE *file)
{
  int i;
  
  for( i=0; i<Nqueue; i++ )
  {
    fprintf( file, "%d\t%f\n", Q[i].index, Q[i].time );
  }
}

/**
  Reads the data of the queue from a file and adds it to the current queue.
  @param[in]  file          The handle of the (already opened) file
  @param[in]  num           The number of items which have to be read.
  @param[in]  time_origin   The current time, to validate input and disregard all items of the past
                            Set it to -infinity if you don't need this feature.
 */
void queue_read(FILE *file, int num, double time_origin)
{
  int i;
  
  // check, if there is enough memory, otherwise enlarge
  if( Nqueue + num > queuelength )
  {
    queuelength += num;
    Q = realloc( Q, queuelength * sizeof( Queue ) );
  }

  // read the items
  for( i=0; i<num; i++ )
  {
    fscanf( file, "%d\t%lg%*s", &Q[Nqueue].index, &Q[Nqueue].time );
    
    // check, if item comes off in the future => otherwise, overwrite item by not stepping counter
    if( Q[Nqueue].time > time_origin )
    {
    
      // check, if it will be the first item of the queue
      if( Q[Nqueue].time < Qtau )
      {
        Qtau = Q[Nqueue].time;
        Qnext = Nqueue;
      }
    
      Nqueue++;
    }
  }
}

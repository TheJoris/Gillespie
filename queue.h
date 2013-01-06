//
// C Interface: Implementation of a queue
//
// Description: 
// Implements a queue which can be used to store indices (and therefore also 
// pointers to other structures). The elements are stored with a time associated
// and there exist methods to retrieve the next element (that is with the smallest
// time) from the queue.
//
// Author: David Zwicker <mail@david-zwicker.de>, (C) 2009
//
//

#include <stdio.h>
#include <stdlib.h>

#define QUEUELENGTH  100      ///< initial length of the queue - will be extended by multiplies of this amount, if neccessary
#define QUEUEMAX     16*65536    ///< maximum length of the queue to prevent freezing of the program

#define INF_POS      1.0/0.0  ///< definition for infinity

/**
  Strucutre holding information of one item on the queue
 */
typedef struct queue_type 
{
  
  int index;    ///< ID to tell apart different items on the queue
  double time;  ///< absolute time, after which the queue is due
  
} Queue;

void queue_init();
void queue_push( const int rid, const double time );
int queue_check( double *time );
int queue_advance( double *time );
int queue_length();
void queue_shift( const double time );
void queue_disregard_randomly( double prob );
void queue_write(FILE *file);
void queue_read(FILE *file, int num, double time_origin);

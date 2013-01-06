/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _STATS_H_
#define _STATS_H_

extern Stats *Xblk, ///< Array for all components with statistical information in current block
             *Xrun; ///< Array for all components with statistical information in total run


/*-----------------------Externally declared functions----------------------*/

void block_init();
void block_acc();
void block_finish(int block);

#endif

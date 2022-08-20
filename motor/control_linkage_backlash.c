/* The code was generated automatically */
/* date: 3/5/2018		 time: 13:30:33 */

#include "control_linkage_backlash.h"
#include <math.h>


void control_linkage_backlash_init( control_linkage_backlash_DATA *data )
{
	data->freq = 10.000000; /* the frequence*/
	data->tick_num = 0; /* the tick number*/ 
	data->data0 = 0.000000;
	data->data1 = 0.000000;
	data->data2 = 0.000000;
	data->data3 = 0.000000;
	data->Vhod = 0.000000;
	data->Vyhod = 0.000000;
}

void control_linkage_backlash_tick( control_linkage_backlash_DATA *data )
{
 	float stack[3];
	stack[0] = data->Vhod;
	data->data0 = stack[0];
	stack[0] = data->data0;
	stack[1] = M_PI/120;
	stack[0] = stack[0]-stack[1];
	stack[1] = 1.000000;
	stack[0] = stack[1]*stack[0];
	data->data2 = stack[0];
	stack[0] = data->data2;
	stack[1] = data->data1;
	stack[0] = ( stack[0] >= stack[1])? 1 : 0;
	stack[1] = data->data2;
	stack[0] = stack[1]*stack[0];
	stack[1] = data->data0;
	stack[2] = -M_PI/120;
	stack[1] = stack[1]-stack[2];
	stack[2] = 1.000000;
	stack[1] = stack[2]*stack[1];
	data->data3 = stack[1];
	stack[1] = data->data3;
	stack[2] = data->data1;
	stack[1] = ( stack[1] <= stack[2])? 1 : 0;
	stack[2] = data->data3;
	stack[1] = stack[2]*stack[1];
	stack[0] = stack[1]+stack[0];
	data->data3 = stack[0];
	stack[0] = data->data3;
	stack[0] = fabs( stack[0] );
	stack[0] = (stack[0] > 0.0001 )? 1 : 0;
	stack[0] = ( stack[0] == 0 ) ? 1 : 0 ;
	stack[1] = data->data1;
	stack[0] = stack[1]*stack[0];
	stack[1] = data->data3;
	stack[0] = stack[1]+stack[0];
	data->data1 = stack[0];
	stack[0] = data->data1;
	data->Vyhod = stack[0];

	data->tick_num++;
}

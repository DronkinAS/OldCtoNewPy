/* The code was generated automatically */
/* date: 3/5/2018		 time: 13:30:33 */

#ifndef ____CONTROL_LINKAGE_BACKLASH____NAME_FILE___
#define ____CONTROL_LINKAGE_BACKLASH____NAME_FILE___

typedef struct struct_control_linkage_backlash_DATA
{
	float freq; /* the frequence*/
	float tick_num; /* the tick number*/ 
	float data0;
	float data1;
	float data2;
	float data3;
	float Vhod;/*Вход*/
	float Vyhod;/*Выход*/
} control_linkage_backlash_DATA;

void control_linkage_backlash_init( control_linkage_backlash_DATA *data );
void control_linkage_backlash_tick( control_linkage_backlash_DATA *data );

#endif
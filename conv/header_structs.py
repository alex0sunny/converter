'''
Created on 15 дек. 2017 г.

@author: Alexey
'''
from ctypes import *
import fields_names as fds

class TxtHeader(Structure):
    _fields_ = [
                 ("txt_header", c_char * 2),
                 ("serial", c_char * 10),
                ]

class BinHeader (Structure):
    _pack_ = 1
    _fields_ = [
                (fds.SIVY_SIZE, c_byte),              #13
                (fds.SIVY_CONFIG_WORD, c_ubyte),        #14 3/14 
                (fds.SIVY_CHANNEL_BIT_MAP, c_byte),     #15 4/15 
                (fds.SIVY_SAMPLE_BYTES, c_byte),        #16 5/16 
                
                (fds.SIVY_BLOCK_SAMPLES, c_uint),     #20 u32 BlockSamples;
                
                (fds.SIVY_SAMPLE_TIME, c_longlong),        #28 s64 SampleTime;
                (fds.SIVY_GPS_TIME, c_longlong),           #36 s64 GPSTime;   
                (fds.SIVY_DRIFT, c_longlong),             #44 s64 Drift;      

                (fds.SIVY_RSVD0, c_uint),             #48 u32 rsvd0;      

                #/* компас и акселерометр */
                (fds.SIVY_ACC_Z, c_short),             #50 s16 acc_z;  /* 11/50 milli-g, where g is gravitational acceleration*/
                (fds.SIVY_ACC_X, c_short),             #52 s16 acc_x;  /* 12/52 */
                (fds.SIVY_ACC_Y, c_short),             #54 s16 acc_y;  /* 13/54 */
                (fds.SIVY_COMP_Z, c_short),            #56 s16 comp_z; /* 14/56 milli-gauss */
                (fds.SIVY_COMP_X, c_short),            #58 s16 comp_x; /* 15/58 */
                (fds.SIVY_COMP_Y, c_short),            #60 s16 comp_y; /* 16/60 */

                #/* Параметры среды: температура, напряжения и пр */
                (fds.SIVY_PRESS_REG, c_uint),        #64 u32 press_reg;      
                (fds.SIVY_POWER_REG, c_ushort),        #66 u16 power_reg;      
                (fds.SIVY_TEMPER_REG, c_short),        #68 s16 temper_reg;     
                (fds.SIVY_HUMIDITY, c_short),          #70 s16 humidity;   // tenth of percent
                (fds.SIVY_H, c_short),                 #72 s16 h;  /* 20/70 height, meter */
                (fds.SIVY_LAT, c_int),               #76 s32 lat;        
                (fds.SIVY_LON, c_int),               #80 s32 lon;        
                ]

'''
Created on 14 дек. 2017 г.

@author: Alexey
'''
TABLE_SHOT_N = "shot_n"
TABLE_SAMPLE_RATE = "sample_rate"
TABLE_SAMPLE_TIME = "sample time"
TABLE_TUNED_DRIFT = "tuned_drift"
TABLE_TEMPER_DRIFT = "temper_drift"

SIVY_SIZE = "size"              #13
SIVY_CONFIG_WORD = "configWord"       #14 
SIVY_CHANNEL_BIT_MAP = "channelBitMap"     #15 
SIVY_SAMPLE_BYTES = "sampleBytes"        #16 
                
SIVY_BLOCK_SAMPLES = "blockSamples"     #20 u32 BlockSamples
                
SIVY_SAMPLE_TIME = "sampleTime"        #28 s64 SampleTime;  
SIVY_GPS_TIME = "gpsTime"           #36 s64 GPSTime;        
SIVY_DRIFT = "drift"             #44 s64 Drift;          

SIVY_RSVD0 = "rsvd0"             #48 u32 rsvd0;      

                #/* компас и акселерометр */
SIVY_ACC_Z = "acc_z"             #50 s16 acc_z;  /* 11/50 milli-g, where g is gravitational acceleration*/
SIVY_ACC_X = "acc_x"             #52 s16 acc_x;  /* 12/52 */
SIVY_ACC_Y = "acc_y"              #54 s16 acc_y;  /* 13/54 */
SIVY_COMP_Z = "comp_z"            #56 s16 comp_z; /* 14/56 milli-gauss */
SIVY_COMP_X = "comp_x"           #58 s16 comp_x; /* 15/58 */
SIVY_COMP_Y = "comp_y"           #60 s16 comp_y; /* 16/60 */

                #/* Параметры среды: температура, напряжения и пр */
SIVY_PRESS_REG = "press_reg"        #64 u32 press_reg;  
SIVY_POWER_REG = "power_reg"        #66 u16 power_reg;  
SIVY_TEMPER_REG = "temper_reg"        #68 s16 temper_reg; 
SIVY_HUMIDITY = "humidity"          #70 s16 humidity;   
SIVY_H = "h"                 #72 s16 h;  /* 20/70 height, meter */
SIVY_LAT = "lat"               #76 s32 lat;        /* 21/74 широта (+ N, - S):   55417872  5541.7872N */
SIVY_LON = "lon"             #80 s32 lon;        /* 22/78 долгота(+ E, - W): -37213760  3721.3760W */

OUTPUT_NAMES = [TABLE_SHOT_N, SIVY_SAMPLE_TIME, SIVY_GPS_TIME, SIVY_DRIFT, SIVY_ACC_Z, SIVY_ACC_X, SIVY_ACC_Y, 
                SIVY_COMP_Z, SIVY_COMP_X, SIVY_COMP_Y, SIVY_PRESS_REG, SIVY_POWER_REG, SIVY_TEMPER_REG, 
                SIVY_H, SIVY_LAT, SIVY_LON, TABLE_SAMPLE_RATE, TABLE_SAMPLE_TIME, TABLE_TUNED_DRIFT,
                TABLE_TEMPER_DRIFT, 'starttime', 'tilt', 'azimuth']

# CH_NAME_MAP = {1: 'ch1', 2: 'BHZ', 3: 'BHN', 4: 'BHE', 5: 'ch5', 6: 'ch6'}
# CH_KEYS = sorted(list(CH_NAME_MAP.keys()))
# CH_NAMES = [CH_NAME_MAP[ch_key] for ch_key in CH_KEYS]
# CH_NUM_MAP = {ch_name: ch_key for ch_name, ch_key in zip(CH_NAMES, CH_KEYS)}

DEVS_CHANS_MAP = {'default': {'network': 'CH', 'station': 'LKBD', 'location': 'EN',
                              'chans': ['ch1', 'BHZ', 'BHN', 'BHE', 'ch5', 'ch6']}}


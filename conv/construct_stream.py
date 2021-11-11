'''
Created on 16 дек. 2017 г.

@author: Alexey
'''
import fields_names as fds
from obspy.core.utcdatetime import UTCDateTime
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
from obspy.core.trace import Trace, Stats
from obspy.core import Stream
import logging
import numpy as np

from conv_util import get_ch_num

logger = logging.getLogger(__name__)

def construct_trace(trace_data, user_stats, chan, trace_no=1, rec_no=1):
    stats = dict(channel=str(chan), npts=trace_data.size,
                 sampling_rate=user_stats['sample_rate'], mseed={'dataquality': 'D'})
    dev_dic_entry = fds.DEVS_CHANS_MAP[user_stats['serial']].copy()
    del dev_dic_entry['chans']
    stats.update(dev_dic_entry)
    stats['starttime'] = UTCDateTime(user_stats[fds.SIVY_SAMPLE_TIME] / (10. ** 9))
    stats['user_stats'] = user_stats

    stats['station'] = fds.DEVS_CHANS_MAP[user_stats['serial']]['station']  # user_stats['serial'][-5:]

    theader = SEGYTraceHeader()
    theader.trace_sequence_number_within_line = trace_no 
    theader.trace_sequence_number_within_segy_file = trace_no 
    theader.original_field_record_number = rec_no
    theader.trace_number_within_the_original_field_record = get_ch_num(chan, user_stats['serial'])

    theader.trace_identification_code = 1
    theader.data_use = 1
    theader.datum_elevation_at_receiver_group = user_stats[fds.SIVY_H]
    theader.datum_elevation_at_source = int(user_stats['src_height'])
    theader.scalar_to_be_applied_to_all_coordinates = -10000  #10^4 multiplying factor
    theader.group_coordinate_x = int(user_stats[fds.SIVY_LON] * 10**5)
    theader.group_coordinate_y = int(user_stats[fds.SIVY_LAT] * 10**5)
    theader.coordinate_units = 3    # 1 -- meters, 2 -- arc seconds, 3 -- decimal degrees, 4 -- DMS
    theader.time_basis_code = 2 #utc datetime
    
    serial = user_stats['serial'].strip().encode()[-4:0]
    if len(serial) < 4:
        serial = (4-len(serial)) * b' ' + serial
    theader.shotpoint_number = 0 #np.frombuffer(serial, dtype=np.int32)[0]
    #logger.debug(theader)

    #trace = Trace( data=trace_data, header = stats )
    trace = Trace()
    trace.stats = Stats()
    trace.stats.update(stats)
    trace.stats.segy = {}
    trace.stats.segy.trace_header = theader
    trace.data = trace_data
    #logger.debug('trace start time:' + str(trace.stats.starttime))
    return trace


def construct_stream(chs_data, user_stats, rec_no=1):
    st = Stream()
    chans = user_stats['chans'].split(' ')
    n_of_chans = len(chans)
    for i in range(n_of_chans):
        ch_name = chans[i]
        trace = construct_trace(chs_data[i], user_stats, ch_name, n_of_chans * rec_no + i, rec_no)
        trace.stats.sampling_rate = trace.stats['user_stats'][fds.TABLE_SAMPLE_RATE]
        st += trace
    return st


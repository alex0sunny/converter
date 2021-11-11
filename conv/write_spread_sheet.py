'''
Created on 20 дек. 2017 г.

@author: Alexey
'''
import csv
import logging
import fields_names as fds
from obspy.core.stream import Stream

logger = logging.getLogger(__name__)


def write_spread_sheet(fp, st):
    writer = csv.writer(fp, delimiter=';')
    output_names = fds.OUTPUT_NAMES.copy()
    if not all([angle in st[0].stats['user_stats'].keys() for angle in ['tilt', 'azimuth']]):
        output_names.remove('tilt')
        output_names.remove('azimuth')
    writer.writerow(output_names)
    chan = st[0].stats.channel
    sst = st.select(channel=str(chan))
    for tr in sst:
        user_stats = tr.stats['user_stats']
        #logger.debug('user_stats:' + str(user_stats))
        row = []
        for oname in output_names:
            if oname == fds.TABLE_TUNED_DRIFT and fds.TABLE_TUNED_DRIFT not in \
                    user_stats.keys():
                row.append(0)
                continue
            if oname == fds.TABLE_TEMPER_DRIFT and fds.TABLE_TEMPER_DRIFT not in \
                    user_stats.keys():
                row.append(0)
                continue
            if oname == fds.SIVY_TEMPER_REG:
                row.append(format(user_stats[oname], '4f'))
                continue
            if oname == 'starttime':
                starttime = tr.stats.starttime
                time_str = starttime.datetime.strftime('%Y.%m.%d %H:%M:%S.%f')[:-3]
                row.append(time_str)
                continue
            if oname == fds.TABLE_SAMPLE_RATE:
                row.append(str(int(tr.stats.sampling_rate)))
                continue
            if oname in ['tilt', 'azimuth'] and oname not in user_stats.keys():
                logger.debug('no tilt and azimuth in user_stats')
                row.append('0')
                continue
            row.append(user_stats[oname])
        writer.writerow(row)
    return None


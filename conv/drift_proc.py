'''
Created on 20 дек. 2017 г.

@author: Alexey
'''
import pickle

from obspy import Trace
from obspy.core.stream import Stream, read

from fields_names import SIVY_DRIFT, SIVY_GPS_TIME, SIVY_SAMPLE_TIME, \
    TABLE_TUNED_DRIFT, SIVY_TEMPER_REG, TABLE_TEMPER_DRIFT
# from _collections import OrderedDict
from copy import deepcopy
import logging
from obspy.core.utcdatetime import UTCDateTime
import numpy as np
from conv_util import get_times_map, hFile_gen, get_st_head
import os, gc

logger = logging.getLogger(__name__)


def proc_drift(st_chan, root=None, pv=None, drift_map=None):
    sync_map = {}  # OrderedDict()
    for tr in st_chan:
        user_stats = tr.stats['user_stats']
        gps_time = user_stats[SIVY_GPS_TIME]
        drift = user_stats[SIVY_DRIFT]
        if gps_time == 0 or drift == 0:
            continue
        sync_map[gps_time] = drift
    if not sync_map:
        logger.warning('no sync data')
        return st_chan
    logger.info('sync data is present:' + str(len(sync_map)) + 'pts')
    if drift_map:
        max_time = max(sync_map.keys())
        for gps_time in drift_map:
            if gps_time < max_time + 3600:
                continue
            drift = drift_map[gps_time]
            logger.debug('additional pt will be used:' + str(UTCDateTime(gps_time / (10 ** 9))) +
                         ' drift, s:' + str(drift / (10 ** 9)))
            sync_map[gps_time] = drift
    # logger.debug('sync data:' + str(sync_map))
    for gps_time_ns in sync_map:
        gps_time = UTCDateTime()
        gps_time._ns = gps_time_ns
        # logger.debug('gpsTime=' + str(gps_time) + ' drift=' + str(sync_map[gps_time_ns]))
    # exit(1)
    if len(sync_map) == 1:
        drift = list(sync_map.items())[0][1]
        for tr in st_chan:
            tr.stats['user_stats'][TABLE_TUNED_DRIFT] = drift
            tr.stats['user_stats']['tuned_time'] = tr.stats['user_stats'][SIVY_SAMPLE_TIME] + drift
        return st_chan
    sync_list = []  #list(sync_map.items())
    sync_keys = sorted(list(sync_map.keys()))
    for sync_key in sync_keys:
        sync_list.append([sync_key, sync_map[sync_key]])
    temp_list = deepcopy(sync_list)
    line_list = []

    while len(temp_list) > 1:
        [t1, drift1] = temp_list.pop(0)
        [t2, drift2] = temp_list[0]
        # logger.debug('t1=' + str(UTCDateTime(t1 / (10**9))) + ' t2=' + str(UTCDateTime(t2 / (10 ** 9))) +
        #              '\nd1=' + str(drift1) + ' d2=' + str(drift2))
        k = (drift2 - drift1) / (t2 - t1)
        d0 = int(drift2 - k * t2)
        line_list.append([t1, {'k': k, 'd0': d0}])

    len_secs = st_chan[-1].stats.starttime - st_chan[0].stats.starttime
    for tr in st_chan:
        sample_time = tr.stats['user_stats'][SIVY_SAMPLE_TIME]
        cur_secs = tr.stats.starttime - st_chan[0].stats.starttime
        if root:
            #logger.debug('set progress:' + str(cur_secs * 100 / len_secs))
            pv.set(cur_secs * 100 / len_secs)
            root.update()
        #logger.debug('sample time:' + str(sample_time))
        if len(sync_list) == 1 or sample_time < sync_list[0][0]:
            drift = sync_list[0][1]
        else:
            while len(sync_list) > 1:
                t1 = sync_list[0][0]
                t2 = sync_list[1][0]
                if t1 < sample_time < t2:
                    line_key = sync_list[0][0]
                    line_item = []  # line_list[line_key]
                    #logger.debug('line_key:' + str(line_key))
                    for cur_key, cur_item in line_list:
                        #logger.debug('cur_key:' + str(cur_key))
                        if cur_key == line_key:
                            line_item = cur_item
                            break
                    #logger.debug('line item:' + str(line_item))
                    k = line_item['k']
                    d0 = line_item['d0']
                    drift = int(d0 + k * sample_time)
                    #logger.debug('calc drift,ns:' + str(drift) + \
                    #             ' gps drift:' + str(tr.stats['user_stat'][SIVY_DRIFT]))
                    break
                sync_list.pop(0)
            if len(sync_list) == 1:
                drift = sync_list[0][1]
        tr.stats['user_stats'][TABLE_TUNED_DRIFT] = drift
        # for trace in st_chan:
        #     cur_sample_time = trace.stats['user_stats'][SIVY_SAMPLE_TIME]
        #     if cur_sample_time < sample_time:
        #         continue
        #     if cur_sample_time > sample_time:
        #         break
        #     #trace.stats['user_stats']['tuned_time'] = sample_time + drift
        #     trace.stats['user_stats'][TABLE_TUNED_DRIFT] = drift
    return st_chan


def proc_temper_drift(st_chan):

    time1 = 0
    time2 = 0
    gps_time1 = 0
    for tr in st_chan:
        stats = tr.stats
        gps_time = stats['user_stats'][SIVY_GPS_TIME] / 10**9
        drift = stats['user_stats'][SIVY_DRIFT]
        starttime = tr.stats.starttime
        #logger.debug('gps_time:' + str(UTCDateTime(gps_time).datetime) + ' starttime:' +
        #             str(starttime.datetime))

        # In case of first trace synchronization may be more distant.
        if starttime - 60 < gps_time < starttime or tr == st_chan[0] and gps_time > 0 and \
                starttime - 3000 < gps_time < starttime:
            logger.debug('gps pt:' + str(UTCDateTime(gps_time)))
            if gps_time1:
                logger.debug('prev pt:' + str(UTCDateTime(gps_time1)))
                time2 = tr.stats.endtime
                if not tr.stats.npts:
                    time2 = tr.stats.starttime + 60 - tr.stats.delta
                if time2 > time1 + 300:
                    st_sliced = Stream()
                    for tr in st_chan:
                        if time1 <= tr.stats.starttime < time2:
                            st_sliced += tr
                        if tr.stats.starttime > time2:
                            break
                    if st_sliced[0].stats['user_stats'][SIVY_TEMPER_REG] != \
                            st_sliced[-1].stats['user_stats'][SIVY_TEMPER_REG]:
                        proc_drift_sub(st_sliced)
            gps_time1 = gps_time
            #logger.debug('gps_time1:' + str(UTCDateTime(gps_time1).datetime))
            time1 = starttime
            time2 = 0
    tuned_drift = 0
    temper_drift = 0
    for tr in st_chan:
        stats = tr.stats
        if TABLE_TUNED_DRIFT in stats['user_stats'].keys():
            tuned_drift = stats['user_stats'][TABLE_TUNED_DRIFT]
        else:
            tuned_drift = stats['user_stats'][SIVY_DRIFT]
            stats['user_stats'][TABLE_TUNED_DRIFT] = tuned_drift
        if TABLE_TEMPER_DRIFT in stats['user_stats'].keys():
            temper_drift = stats['user_stats'][TABLE_TEMPER_DRIFT]
        else:
            temper_drift = tuned_drift
            stats['user_stats'][TABLE_TEMPER_DRIFT] = tuned_drift
        #stats['user_stats']['tuned_time'] = stats.starttime + temper_drift / 10**9
        stats['user_stats']['tuned_time'] = stats.starttime + tuned_drift / 10 ** 9
    return st_chan


def proc_drift_sub(st):
    t1 = st[0].stats['user_stats'][SIVY_TEMPER_REG]
    delta = 0
    logger.debug('chunk to process temper drift:' + st.__str__(extended=True))
    for tr, tr_prev in zip(st[1:], st[:-1]):
        t2 = tr.stats['user_stats'][SIVY_TEMPER_REG]
        delta_sub = t2 - t1
        if tr.stats.starttime - tr_prev.stats.starttime > 61:
            # logger.warning('break detected:' + str(tr_prev.stats.starttime) +
            #                '-' + str(tr.stats.starttime))
            delta_sub *= round((tr.stats.starttime - tr_prev.stats.starttime) / 60)
        delta += delta_sub
        # logger.debug('t2:' + str(t2) + ' t1:' + str(t1) + ' cumulative delta t:' + str(delta))
    drift0 = st[0].stats['user_stats'][SIVY_DRIFT]
    delta_d = st[-1].stats['user_stats'][SIVY_DRIFT] - drift0
    #logger.debug('delta_d:' + str(delta_d))
    if delta:
        k = delta_d / delta
        delta = 0
        for tr, tr_prev in zip(st[1:], st[:-1]):
            t2 = tr.stats['user_stats'][SIVY_TEMPER_REG]
            delta_sub = t2 - t1
            if tr.stats.starttime - tr_prev.stats.starttime > 61:
                # logger.warning('break detected:' + str(tr_prev.stats.starttime) +
                #                '\t-\t' + str(tr.stats.starttime))
                delta_sub *= round((tr.stats.starttime - tr_prev.stats.starttime) / 60)
            delta += delta_sub
            tr.stats['user_stats'][TABLE_TEMPER_DRIFT] = drift0 + int(k * delta)
            # logger.debug('starttime:' + str(tr.stats.starttime) + ' drift0:' + str(drift0) + ' k:' +
            #              str(k) + ' tuned drift:' + str(tr.stats['user_stats'][TABLE_TEMPER_DRIFT]))
    return st


def tune_drift(st, root=None, pv=None, drift_map=None):
    logger.info('proc linear drift..')
    proc_drift(st, root, pv, drift_map)
    logger.info('proc temper drift..')
    st = proc_temper_drift(st)
    logger.info('drift processed')
    return st


def tune_starttime(st_head):
    if 'tuned_time' in st_head[0].stats.user_stats.keys():
        logger.debug('initial tuned starttime:' + str(st_head[0].stats.user_stats.tuned_time))
    starttime = st_head[0].stats.starttime
    if 'tuned_time' in st_head[0].stats.user_stats.keys():
        starttime = st_head[0].stats.user_stats.tuned_time
    # starttime = UTCDateTime(starttime.datetime) + 1    # nullify ns
    # starttime.microsecond = 0   # nullify ms
    # logger.debug('corrected starttime:' + str(starttime))
    st_head.sort()
    chan = st_head[0].stats.channel
    nexttime = starttime
    for tr in st_head:
        if 'tuned_time' in st_head[0].stats.user_stats.keys():
            logger.debug('initial tuned starttime:' + str(st_head[0].stats.user_stats.tuned_time))
        starttime = st_head[0].stats.starttime
    return st_head


def interpolate_stream_core(st_chan, starttime=None, endtime=None):
    logger.debug('st_chan before interpolation:' + str(st_chan) + '\nstarttime:' + str(starttime))
    sr_ref = int(st_chan[0].stats.sampling_rate)
    delta = 1 / sr_ref
    npts_ref = sr_ref * 60
    if not starttime:
        starttime = st_chan[0].stats.starttime
    if starttime > st_chan[-1].stats.endtime - 1:
        logger.warning('the starttime ' + str(starttime) +
                       ' is too close to the end of stream:' +
                       str(st_chan[-1].stats.endtime))
        return Stream()
    if len(st_chan) < 3:
        #logger.debug('st_chan:' + str(st_chan))
        st_chan.merge(method=-1, misalignment_threshold=.5)
        #logger.debug('after merging:' + str(st_chan))
        st_out = st_chan.interpolate(sampling_rate=sr_ref, starttime=starttime,
                                     method='lanczos', a=20)
        #logger.debug('st_out:' + str(st_out))
        return st_out
    if st_chan[-1].stats.endtime - starttime < 120:
        st_sliced = st_chan[-3:]
        sr = npts_ref / (st_sliced[-1].stats.starttime -
                         st_sliced[-2].stats.starttime)
        for tr in st_sliced:
            tr.stats.sampling_rate = sr
        st_sliced.interpolate(sampling_rate=sr, method='lanczos', a=20)
        tr = st_sliced[-2].copy()
        tr.data = np.append(tr.data, st_sliced[-1].data)
        st_out = Stream()
        next_time = starttime
        if tr.stats.endtime - next_time > 60:
            logger.debug('next_time:' + str(next_time) + '\nst_sliced:' +
                         str(st_sliced) + '\ntr:' + str(tr))
            st_out += tr.copy().interpolate(starttime=next_time, sampling_rate=sr_ref,
                                            method='lanczos', a=20, npts=npts_ref)
            next_time += 60
        st_out += tr.copy().interpolate(starttime=next_time, sampling_rate=sr_ref,
                                        method='lanczos', a=20)
        return st_out
    logger.debug('st_chan:' + (st_chan[:8] + st_chan[-8:]).__str__(extended=True))
    st_sliced = st_chan.slice(starttime - 61)
    logger.debug('st_sliced:\n' + st_sliced[:12].__str__(extended=True))
    if st_sliced[0].stats.npts < npts_ref:
        st_sliced = st_sliced[1:]
    if st_sliced[-1].stats.npts < npts_ref:
        st_sliced = st_sliced[:-1]
    next_time = starttime
    st_out = Stream()
    logger.debug('st_sliced after cut:\n' + st_sliced[:12].__str__(extended=True))
    for tr, tr_next in zip(st_sliced, st_sliced[1:]):
        while True:
            tr_cur = tr.copy()
            #logger.debug('tr_cur:' + str(tr_cur))
            logger.debug('slow interpolation next time:' + str(next_time))
            if next_time < tr.stats.starttime:
                #logger.debug('go to next time')
                next_time += 60
            if next_time > tr_next.stats.starttime:
                #logger.debug('go to next trace')
                break
            sr = npts_ref / (tr_next.stats.starttime - tr_cur.stats.starttime)
            tr_cur.stats.sampling_rate = sr
            tr_cur.data = np.append(tr_cur.data, tr_next.data)
            tr_cur.data = np.append(tr_cur.data, tr_cur.data[-1])  # crutch
            logger.debug('tr_cur:' + str(tr_cur) + '\nstarttime:' + str(next_time))
            tr_cur.interpolate(sampling_rate=sr_ref, starttime=next_time,
                               npts=npts_ref, method='lanczos', a=20)
            st_out += tr_cur
            #logger.debug('st_out:' + st_out[-8:].__str__(extended=True))
            next_time += 60
        if endtime and next_time > endtime + delta / 10:
            logger.debug('break out of cycle, next_time > endtime - delta / 10,\nnext_time:' + str(next_time) +
                         '\nendtime - delta / 10:' + str(endtime - delta / 10))
            break
    return st_out


def interpolate_stream(st, starttime=None):
    logger.debug('stream chunk to interpolate:' + st.__str__(extended=False))
    starttime_ = st[0].stats.starttime
    logger.debug('starttime_:' + str(starttime_))
    chans = [tr.stats.channel for tr in st.slice(starttime_, starttime_ + 1)]
    chans = st[0].stats['user_stats']['chans'].split(' ')
    #st_out = Stream()
    logger.debug('chans:' + str(chans))
    #logger.debug('st before interpolation\n:' + (st[:8] + st[-8:]).__str__(extended=True))
    st_out = interpolate_stream_core(st.select(channel=chans[0]), starttime)
    logger.debug('after first channel interpolation:\n' + str(st_out) + '\n\n' +
                 (st_out[:8] + st_out[-8:]).__str__(extended=True))
    endtime = st_out[-1].stats.endtime + st_out[-1].stats.delta
    for chan in chans[1:]:
        st_temp = interpolate_stream_core(st.select(channel=chan), starttime,
                                          endtime=endtime)
        logger.debug('after channel interpolation:\n' + str(st_temp) + '\n\n' +
                     (st_temp[:8] + st_temp[-8:]).__str__(extended=True))
        st_out += st_temp
    for tr in st_out:
        tr.data = tr.data.astype('float32')
    return st_out


def gen_stream_chunk(starttime=None, threshold=.01, n_of_trs=120):
    st = Stream()
    for _ in range(n_of_trs):
        if st:
            starttime = st[-1].stats.endtime + st[-1].stats.delta * (1 + threshold)
        elif not starttime:
            starttime = UTCDateTime(2000, 1, 1)
        tr = Trace(np.zeros(60000, dtype='float32'))
        tr.stats.station = 'ND000'
        tr.stats.network = 'RU'
        tr.stats.channel = 'X'
        tr.stats.starttime = starttime
        tr.stats.delta = 1 / 1000
        secs = [(time._get_second() + time._get_microsecond() / 1000000)
                for time in tr.times('UTCDateTime')]
        tr.data = np.array([min(val, 60 - val) for val in secs], dtype='float32')
        st += tr
        tr = tr.copy()
        tr.stats.channel = 'Y'
        tr.data = np.multiply(tr.data, 1000).astype('float32')
        st += tr
    return st


def interpolate_h(hFile, root=None, pv=None):
    st_head = get_st_head(hFile, None, root, pv)
    gen_sts = hFile_gen(hFile, root=root, pv=pv)
    starttime = None
    st_prev = Stream()
    for st in gen_sts:
        logger.debug('current chunk:\n' + (st[:8] + st[-8:]).__str__(extended=True))
        # * Interpolate next chunk until:
        #   ** End of original stream is reached -- get out of 'for' cycle.
        #   ** The gap to the end of chunk is less than an hour -- go to the next chunk
        #      -- next 'for' iteration.
        while True:
            if not starttime:
                starttime = UTCDateTime(st[0].stats.starttime.datetime) + 1    # get out of ns
                starttime.microsecond = 0
            logger.debug('st_prev:' + str(st_prev))
            st_out = interpolate_stream(st_prev + st, starttime=starttime)
            #st_cpy = st_out.copy().sort(keys=['starttime'])
            #logger.debug('add waveforms to ds:\n' + (st_cpy[:8] + st_cpy[-8:]).__str__(extended=True))
            hp = open(os.path.splitext(hFile)[0][:-1] + 'p.pickle', 'ab')
            #gc.collect()
            pickle.dump(st_out, hp, protocol=4)
            hp.close()
            #ds.add_waveforms(st_out, tag='processed')
            logger.debug('st_out:' + str(st_out))
            starttime = st_out[-1].stats.endtime + st_out[-1].stats.delta
            logger.debug('interpolating stream, next time:' + str(starttime.datetime))
            if starttime > st_head[-1].stats.starttime - 2 * st_head[-1].stats.delta and \
                    starttime > st[-1].stats.endtime - 1:
                logger.info('Processing is finished. Next time:' + str(starttime.datetime) +
                            ' orig stream last start time:' + str(st_head[-1].stats.starttime.datetime)
                            + ' orig stream end time:' + str(st[-1].stats.endtime.datetime))
                break
            if starttime + 1000 > st[-1].stats.starttime:
                if abs(st[-1].stats.starttime - st_head[-1].stats.starttime) < 60:
                    logger.debug('The current chunk ' + str(st) + ' is the last. Last starttime:' +
                                 str(st_head[-1].stats.starttime.datetime))
                    continue
                logger.debug('End of trace is approaching. Next time:' + str(starttime.datetime) +
                             '. End of chunk:' + str(st[-1].stats.endtime) + '. So go to the next chunk.')
                break
        logger.debug('last st_prev:\n' + st_prev.__str__(extended=True))
        logger.debug('st before st_prev:\n' + st[-20:].__str__(extended=True))
        slice_time = starttime - 2
        if st:
            slice_time = min(starttime, st[-1].stats.starttime) - 2
        logger.debug('slice_time:' + str(slice_time))
        st_prev = (st_prev + st).slice(starttime=slice_time)
    #ds = None

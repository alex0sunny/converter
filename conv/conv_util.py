'''
Created on 9 янв. 2018 г.

@author: Alexey
'''
import csv
import json
import pickle
from os.path import expanduser

import re
import traceback
from time import sleep

import sys

from read_cf import ControlRow
from obspy import Stream, read, Trace
from copy import deepcopy
import logging
from obspy.core.utcdatetime import UTCDateTime
import gc
import numpy as np
from obspy.io.segy.core import _read_su
import os
from obspy.io.segy.segy import SEGYBinaryFileHeader, SEGYTraceHeader
from _io import BufferedWriter
import fields_names as fds
from obspy.core.util import AttribDict
import gc

logger = logging.getLogger('conv_util')


def merge_stream(st, root=None, pv=None):  # returns a maerged copy of stream, may update progress bar
    st_out = Stream()
    starttime = st[0].stats.starttime
    st_duration = st[-1].stats.endtime - starttime
    for i in range(int(st_duration / 600 + 1)):
        logger.debug('start time:' + str(starttime))
        logger.debug('endtime:' + str(starttime + 600))
        st_out += deepcopy(st.slice(starttime, starttime + 600))
        last_stats = st_out[-1].stats
        starttime = last_stats.endtime + last_stats.delta
        st_out._cleanup(misalignment_threshold=.5)
        if root:
            pv.set(int(i / st_duration * 60000))
            root.update()
    return st_out


def extract_from_stream(st, cf_data, shot_len=10, root=None, pv=None):
    ''' * Buffer chunk of data from the overall stream.
        * If buffer is empty, next chunk.
        * If buffer is after last shot, return.
        * Get the next interval from control file.
        * If interval lies after overall stream, return from function.
        * While chunk is empty or lies before the interval, get next chunk.
        * If interval lies before chunk, get next interval.
        * Get data from chunk and append to result.
        * Next interval.'''
    delta = st[0].stats.delta
    npts = int(shot_len // delta) + 1
    logger.debug('npts for shot:' + str(npts))
    step_chunk = 1000
    controlRow = ControlRow()
    st_out = Stream()
    n_of_rows = len(cf_data)
    starttime = st[0].stats.starttime
    endtime = st[-1].stats.endtime
    st_chunk = st.slice(starttime - 1, starttime + step_chunk)
    for i in range(n_of_rows):
        controlRow = cf_data[i]
        shot_time = controlRow.date_time
        if shot_time > endtime:  # end of overall stream is reached
            logger.warning('no data for shots from ' + str(controlRow.date_time) +
                           ' shot n:' + str(controlRow.shot_number))
            break
        if shot_time < st_chunk[0].stats.starttime:
            continue  # shot is before chunk, next shot
        if shot_time + shot_len > st_chunk[-1].stats.endtime:  # shot is after chunk, next chunk
            logger.debug('next chunk for shot time:' + str(shot_time) + ' endtime:' +
                         str(st_chunk[-1].stats.endtime))
            st_chunk = st.slice(shot_time - 1, shot_time + step_chunk)
            if st_chunk:
                th = st_chunk[-1].stats.segy.trace_header
                logger.debug('ch num:' +
                             str(th.trace_number_within_the_original_field_record))
        starttime = shot_time + step_chunk
        while not st_chunk and starttime < endtime:  # empty chunk, find next chunk
            logger.debug('Empty chunk before time ' + str(starttime) + '. Next chunk.')
            st_chunk = st.slice(starttime - 30, starttime + step_chunk)
            starttime += step_chunk
        if not st_chunk:  # still no appropriate chunk, return from function
            logger.debug('No chunk is found for shottime ' + str(shot_time))
            break
        if st_chunk[0].stats.starttime > shot_time:
            logger.warning('No chunk is found for shot time ' + str(shot_time) + '.')
            continue  # current chunk is after shot

        st_shot = st_chunk.slice(shot_time, shot_time + shot_len + delta)
        st_shot.merge()
        # ._cleanup(misalignment_threshold=.5)
        # st_shot = st.slice(shot_time, shot_time + shot_len)\
        #                    ._cleanup(misalignment_threshold=.5)

        for tr in st_shot:
            logger.debug('npts:' + str(npts))
            tr.data = tr.data[:npts]
            tr.stats[fds.TABLE_SHOT_N] = controlRow.shot_number
            # th = tr.stats.segy.trace_header
            # logger.debug('trace header time:' + str(th.hour_of_day) + ':' + str(th.minute_of_hour) \
            #             + ':' + str(th.second_of_minute))
        logger.debug('record n:' + str(i) + ' shot n:' +
                     str(st_shot[0].stats[fds.TABLE_SHOT_N]) + ' shot time:' +
                     str(shot_time))
        st_out += st_shot
        # st_out._cleanup(misalignment_threshold=.5)
        if root:
            pv.set(int(i / n_of_rows * 100))
            root.update()
    return st_out


'''def set_shot_numbers(st,cf_data,root=None,pv=None):
    controlRow = ControlRow()
    j_max = len(cf_data)
    j = 0
    n_of_trs = len(st)
    tr = st[0]
    controlRow = cf_data[0]
    limit_time = st[-1].stats.endtime
    for i in range(n_of_trs):
        if j==j_max: break
        if root:
            pv.set(int(i/n_of_trs*100))
        controlRow = cf_data[j]
        tr = st[i]
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime
        controltime = controlRow.date_time
        if controltime > limit_time: break
        if controltime>starttime:
            if controltime<endtime:
                tr.stats['user_stat'][TABLE_SHOT_N] = controlRow.shot_number
                logger.debug('set shot number ' + str(controlRow.shot_number) + \
                    ' for trace time ' + str(starttime))
            else: continue
        else:
            while controltime<=starttime:
                j += 1
                if j==j_max: break
                controlRow = cf_data[j]
                controltime = controlRow.date_time    
    return st'''


def prepare_for_segy(st, ofile, dir_out, root=None, pv=None):
    '''* Read frist trace.
       * Divide if it is too long.
       * Write trace as segy file so the binary header will be generated.
       * Read binary header from the just created file.
       * Set number of traces in binary header according to the the divider.
       * Overwrite the trace to the segy file one more time.
       * Read next trace from stream.
       * Divide if it's too long.
       * Write the trace to the tempory su.
       * Write su to the end of the target segy.
       * Repeat reading and writing to the segy file for all remaining traces.
       * Read the resulted su file.
       * '''
    th = st[0].stats.segy.trace_header
    logger.debug('trace header time before writing:' + str(th.hour_of_day) + ':' +
                 str(th.minute_of_hour) + ':' + str(th.second_of_minute))

    logger.debug('prepare segy stream')
    temp_segy = dir_out + '/temp.segy'
    if os.path.exists(ofile):
        os.remove(ofile)
    if os.path.exists(temp_segy):
        os.remove(temp_segy)
    fp = open(ofile, 'ab')
    fp.close()
    st_len = len(st)
    st_temp = Stream()
    max_len = 2 ** 15 - 1
    for i in range(st_len):
        trace = st[i]
        k_trace = 1
        if len(trace) > max_len:
            k_trace = int(1 + len(trace) / max_len)
            st_temp += trace / k_trace
        else:
            st_temp += trace
        logger.debug('last trace stats:' + str(st_temp[-1].stats))
        if i == 0:
            st_temp.write(ofile, format="SEGY", data_encoding=5)
            st_temp.stats.binary_file_header.number_of_data_traces_per_ensemble \
                = k_trace * st_len
            st_temp.sort(keys=['starttime']).write(ofile, format='SEGY',
                                                   data_encoding=5)  # write with completely filled binary header
            np.fromfile(ofile, 'uint8')[:3600].tofile(ofile)  # overwrite file so that it will only contain headers
            fp = open(ofile, 'ab')
        if st_temp[-1].stats.endtime - st_temp[0].stats.starttime >= 3600 \
                or i == st_len - 1:
            st_temp = st_temp.sort(keys=['starttime'])
            st_temp.write(temp_segy, format='SEGY', data_encoding=5)
            th = st_temp[0].stats.segy.trace_header
            st_temp.clear()
            logger.debug('before appending..')
            np.fromfile(temp_segy, 'uint8')[3600:].tofile(fp)
            logger.debug('after appending')
        if root:
            pv.set(int(i / st_len * 100))
            root.update()
    fp.close()
    if os.path.exists(temp_segy):
        os.remove(temp_segy)
    return ofile


def get_times_map(st):  # returns time pts from a merged stream
    channel = st[0].stats.channel
    st_sliced = st.select(channel=channel)
    gaps = st_sliced.get_gaps(min_gap=st[0].stats.delta)
    times_map = {}
    starttime = st_sliced[0].stats.starttime
    for gap in gaps:
        endtime = gap[4]
        times_map[starttime._get_ns()] = endtime
        starttime = gap[5]
    times_map[starttime._get_ns()] = st_sliced[-1].stats.endtime
    max_len = 3600 * 4
    for starttime in sorted(times_map.keys()):
        endtime = times_map[starttime]
        limit_time = UTCDateTime(starttime / 10 ** 9) + max_len
        while endtime > limit_time:
            times_map[starttime] = limit_time
            times_map[limit_time._get_ns()] = endtime
            starttime = limit_time._get_ns()
            limit_time += max_len
    return times_map


def get_ch_num(ch_name, serial):
    return fds.DEVS_CHANS_MAP[serial]['chans'].index(ch_name) + 1


def save_to_txt_core(st, channels, ofile, root=None, pv=None):
    fp = open(ofile, 'a')
    sts = []
    for channel in channels:
        st_channel = st.select(channel=channel).copy().merge()
        sts.append(st_channel)
    data2d = []
    for st_channel in sts:  # iterate over channels
        data = np.transpose([st_channel[0].data])
        for tr in st_channel[1:]:
            data = np.append(data, np.transpose([tr.data]), axis=0)
        if data2d == []:
            data2d = data
        else:
            # logger.debug('data:' + str(data))
            # logger.debug('data2d:' + str(data2d))
            data2d = np.append(data2d, data, axis=1)
    np.savetxt(fp, data2d, delimiter='\t', fmt='%13.10f')
    fp.close()
    # exit()

    # for i in range(len(sts[0])):
    #     trs = []
    #     data = np.array([sts[0][i].data])
    #     for st_channel in sts[1:]:
    #         data = np.append(data, [st_channel[i]], axis=0)
    #     data = np.transpose(data)
    #     np.savetxt(fp, data, delimiter='\t', fmt='%13.10f')
    # fp.close()


def save_to_txt(st, channels, ofile, root=None, pv=None):
    st.sort()
    if os.path.exists(ofile):
        os.remove(ofile)
    fp = open(ofile, 'a')
    for channel in channels[:-1]:
        fp.write(channel + '\t')
    fp.write(channels[-1] + '\n')
    fp.close()
    starttime = st[0].stats.starttime
    step = 1000
    i_max = (st[-1].stats.endtime - starttime) // step
    i_max = max(1, i_max)
    i = 0
    pv.set(0)
    root.update()
    st_temp = Stream()
    for tr in st:
        start_time = tr.stats.starttime
        if start_time > starttime + step:
            # st_temp = st_temp.merge()
            starttime = start_time
            print('st_temp before saving to txt:')
            print(st_temp)
            save_to_txt_core(st_temp, channels, ofile)
            st_temp.clear()
            i += 1
            if root:
                pv.set(int(i / i_max * 100))
                root.update()
        st_temp += tr
        if tr == st[-1]:
            break
    if len(st_temp) > 0:
        save_to_txt_core(st_temp, channels, ofile)

    '''while starttime < st[-1].stats.endtime:
        endtime = starttime + step
        st_sliced = st.slice(starttime=starttime, endtime=endtime).copy().merge()
        #logger.debug('stream merged:' + str(st_sliced))
        save_to_txt_core(st_sliced, channels, ofile)
        starttime = st_sliced[-1].stats.endtime + st_sliced[-1].stats.delta
        i += 1
        if root:
            pv.set(int(i / i_max * 100))
            root.update()'''


def save_to_mseed(st, ofile):
    for tr in st:
        tr.stats.station = tr.stats['serial'][-5:]
    st.write(ofile)
    for tr in st:
        tr.stats.station = tr.stats['serial']


def get_ofile_base0(st, mseed=False):
    filetime = st[0].stats.starttime
    filetime += round(filetime.microsecond / 10 ** 6)
    first_part = st[0].stats['user_stats']['serial']
    if mseed:
        first_part = st[0].stats.station
    return first_part + '_' + filetime.datetime.strftime('%y%m%d_%H%M%S')


def get_ofile_base(st, dir_out, mseed=False):
    return dir_out + '/' + get_ofile_base0(st, mseed)


def get_dir_files(ext, dir_out):
    return [dir_out + '/' + x for x in os.listdir(dir_out) if x.endswith('.' + ext)]


def restore_user_stat(st):
    for tr in st:
        stats = tr.stats
        user_stat = {}
        for key in stats:
            if key in fds:
                user_stat[key] = stats[key]
                stats.pop(key)
        stats.update(user_stat)


def conv_attrib_dict(attrib_dict):
    ret_val = {}
    for k, v in attrib_dict.items():
        if isinstance(v, AttribDict):
            ret_val[k] = conv_attrib_dict(v)
            continue
        if isinstance(v, bytes):
            continue
        # logger.debug('key:' + str(k) + '\nvalue:' + str(v))
        ret_val[k] = v
    return ret_val


def temp_sac_conv(dir_out):
    mseed_files = get_dir_files('mseed', dir_out)
    st_mseed = Stream()
    for mseed_file in mseed_files:
        logger.debug('try to read file:' + mseed_file)
        gc.collect()
        sleep(1)
        st_mseed = read(mseed_file)
        st_mseed.sort()
        starttime = st_mseed[0].stats.starttime
        step = 3600 * 2
        while starttime < st_mseed[-1].stats.endtime:
            st_sliced = st_mseed.slice(starttime, starttime + step).merge()
            st_sliced.write(dir_out + '/' + st_sliced[0].stats.station + '_' +
                            st_sliced[0].stats.starttime.datetime.strftime('%y%m%d_%H%M%S')
                            + '.sac')
            logger.debug('stream written:' + str(st_sliced))
            starttime += step
            logger.debug('next starttime:' + str(starttime))


def create_segy(st0, ofile, number_of_traces):
    tr = st0[0]
    max_len = 2 ** 15 - 1
    k_tr = 1
    if tr.stats.npts > max_len:
        k_tr = tr.stats.npts // max_len + 1
        tr = (tr / k_tr)[0]
    st0 = Stream() + tr
    logger.debug('st0:' + str(st0))
    st0.write('temp.segy', data_encoding=5, format='SEGY')
    st0.stats.binary_file_header.number_of_data_traces_per_ensemble = int(min(number_of_traces, max_len))
    # logger.debug('number of traces:' + str(number_of_traces))
    st0.write('temp.segy', data_encoding=5)
    fp = open('temp.segy', 'rb')
    bytes = fp.read(3600)
    fp.close()
    os.remove('temp.segy')
    fp = open(ofile, 'wb')
    fp.write(bytes)
    fp.close()


def append_to_segy(st, ofile):
    st_temp = Stream()
    max_len = 2 ** 15 - 1
    for tr in st:
        if len(tr) > max_len:
            k_tr = int(1 + len(tr) / max_len)
            st_temp += tr / k_tr
        else:
            st_temp += tr
    st_temp.write('temp.segy', data_encoding=5)
    st_temp = read('temp.segy', unpack_trace_headers=True)
    # logger.debug('st_temp last trace data:' + str(st_temp[-1].data))
    fi = open('temp.segy', 'rb')
    fo = open(ofile, 'ab')
    fi.seek(3600)
    fo.write(fi.read())
    fi.close()
    fo.close()
    # logger.debug('temp.segy size:' + str(os.path.getsize('temp.segy')))
    # logger.debug(ofile + ' size:' + str(os.path.getsize(ofile)))
    # st_temp = read(ofile)
    # logger.debug('data in the last trace of ' + ofile + ':\n' + str(st_temp[-1].data))
    # exit(1)
    # os.remove('temp.segy')


def hFile_gen(hFile, root=None, pv=None):
    logger.debug('first yield')
    hp = open(hFile, 'rb')
    hp.seek(0, os.SEEK_END)
    file_len = hp.tell()
    hp.close()
    hp = open(hFile, 'rb')
    while True:
        st = pickle.load(hp)
        logger.debug('yield st:' + str(st))
        yield st
        if pv and root:
            pv.set(100 * hp.tell() / file_len)
            root.update()
        if hp.tell() == file_len:
            hp.close()
            break


def cfile_gen(hpFile, cfdata, shot_len):

    # Get 2-hour data chunk beginning from the first control row.
    # For each row check if row lies within the chunk.
    # If row within chunk add data to output stream.
    # Update output stream with headers.
    # Update stream with control row data.
    # If row does not fall into chunk, get next chunk.
    # If row is after the end of st_head return from function.
    delta = None
    npts_ref = None
    shot_st_gen = superimposed_st_gen(hpFile)

    try:
        st = next(shot_st_gen)
        for control_row in cfdata:
            if not control_row:
                logger.info('Empty row means next pass. Return empty stream.')
                yield Stream()
                continue
            cf_time = control_row.date_time
            logger.debug('cf_time:' + str(cf_time))
            if not delta:
                delta = st[0].stats.delta
            while st[-1].stats.endtime < cf_time + shot_len:
                st = next(shot_st_gen)
            #logger.debug('st to cut from:' + str(st))
            st_sliced = st.copy()
            npts_offset = round((cf_time - st[0].stats.starttime) / delta)
            #st_sliced = st.slice(cf_time, cf_time + shot_len - delta / 2)
            if not npts_ref:
                npts_ref = round(shot_len * st_sliced[0].stats.sampling_rate)
            for tr in st_sliced:
                tr.data = tr.data[npts_offset: npts_offset + npts_ref]
                tr.stats.starttime += delta * npts_offset
            # if st_sliced[0].stats.npts > npts_ref:
            #     logger.debug('npts ' + str(st_sliced[0].stats.npts) + ' is more than npts_ref ' + str(npts_ref))
            #     for tr in st_sliced:
            #         tr.data = tr.data[:-1]
            for tr in st_sliced:
                k = 10 ** 5
                tr.stats.segy.trace_header.source_coordinate_x = int(control_row.src_lon * k)
                tr.stats.segy.trace_header.source_coordinate_y = int(control_row.src_lat * k)
                tr.stats.segy.trace_header.datum_elevation_at_source = int(control_row.src_height)
                tr.stats.segy.trace_header.shotpoint_number = control_row.shot_number
            st_out = Stream()
            max_len = 2 ** 15 - 1
            for tr in st_sliced:
                npts = tr.stats.npts
                if npts > max_len:
                    k_tr = npts // max_len + 1
                    st_out += tr / k_tr
                else:
                    st_out += tr
            yield st_out
    except StopIteration:
        pass

    shot_st_gen = None


def calc_n_of_tr(st_head, cf_data, shot_len):
    delta = st_head[0].stats.delta
    st_channel = st_head.select(channel=st_head[0].stats.channel)
    npts = shot_len // delta
    k_tr = 1
    max_len = 2 ** 15 - 1
    if npts > max_len:
        k_tr = npts // max_len + 1
    cf_times = [control_row.date_time for control_row in cf_data if control_row]
    n_of_tr = 0
    for cf_time in cf_times:
        if st_channel[0].stats.starttime < cf_time < st_channel[-1].stats.starttime - shot_len:
            n_of_tr += k_tr
            # logger.debug('cf_time ' + str(cf_time) + ' is in range ' + str(st_channel[0].stats.starttime) +
            #              ' - ' + str(st_channel[-1].stats.starttime - shot_len) +
            #              ', n_of_tr:' + str(n_of_tr))
    return n_of_tr


def slice_hfile(hFile, starttime, endtime, foffset=0):
    st_res = Stream()
    hp = open(hFile, 'rb')
    hp.seek(0, os.SEEK_END)
    eof_pos = hp.tell()
    hp.seek(foffset, 0)
    while hp.tell() < eof_pos:
        cur_offset = hp.tell()
        logger.debug('current offset:' + str(cur_offset))
        sth = pickle.load(hp)
        if sth[-1].stats.endtime < starttime:
            foffset = cur_offset
            continue
        if sth[0].stats.starttime > endtime:
            break
        for tr in sth:
            if tr.stats.starttime < starttime < tr.stats.endtime or \
                    tr.stats.starttime < endtime < tr.stats.endtime or \
                    starttime < tr.stats.starttime < endtime or \
                    starttime < tr.stats.endtime < endtime:
                st_res += tr

        # [f, l] = [sth[0].stats.starttime-1 < dt < sth[-1].stats.endtime+1
        #           for dt in [starttime, endtime]]
        # logger.debug('sth:' + str(sth) + '\nstarttime:' + str(starttime) + ' endtime:' + str(endtime) +
        #              ' f, l:' + str(f) + ' ' + str(l))
        # if l:
        #     foffset = cur_offset
        # if f or l or starttime < sth[0].stats.starttime < endtime or \
        #         starttime < sth[-1].stats.starttime < endtime:
        #     for tr in sth:
        #         [ftr, ltr] = [starttime - 1 < dt < endtime + 1
        #                       for dt in [tr.stats.starttime, tr.stats.endtime]]
        #         # logger.debug('ftr:' + str(ftr) + ' ltr:' + str(ltr) + '\ntr:' + str(tr))
        #         if ftr or ltr:
        #             # logger.debug('append tr to st_res:' + str(tr) + '\nftr=' + str(ftr) + ' ltr=' +
        #             #              str(ltr) + ' starttime=' + str(starttime) + ' endtime=' + str(endtime))
        #             st_res += tr
        if sth[0].stats.starttime > endtime + 1:
            break
    hp.close()
    logger.debug('return st_res:' + st_res.__str__(extended=False) + '\nfoffset:' + str(foffset))
    return st_res, foffset


def TrGen(hFile, channel=None):
    hp = open(hFile, 'rb')
    hp.seek(0, os.SEEK_END)
    eof_pos = hp.tell()
    hp.seek(0)
    while hp:
        st = pickle.load(hp)
        if hp.tell() == eof_pos:
            hp.close()
            hp = None
        logger.debug('next st, starttime:' + str(st[0].stats.starttime))
        chans = st[0].stats['user_stats']['chans'].split(' ')
        st_temp = st.select(channel=chans[0])
        for chan in chans[1:]:
            st_temp += st.select(channel=chan)
        st = st_temp.sort(keys=['starttime'])
        for tr in st:
            if channel and channel != tr.stats.channel:
                pass
            yield tr


def get_st_head(hfile, channel=None, root=None, pv=None):
    st_head = Stream()
    hp = open(hfile, 'rb')
    hp.seek(0, os.SEEK_END)
    eof_pos = hp.tell()
    hp.seek(0)
    while hp.tell() < eof_pos:
        st = pickle.load(hp)
        if root and pv:
            pv.set(int(hp.tell() / eof_pos * 100))
            root.update()
        logger.debug('next st, starttime:' + str(st[0].stats.starttime))
        chans = st[0].stats['user_stats']['chans'].split(' ')
        st_temp = st.select(channel=chans[0])
        for chan in chans[1:]:
            st_temp += st.select(channel=chan)
        st = st_temp.sort(keys=['starttime'])
        for tr in st:
            if channel and channel != tr.stats.channel:
                pass
            tr.data = np.zeros(0, 'float32')
        st_head += st
    hp.close()
    return st_head


def fre(sts_group, st_head):
    logger.debug('st_head:' + str(st_head))
    stats1 = sts_group[-1][-1].stats
    ustats1 = stats1['user_stats']
    stats2 = st_head[0].stats
    ustats2 = stats2['user_stats']
    if stats1.station == stats2.station and \
            stats1.sampling_rate == stats2.sampling_rate and ustats1['chans'] == ustats2['chans']:
        # and ustats1[fds.SIVY_GPS_TIME] == ustats2[fds.SIVY_GPS_TIME]:
        sts_group[-1] += st_head
    else:
        sts_group.append(st_head)
        logger.debug('append stream:' + str(st_head))
    return sts_group


def format_exception(e):
    exception_list = traceback.format_stack()
    exception_list = exception_list[:-2]
    exception_list.extend(traceback.format_tb(sys.exc_info()[2]))
    exception_list.extend(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1]))

    exception_str = "Traceback (most recent call last):\n"
    exception_str += "".join(exception_list)
    # Removing the last \n
    exception_str = exception_str[:-1]

    return exception_str


def check_next_st(st_cur, st_next):
    return st_cur and st_next and (st_next[0].stats.starttime > st_cur[0].stats.starttime + 60 * 60 * 2 or
                                   st_next[0].stats.starttime.hour > st_cur[0].stats.starttime.hour + 1)


def superimposed_st_gen_sub(hfile):
    with open(hfile, 'rb') as hp:
        hp.seek(0, os.SEEK_END)
        eof_pos = hp.tell()
        hp.seek(0)
        st_prev = Stream()
        n_of_chs = 0
        chans = None
        while hp.tell() < eof_pos:
            st = pickle.load(hp)
            if not chans:
                chans = st[0].stats['user_stats']['chans'].split(' ')
                n_of_chs = len(chans)
            st_out = Stream()
            for chan in chans:
                st_out += st.select(channel=chan)
            st_out.sort(keys=['starttime'])
            yield st_prev + st_out
            st_prev = st_out[-n_of_chs:]


def superimposed_st_gen(hfile):
    st_gen = superimposed_st_gen_sub(hfile)
    n_of_chs = 0
    for st in st_gen:
        if not n_of_chs:
            n_of_chs = st[0].stats.user_stats.n_of_chs
        st_out = Stream()
        for tr, tr_next in zip(st, st[n_of_chs:]):
            st_out += tr
            logger.debug('adding tr:' + str(tr))
            st_out += tr_next
            logger.debug('adding next tr:' + str(tr_next) +
                         '\ncurrent st_out:' + st_out.__str__(extended=True))
            if not len(st_out) % (2 * n_of_chs):
                st_out.merge()
                logger.debug('yield st:' + st_out.__str__(extended=True))
                yield st_out
                st_out = Stream()
    if st_out:
        yield st_out.merge()


def get_log_map(log_file):
    #log_file = 'D:/converter_data/LOG.TXT'
    drift_map = {}
    with open(log_file) as fp:
        reader = csv.reader(fp, delimiter=' ')
        for row in reader:
            if len(row) < 4:
                continue
            if row[-1] == 's' and row[-3] == 'drift:':
                drift = int(row[-2])
            drift = row[-2]
            drift_str = row[-1]
            pattern = re.compile(r'^.*\ds$')
            if drift != 'drift' or not pattern.match(drift_str):
                continue
            date_str = row[0] + ' ' + row[1]
            gps_time = UTCDateTime.strptime(date_str, '%d-%m-%y %H:%M:%S')
            drift = int(drift_str[:-1])
            drift_map[gps_time._ns] = drift * 10 ** 9
    return drift_map


def get_new_map(log_file):
    drift_map = {}
    with open(log_file) as fp:
        reader = csv.reader(fp, delimiter=' ')
        for row in reader:
            if len(row) < 4 or row[-1] != 's' or row[-3] != 'drift:':
                continue
            drift = int(row[-2])
            date_str = row[0] + ' ' + row[1]
            gps_time = UTCDateTime.strptime(date_str, '%y-%m-%d %H:%M:%S')
            drift_map[gps_time._ns] = drift * 10 ** 9
    return drift_map


def load_chans_csv(file_path):
    rows = []
    if os.path.exists(file_path):
        with open(file_path, newline='') as csvfile:
            reader = csv.reader(csvfile, csv.excel, delimiter=';')
            for row in reader:
                rows.append(row)
    return rows


def save_chans_csv(rows, file_path):
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect=csv.excel, delimiter=';')
        writer.writerows(rows)


def rows_to_devmap(rows):
    devmap = {}
    for row in rows:
        dev = row[0]
        devmap[dev] = {}
        if not dev:
            break
        devmap[dev]['network'] = row[1]
        devmap[dev]['station'] = row[2]
        devmap[dev]['location'] = row[3]
        devmap[dev]['chans'] = row[4:]
    return devmap


def devmap_to_rows(rowmap):
    row = ['default']
    default_entry = rowmap['default']
    logger.debug('default entry:' + str(default_entry))
    row.append(default_entry['network'])
    #row.append(default_entry['station'])
    row.append('')
    row.append(default_entry['location'])
    for cell_val in default_entry['chans']:
        row.append(cell_val)
    rows = [row]
    devs = sorted(list(rowmap.keys()))
    devs.remove('default')
    for dev in devs:
        row = [dev]
        row_entry = rowmap[dev]
        row.append(row_entry['network'])
        row.append(row_entry['station'])
        row.append(row_entry['location'])
        for cell_val in row_entry['chans']:
            row.append(cell_val)
        rows.append(row)
    return rows


def get_app_dir():
    home = expanduser("~")
    app_dir = os.path.join(home, '.conv').replace('\\', '/')
    return app_dir


chans_csv_file = get_app_dir() + '/chans.csv'


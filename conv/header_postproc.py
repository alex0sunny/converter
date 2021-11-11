'''
Created on 15 дек. 2017 г.

@author: Alexey
'''
import BitVector as bv
import fields_names as fds
# from _collections import OrderedDict
from obspy.core.utcdatetime import UTCDateTime
import logging
from conv_util import save_chans_csv, devmap_to_rows, chans_csv_file

logger = logging.getLogger(__name__)


def int_to_codes(int_val):
    bit_data = bv.BitVector(intVal=int_val, size=32)[1:31].reverse()
    chs_bit_data = []
    for i in range(0, 30, 5):
        chs_bit_data.append(bit_data[i:i + 5])

    chs_codes_dics = []
    # gains = [1, 2, 4, 8, 12]
    gains = {0: 1, 1: 1, 2: 2, 4: 4, 5: 8, 6: 12}
    for ch_bit_data in chs_bit_data:
        ch_dic = {}
        ch_dic['enable'] = ch_bit_data[0:1].int_val()
        ch_dic['range'] = ch_bit_data[1:2].int_val()
        gain_bit_data = bv.BitVector(size=4)
        gain_bit_data[1:4] = ch_bit_data[2:5].reverse()
        gain_code = gain_bit_data.int_val()
        # logger.debug('gain code:' + str(gain_code))
        # logger.debug('gain bit data:' + str(gain_bit_data) + ' gain:' + str(gains[gain_code]))
        ch_dic['gain'] = gains[gain_code]
        chs_codes_dics.append(ch_dic)

    chs_ks = []  # these ks are still not used currently
    for ch_dic in chs_codes_dics:
        k = 0.
        if ch_dic['enable'] == 1:
            k = 1.024 / ch_dic['gain']
            if ch_dic['range'] == 1:
                k *= 12.
            else:
                k *= 2.
            # logger.debug('k=' + str(k) + ' k to append:' + str(k*2))
        chs_ks.append(k * 2.)  # last multiplication is somewhat mysterious
    return chs_codes_dics


def byte_to_bits(byte_val):
    bit_data = bv.BitVector(intVal=byte_val, size=8).reverse()
    chs_en = []
    for i in range(8):
        ch_en = bit_data[i:i + 1].int_val()
        chs_en.append(ch_en)
    return chs_en


def get_chs_inds(bit_list):
    inds = []
    for i, x in enumerate(bit_list):
        if (1 == x):
            inds.append(i)
    return inds


def get_chs_dic(checked_mask, ranges_mask, ch_codes_dics, serial, status_var=None):
    inds = get_chs_inds(byte_to_bits(checked_mask))
    # logger.debug('inds:' + str(inds))
    ranges = byte_to_bits(ranges_mask)
    chs_ks = []
    for i in inds:
        ch_k = 2.0
        if 1 == ranges[i]: ch_k = 12.0
        ch_k /= ch_codes_dics[i]['gain']
        # logger.debug('gain:' + str(ch_codes_dics[i]['gain']))
        ch_k *= 1.024 / (2 ** 31 - 1)
        chs_ks.append(ch_k * 2.)  # unclear multiplication
    chs_dic = {}
    chs_dic['n_of_chs'] = len(inds)
    if serial not in fds.DEVS_CHANS_MAP:
        if status_var:
            status_var.set(serial + ': NO DEVICE DATA!')
        fds.DEVS_CHANS_MAP[serial] = fds.DEVS_CHANS_MAP['default'].copy()
        fds.DEVS_CHANS_MAP[serial]['station'] = serial[-5:]
        save_chans_csv(devmap_to_rows(fds.DEVS_CHANS_MAP), chans_csv_file)
    else:
        if status_var.get()[:len(serial)] != serial:
            dev_entry = fds.DEVS_CHANS_MAP[serial]
            status_var.set(serial + ': ' + dev_entry['network'] + '/' + dev_entry['station'] + '/' +
                           dev_entry['location'])
    ch_ref_names = fds.DEVS_CHANS_MAP[serial]['chans']
    chans = [ch_ref_names[i] for i in inds]
    chs_dic['chans'] = ' '.join(str(chan) for chan in chans)
    # logger.debug('chans:' + chs_dic['chans'])
    # exit(1)
    chs_dic['ks'] = ' '.join(str(ch_k) for ch_k in chs_ks)
    # logger.debug('ks:' + str(chs_ks))
    # exit()
    return chs_dic


def header_to_user_dic(binHeader):
    udic = {}
    '''for fname, ftype in binHeader._fields_:
        if fname in fds.OUTPUT_NAMES:
            udic[fname] = getattr(binHeader,fname)'''
    for fname, _ in binHeader._fields_:
        udic[fname] = getattr(binHeader, fname)
        #logger.debug(f'field name: {fname} value:{udic[fname]}')

    udic[fds.SIVY_LAT] /= 10. ** 7
    udic[fds.SIVY_LON] /= 10. ** 7
    udic[fds.SIVY_TEMPER_REG] /= 10.
    udic[fds.TABLE_SAMPLE_RATE] = round(getattr(binHeader, fds.SIVY_BLOCK_SAMPLES) / 60)
    if udic[fds.SIVY_SAMPLE_TIME] < 0:
        udic[fds.SIVY_SAMPLE_TIME] = 0
    sivy_sample_time = UTCDateTime(getattr(binHeader, fds.SIVY_SAMPLE_TIME) / (10. ** 9))
    udic[fds.TABLE_SAMPLE_TIME] = \
        sivy_sample_time.datetime.strftime('%Y.%m.%d %H:%M:%S.%f')[:-3]
    return udic


def header_to_chs_stat(binHeader, serial, status_var=None):
    checked_mask = getattr(binHeader, fds.SIVY_CHANNEL_BIT_MAP)
    ranges_mask = getattr(binHeader, fds.SIVY_CONFIG_WORD)
    gains_mask = getattr(binHeader, fds.SIVY_RSVD0)
    ch_codes_dics = int_to_codes(int(gains_mask))
    chs_stat = get_chs_dic(checked_mask, ranges_mask, ch_codes_dics, serial, status_var)
    chs_stat['n_of_counts'] = int(getattr(binHeader, fds.SIVY_BLOCK_SAMPLES))
    chs_stat['capacity'] = int(getattr(binHeader, fds.SIVY_SAMPLE_BYTES) / chs_stat['n_of_chs'])
    chs_stat['n_of_samples'] = int(chs_stat['n_of_counts'] * chs_stat['n_of_chs'])
    return chs_stat


def header_to_user_stat(binHeader, serial, status_var=None):
    user_stats = header_to_user_dic(binHeader)
    chs_stats = header_to_chs_stat(binHeader, serial, status_var)
    user_stats.update(chs_stats)
    return user_stats

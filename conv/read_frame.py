'''
Created on 15 дек. 2017 г.

@author: Alexey
'''
from time import sleep
from tkinter import messagebox

from conv_util import format_exception
from header_structs import *
from header_postproc import header_to_user_dic, header_to_chs_stat, header_to_user_stat
import logging
from read_frame_data import read_frame_data
import numpy as np
from construct_stream import construct_stream
import os
from rotation import rotate

logger = logging.getLogger(__name__)


def read_frame(fp, rec_no, rotation, status_var=None):
    txtHeader = TxtHeader()
    if fp.readinto(txtHeader) < 2:
        return None
    binHeader = BinHeader()
    if fp.readinto(binHeader) < 10:
        return None

    try:
        nd = txtHeader.txt_header.decode('utf8')
    except UnicodeDecodeError as e:
        logger.error(format_exception(e))
        return None
    if nd != 'NS' and nd != 'ND':
        error_mes = 'Incorrect header bytes:' + nd + '. Expected:NS or ND.\nFile:' + fp.name + \
                    '\nsize:' + str(os.stat(fp.name).st_size) + ' position:' + str(fp.tell()) + '.'
        #messagebox.showerror('Error', error_mes)
        logger.error(error_mes)
        return None
        exit(1)
    serial = txtHeader.serial.decode()

    user_stats = header_to_user_stat(binHeader, serial, status_var)
    user_stats['serial'] = serial
    user_stats['src_lon'] = 0.
    user_stats['src_lat'] = 0.
    user_stats['src_height'] = 0.
    user_stats[fds.TABLE_SHOT_N] = 0
    #logger.debug("user stat:" + str(user_stat))
    #exit()
    n_of_chs = user_stats['n_of_chs']
    data_size = user_stats['n_of_samples'] * user_stats['capacity']
    chs_data_int = []
    if user_stats['capacity'] == 3:
        chs_data_int = read_frame_data(fp, data_size, n_of_chs)
    else:
        if user_stats['capacity'] == 4:
            chs_data_int = read_frame_data(fp, data_size, n_of_chs, i24=False)
        else:
            logger.error('unexpected capacity:' + str(user_stats['capacity']))
            return None
    if chs_data_int.size == 0:
        return None
    chs_data = []
    ks = [float(k) for k in user_stats['ks'].split(' ')]
    #logger.debug('ks:' + str(ks))
    #exit()
    for i in range(n_of_chs):
        ch_data = np.multiply(chs_data_int[i], ks[i]).astype(np.float32)
        chs_data.append(ch_data) 
    #logger.debug("channels data:" + str(chs_data))
    
    stream = construct_stream(chs_data, user_stats, rec_no)
    if rotation:
        stream = rotate(stream)
    return stream


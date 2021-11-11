'''
Created on 19 дек. 2017 г.

@author: Alexey
'''

import csv
from tkinter import messagebox

from obspy import UTCDateTime
import logging
import os

logger = logging.getLogger(__name__)


class ControlRow:
    shot_number = 0
    date_time = UTCDateTime(0)
    src_lon = 0.
    src_lat = 0.
    src_height = 0.


# 'D:/temp/sivy/in/Giorgos/shotforconv_new.txt'
def read_cf(file_path):
    while True:
        # remove trailing newlines from file
        fp = open(file_path, 'r')
        fp.seek(0, os.SEEK_END)
        eof_pos = fp.tell()
        fp.seek(eof_pos-1)
        if fp.read() == '\n':
            fp.seek(0)
            content = fp.read()[:-1]
            fp.close()
            fp = open(file_path, 'w')
            fp.write(content)
            fp.close()
        else:
            fp.close()
            break
    with open(file_path) as fp:
        reader = csv.reader(fp, delimiter='\t')
        control_rows = []
        for row in reader:
            if not row:
                control_rows.append([])
                continue
            date_str = row[2]
            date_list = date_str.split('/')
            date_list.reverse()
            date_str = '-'.join(date_list)
            dt = UTCDateTime(date_str + 'T' + row[1])
            control_row = ControlRow()
            control_row.shot_number = int(row[0])
            control_row.date_time = dt
            if len(row) > 3:
                control_row.src_lon = float(row[3])
            if len(row) > 4:
                control_row.src_lat = float(row[4])
            if len(row) > 5:
                control_row.src_height = float(row[5])
            # logger.debug('control row:' + str(control_row.shot_number) + ' ' + str(control_row.date_time))
            error_mes = ''
            if not -180 < control_row.src_lat < 180:
                error_mes = 'src_lat should be between in degrees, between -180 and 180.\n' \
                      'The lattitude value in control file is ' + str(control_row.src_lat)
            if not -180 < control_row.src_lon < 180:
                error_mes = 'src_lon should be between in degrees, between -180 and 180.\n' \
                      'The longitude value in control file is ' + str(control_row.src_lon)
            if not -2**31 < control_row.src_height < 2**31:
                error_mes = 'Height value in control file ' + str(control_row.src_height) + \
                            ' is out of range.'
            if error_mes:
                logger.error(error_mes)
                messagebox.showerror('Error', error_mes)
                exit(1)
            control_rows.append(control_row)
            #logger.debug('n of rows:' + str(len(control_rows)))
    return control_rows

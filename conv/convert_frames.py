'''
Created on 18 дек. 2017 г.

@author: Alexey
'''
import json
import pickle
import logging
import tkinter
from tkinter import filedialog

#import tk
from functools import reduce
from operator import concat
from trace import Trace

from read_frame import read_frame
from obspy.core.stream import Stream, read

from write_spread_sheet import write_spread_sheet
from copy import deepcopy
import numpy as np
import os
from tkinter import *
from tkinter.constants import LEFT, END
import configparser
import re
from tkinter.ttk import Progressbar
from conv_util import extract_from_stream, prepare_for_segy, get_times_map, save_to_txt, \
    save_to_mseed, get_ofile_base, get_dir_files, temp_sac_conv, get_ofile_base0, hFile_gen, create_segy, \
    append_to_segy, save_to_txt_core, cfile_gen, calc_n_of_tr, slice_hfile, TrGen, get_st_head, fre, get_log_map, \
    get_new_map, get_app_dir, load_chans_csv, rows_to_devmap, get_ch_num, devmap_to_rows, chans_csv_file, save_chans_csv
from read_cf import read_cf
from drift_proc import tune_drift, interpolate_h, tune_starttime
from obspy.core.utcdatetime import UTCDateTime
import fields_names as fds
from fields_names import DEVS_CHANS_MAP
from time import sleep

import gc
import pickle

logging.basicConfig(format='%(levelname)s %(asctime)s %(funcName)s %(filename)s:%(lineno)d '
                           '%(message)s',
                    level=logging.DEBUG)
# filename=os.getcwd() + '/log.txt',
# filemode='w')
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logger = logging.getLogger(__name__)

app_dir = get_app_dir()
if not os.path.isdir(app_dir):
    os.mkdir(app_dir)
config_file = os.path.join(app_dir, 'config.ini').replace('\\', '/')

# chans_csv_file = get_app_dir() + '/chans.csv'

fds.DEVS_CHANS_MAP = rows_to_devmap(load_chans_csv(chans_csv_file))
if 'default' not in fds.DEVS_CHANS_MAP:
    fds.DEVS_CHANS_MAP['default'] = {}
    fds.DEVS_CHANS_MAP['default']['chans'] = ['ch1', 'BHZ', 'BHN', 'BHE', 'ch5', 'ch6']
    fds.DEVS_CHANS_MAP['default'].update({'network': 'CH', 'station': 'LKBD', 'location': 'EN'})

ok_pressed = -1
gps_pt = {}

rotation = True


class DialogData:
    dir_in = ''
    use_file = False
    dir_out = ''
    cb = False
    cfile = ''
    corr = False
    shot_len = 10
    cb_txt = False
    cb_bin = True
    cb_segy = False
    cb_mseed = False
    cb_sac = False


dd = DialogData()


class MyDialog:
    use_file = False

    def __init__(self, parent):

        top = self.top = Toplevel(parent)
        top.resizable(width=False, height=False)

        Label(top, text="input:").grid(row=0, columnspan=2)
        self.e = Entry(top, width=50)
        self.e.grid(row=0, column=2, columnspan=4, sticky=tkinter.W)
        sel_button = Button(top, text='file', command=self.select_in_file)
        sel_button.grid(row=0, column=6)
        sel_button = Button(top, text='dir_in', command=self.select_input)
        sel_button.grid(row=0, column=7)

        Label(top, text="output:").grid(row=1, columnspan=2)
        self.check_drift_var = IntVar()
        self.cb_corr = Checkbutton(top, text='drift corr', variable=self.check_drift_var)
        self.cb_corr.grid(row=1, column=2)
        self.eo = Entry(top, width=40)
        self.eo.grid(row=1, column=3, columnspan=4, sticky=tkinter.W)
        sel_button2 = Button(top, text='dir_out', command=self.select_output)
        sel_button2.grid(row=1, column=7)

        self.check_var = IntVar()
        self.check_button = Checkbutton(top, text='control file', variable=self.check_var)
        self.check_button.grid(row=2, columnspan=2)
        self.ec = Entry(top, width=50)
        self.ec.grid(row=2, column=2, columnspan=4)
        sel_cf = Button(top, text='select', command=self.select_file)
        sel_cf.grid(row=2, column=6)
        chans_button = Button(top, text='chans', command=self.set_chans)
        chans_button.grid(row=2, column=7)


        # chans_cf = Button(top, text='chans')
        # chans_cf.grid(row=2, column=7)


        self.check_txt_var = IntVar()
        self.cb_txt = Checkbutton(top, text='txt', variable=self.check_txt_var)
        self.cb_txt.grid(row=3, column=2)

        self.bcheck_var = IntVar()
        self.cb_bin = Checkbutton(top, text='bin', variable=self.bcheck_var)
        self.cb_bin.grid(row=3, column=3)

        self.check_segy_var = IntVar()
        self.cb_segy = Checkbutton(top, text='segy', variable=self.check_segy_var)
        self.cb_segy.grid(row=3, column=4)

        self.check_mseed_var = IntVar()
        self.cb_mseed = Checkbutton(top, text='mseed', variable=self.check_mseed_var)
        self.cb_mseed.grid(row=3, column=5)

        self.check_sac_var = IntVar()
        self.cb_sac = Checkbutton(top, text='sac', variable=self.check_sac_var)
        self.cb_sac.grid(row=3, column=6)

        b = Button(top, text="OK", command=self.ok)
        b.grid(row=3, column=7)

        Label(top, text='shot,s:').grid(row=3)
        self.shot_len = Entry(top, width=3)
        self.shot_len.grid(row=3, column=1, sticky=tkinter.W)

        config_parser = configparser.ConfigParser()
        if os.path.exists(config_file):
            config_parser.read(config_file)
        config = config_parser['DEFAULT']
        global ok_pressed
        if ok_pressed == -1:
            if 'gps_time' in config and 'drift' in config:
                gps_pt[UTCDateTime.strptime(config['gps_time'], '%d.%m.%Y %H:%M:%S')._ns] = \
                    int(float(config['drift']) * 10 ** 9)
            else:
                config['gps_time'] = UTCDateTime(0).strftime('%d.%m.%Y %H:%M:%S')
                config['drift'] = str(float(0))
                with open(config_file, 'w') as cfile:
                    config_parser.write(cfile)
            if 'input' in config:
                input = config['input']
                self.e.delete(0, END)
                self.e.insert(0, input)
            if 'output' in config:
                output = config['output']
                self.eo.delete(0, END)
                self.eo.insert(0, output)
            if 'control_file' in config:
                control_file = config['control_file']
                self.ec.delete(0, END)
                self.ec.insert(0, control_file)
            if 'cb' in config and config['cb'] == 'true':
                self.check_button.select()
            if 'cb_corr' in config and config['cb_corr'] == 'true':
                self.cb_corr.select()
            else:
                self.cb_corr.deselect()
            if 'shot_len' in config:
                shot = config['shot_len']
                self.shot_len.delete(0, END)
                self.shot_len.insert(0, shot)
            if 'cb_txt' in config and config['cb_txt'] == 'true':
                logger.debug('initially select txt checkbox')
                self.cb_txt.select()
            else:
                logger.debug('initially deselect txt checkbox')
                self.cb_txt.deselect()
            if 'cb_bin' in config and config['cb_bin'] == 'true':
                logger.debug('initially select bin checkbox')
                self.cb_bin.select()
            else:
                logger.debug('initially deselect bin checkbox')
                self.cb_bin.deselect()
            if 'cb_mseed' in config and config['cb_mseed'] == 'true':
                self.cb_mseed.select()
            else:
                self.cb_mseed.deselect()
            if 'cb_sac' in config and config['cb_sac'] == 'true':
                self.cb_sac.select()
            else:
                self.cb_sac.deselect()
            if 'cb_segy' in config and config['cb_segy'] == 'true':
                self.cb_segy.select()
            else:
                self.cb_segy.deselect()
            '''logger.debug('bcheck_var:' + str(self.bcheck_var.get()) + ' check_segy_var:' + self.check_segy_var.get() \ 
                         + ' check_mseed_var:' + str(self.check_mseed_var.get()) + ' check_txt_var:' \
                         + self.check_txt_var)'''
            if self.bcheck_var.get() == 0 and self.check_segy_var.get() == 0 and self.check_mseed_var.get() == 0 \
                    and self.check_txt_var.get() == 0 and self.check_sac_var.get() == 0:
                logger.debug('select bin checkbox since no other boxes are checked')
                self.cb_bin.select()
            if self.shot_len.get() == '': self.shot_len.insert(0, '10')
        else:
            # global ok_pressed
            self.from_dd(dd)
            ok_pressed = -1


    def to_dd(self):
        dd = DialogData()
        dd.dir_in = self.e.get()
        dd.dir_out = self.eo.get()
        dd.cb = self.check_var.get()
        dd.shot_len = int(self.shot_len.get())
        dd.cfile = self.ec.get()
        dd.corr = self.check_drift_var.get()
        dd.cb_txt = self.check_txt_var.get()
        dd.cb_bin = self.bcheck_var.get()
        dd.cb_segy = self.check_segy_var.get()
        dd.cb_mseed = self.check_mseed_var.get()
        dd.cb_sac = self.check_sac_var.get()
        if dd.cb:
            dd.shot = int(self.shot_len.get())
        dd.use_file = os.path.isfile(dd.dir_in)
        return dd

    def from_dd(self, dd):
        self.e.delete(0, END)
        self.e.insert(0, dd.dir_in)
        self.eo.delete(0, END)
        self.eo.insert(0, dd.dir_out)
        self.ec.delete(0, END)
        self.ec.insert(0, dd.cfile)
        if dd.cb:
            self.check_button.select()
        else:
            self.check_button.deselect()
        if dd.corr:
            self.cb_corr.select()
        else:
            self.cb_corr.deselect()
        self.shot_len.delete(0, END)
        self.shot_len.insert(0, str(dd.shot_len))
        if dd.cb_txt:
            logger.debug('initially select txt checkbox')
            self.cb_txt.select()
        else:
            logger.debug('initially deselect txt checkbox')
            self.cb_txt.deselect()
        if dd.cb_bin:
            logger.debug('initially select bin checkbox')
            self.cb_bin.select()
        else:
            logger.debug('initially deselect bin checkbox')
            self.cb_bin.deselect()
        if dd.cb_mseed:
            self.cb_mseed.select()
        else:
            self.cb_mseed.deselect()
        if dd.cb_sac:
            self.cb_sac.select()
        else:
            self.cb_sac.deselect()
        if dd.cb_segy:
            self.cb_segy.select()
        else:
            self.cb_segy.deselect()

    def ok(self):
        status_var.set('')
        logger.debug('bin check var:' + str(self.bcheck_var.get()))
        config_parser = configparser.ConfigParser()
        config = config_parser['DEFAULT']
        config['input'] = self.e.get()
        config['output'] = self.eo.get()
        config['control_file'] = self.ec.get()
        config['shot_len'] = self.shot_len.get()
        logger.debug('check var:' + str(self.check_var.get()))
        if self.check_var.get():
            config['cb'] = 'true'
            config['shot_len'] = self.shot_len.get()
        else:
            config['cb'] = 'false'
        logger.debug("value is" + str(self.e.get()))

        if self.check_drift_var.get():
            config['cb_corr'] = 'true'
        else:
            config['cb_corr'] = 'false'

        if self.bcheck_var.get() == 0 and self.check_segy_var.get() == 0 and \
                self.check_mseed_var.get() == 0 and self.check_sac_var.get() == 0 \
                and self.check_txt_var == 0:
            logger.debug('select bin checkbox programmatically')
            self.cb_bin.select()

        if self.bcheck_var.get():
            logger.debug('bin checkbox is selected. cb_bin:' + str(self.bcheck_var.get()))
            config['cb_bin'] = 'true'
        else:
            logger.debug('bin check bos is deselected:' + str(self.bcheck_var.get()))
            config['cb_bin'] = 'false'

        if self.check_segy_var.get():
            config['cb_segy'] = 'true'
        else:
            config['cb_segy'] = 'false'

        if self.check_mseed_var.get():
            config['cb_mseed'] = 'true'
        else:
            config['cb_mseed'] = 'false'

        if self.check_sac_var.get():
            config['cb_sac'] = 'true'
        else:
            config['cb_sac'] = 'false'

        if self.check_txt_var.get():
            config['cb_txt'] = 'true'
        else:
            config['cb_txt'] = 'false'
        if gps_pt:
            config['gps_time'] = UTCDateTime(list(gps_pt.keys())[0] / (10**9)).strftime('%d.%m.%Y %H:%M:%S')
            config['drift'] = str(list(gps_pt.values())[0] / (10**9))
        with open(config_file, 'w') as cfile:
            config_parser.write(cfile)

        global dd
        dd = self.to_dd()

        global ok_pressed
        ok_pressed = 1

        self.top.destroy()

    def select_input(self):
        self.use_file = False
        dir_opt = {}
        initialdir = self.e.get()
        if os.path.isfile(initialdir):
            initialdir = os.path.dirname(initialdir)
        dir_opt['initialdir'] = initialdir
        dir_in = filedialog.askdirectory(**dir_opt)
        if dir_in == '':
            dir_in = dir_opt['initialdir']
        self.e.delete(0, END)
        self.e.insert(0, dir_in)
        return dir_in

    def select_in_file(self):
        self.use_file = True
        # dir_opt = {}
        initialdir = self.e.get()
        if os.path.isfile(initialdir):
            initialdir = os.path.dirname(initialdir)
        # dir_opt['initialdir'] = initialdir
        file_in = filedialog.askopenfilename(initialdir=initialdir)
        # files = filedialog.askopenfilenames(initialdir=initialdir)
        # files = sorted(list(files))
        # if not files:
        if file_in == '':
            logger.error('no input file is selected')
            exit(1)
        self.e.delete(0, END)
        self.e.insert(0, file_in)
        return file_in

    def select_output(self):
        dir_opt = {}
        dir_opt['initialdir'] = self.eo.get()
        dir_out = filedialog.askdirectory(**dir_opt)
        if dir_out == '':
            dir_out = self.eo.get()
        self.eo.delete(0, END)
        self.eo.insert(0, dir_out)
        return dir_out

    def select_file(self):
        FILEOPENOPTIONS = dict(defaultextension='.txt', \
                               filetypes=[('control file', '*.txt'), ('All files', '*.*')])
        FILEOPENOPTIONS['initialdir'] = os.path.dirname(self.ec.get())
        file_name = filedialog.askopenfilename(**FILEOPENOPTIONS)
        if file_name == '':
            file_name = self.ec.get()
        self.ec.delete(0, END)
        self.ec.insert(0, file_name)
        return file_name

    def set_chans(self):
        global ok_pressed
        ok_pressed = 0
        global dd
        dd = self.to_dd()
        logger.debug('dd saved to struct:' + str(dd))
        self.top.destroy()
        # top = Toplevel()
        # top.title("About this application...")
        # msg = tk.Message(top, text='about_message')
        # msg.pack()
        #
        # button = Button(top, text="Dismiss", command=top.destroy)
        # button.pack()


class ChannelsDialog:

    def __init__(self, parent):

        top = self.top = Toplevel(parent)
        self.n_of_chans = 6
        self.n_of_rows = 20
        self.labels = []
        self.labels.append(Label(top, text='devices'))
        self.labels[-1].grid(column=0, row=0)
        self.labels.append(Label(top, text='network'))
        self.labels[-1].grid(column=1, row=0)
        self.labels.append(Label(top, text='station'))
        self.labels[-1].grid(column=2, row=0)
        self.labels.append(Label(top, text='location'))
        self.labels[-1].grid(column=3, row=0)
        for i in range(self.n_of_chans):
            self.labels.append(Label(top, text='chan' + str(i+1)))
            self.labels[-1].grid(column=i+4, row=0)
        self.cells = []
        for j in range(self.n_of_rows):
            for i in range(self.n_of_chans+4):
                self.cells.append(Entry(top, width=10))
                self.cells[-1].grid(row=j+1, column=i, sticky=tkinter.W)
        i = 0
        # logger.debug('devs map:' + str(fds.DEVS_CHANS_MAP))
        chan_rows = devmap_to_rows(fds.DEVS_CHANS_MAP)
        # logger.debug('chan_rows:' + str(chan_rows))
        for row in chan_rows:
            for row_val in row:
                self.cells[i].insert(0, row_val)
                i += 1

        self.b = Button(top, text="OK", command=self.ok)
        self.b.grid(column=self.n_of_chans+3, row=self.n_of_rows+2)

        return

    def ok(self):
        chan_rows = []
        for n_row in range(self.n_of_rows):
            chan_row = []
            serial = self.cells[n_row * (self.n_of_chans+4)].get()
            if not serial:
                break
            chan_row.append(serial)
            for n_chan in range(self.n_of_chans+3):
                cell = self.cells[n_row * (self.n_of_chans+4) + 1 + n_chan]
                chan_row.append(cell.get())
            logger.debug('chan_row:' + str(chan_row))
            if not chan_row[0]:
                break
            chan_rows.append(chan_row)
        save_chans_csv(chan_rows, chans_csv_file)
        fds.DEVS_CHANS_MAP = rows_to_devmap(chan_rows)
        self.top.destroy()


root = Tk()
root.title('NDAS converter 5.89')

pv_files = IntVar()
pb_files = Progressbar(root, variable=pv_files, maximum=100)
pb_files.grid(columnspan=5, sticky='EW')

progress_var = IntVar()
progressbar = Progressbar(root, variable=progress_var, maximum=100)
progressbar.grid(row=1, columnspan=5, sticky='EW')

status_var = tkinter.StringVar()
status_label = Label(textvariable=status_var, width=30).grid(row=2, column=0, columnspan=2, sticky='EW')
status_var.set('')

txt_var = tkinter.StringVar()

exit_var = tkinter.IntVar()
button = Button(root, textvariable=txt_var, command=lambda: exit_var.set(1))
button.grid(row=2, column=3)

button_back = Button(root, text="<< Back", command=lambda: exit_var.set(2))
button_back.grid(row=2, column=2)

# root.pack()
root.resizable(width=False, height=False)

while exit_var.get() != 1:
    txt_var.set('Processing')
    button_back['state'] = 'disable'
    button['state'] = 'disable'
    root.wm_state('iconic')
    exit_var.set(0)
    progress_var.set(0)
    pv_files.set(0)

    ok_pressed = -1
    while ok_pressed != 1:
        # ok_pressed = -1
        logger.debug('ok_pressed=' + str(ok_pressed))
        d = MyDialog(root)
        logger.debug('ok_pressed=' + str(ok_pressed))
        root.wait_window(d.top)
        logger.debug('ok_pressed=' + str(ok_pressed))
        if ok_pressed == -1:
            break
            # exit(1)
        if not ok_pressed:
            d = ChannelsDialog(root)
            root.wait_window(d.top)
            logger.debug('devs chans map:' + str(fds.DEVS_CHANS_MAP))
    if ok_pressed == -1:
        break

    root.wm_state('normal')

    logger.debug('dir_in: ' + dd.dir_in)

    pattern = re.compile(r'[.]\d{2}$')
    # files = filter(lambda x:x.endswith('.siv') or pattern.search(x) or x.endswith('.sivy'),
    #               [dd.dir_in+'/'+x for x in os.listdir(dd.dir_in)])

    files = []
    if dd.use_file:
        files = [dd.dir_in]
    else:
        files = [dd.dir_in + '/' + x for x in os.listdir(dd.dir_in)]
    files = [x for x in files
             if x.endswith('.siv') or x.endswith('.sivy') or pattern.search(x)]
    if not files:
        logger.error('no appropriate files selected')

    for siv_file in files:
        logger.debug('siv file: ' + siv_file)

    pickleFiles = get_dir_files('pickle', dd.dir_out)
    for pickleFile in pickleFiles:
        os.remove(pickleFile)

    # st = Stream()
    # st_merged = Stream()
    hFile = None
    user_dic = {}
    ret_val = Stream()
    j = 0
    j_max = len(files)
    stations = set()
    os.chdir(dd.dir_out)
    # st_head = Stream()
    st_buff = Stream()
    # ds = None
    hp = None
    hFiles = []
    sth = Stream()
    for siv_file in files:
        i = 1
        siv_file_size = os.path.getsize(siv_file)
        logger.debug('add traces from file:' + siv_file)
        srate_ref = None
        with open(siv_file, "rb") as fp:
            # * Read sivy file frame.
            # * If st_buff is empty, not adjacent or corresponds to different station, start new h5 file.
            # * Write wf portion to h5 file.
            # * Write auxiliary data to file.
            # * Save current portion to st_buff
            ret_val = read_frame(fp, i, rotation, status_var)
            while ret_val:
                if st_buff:
                    adjacent = abs(ret_val[0].stats.starttime - st_buff[-1].stats.endtime) <= \
                        st_buff[-1].stats.delta * 1.5
                if not st_buff or not adjacent or \
                        st_buff[0].stats.station != ret_val[0].stats.station or \
                        (st_buff and st_buff[0].stats['user_stats']['chans'] !=
                         ret_val[0].stats['user_stats']['chans']):
                    if sth:
                        #logger.debug('sth:' + str(sth))
                        pickle.dump(sth, hp, protocol=4)
                        sth = Stream()
                    hFile = get_ofile_base0(ret_val) + 'r.pickle'
                    if hp:
                        hp.close()
                        hp = None
                    hp = open(hFile, 'ab')
                    logger.debug('new file creation:' + hFile)
                    hFiles.append(hFile)
                if ret_val and st_buff and adjacent and \
                        ret_val[0].stats.station == st_buff[0].stats.station and \
                        ret_val[0].stats.sampling_rate != st_buff[0].stats.sampling_rate:
                    for tr in ret_val:
                        tr.stats.sampling_rate = st_buff[0].stats.sampling_rate
                if sth and ret_val[0].stats.starttime.hour != sth[0].stats.starttime.hour:
                    logger.debug('sth reaches limit len, dump sth:' + str(sth))
                    pickle.dump(sth, hp, protocol=4)
                    hp.close()
                    hp = open(hFile, 'ab')
                    sth = Stream()
                sth += ret_val
                st_buff = ret_val
                ret_val = read_frame(fp, i, rotation, status_var)
                progress_var.set(int(fp.tell() / siv_file_size * 100))
                root.update()

        stations.add(st_buff[0].stats.station)
        j += 1
        pv_files.set(int(j / j_max * 30))
        root.update()
    if sth:
        logger.debug('final dump sth:' + str(sth))
        hp.tell()
        pickle.dump(sth, hp, protocol=4)
        hp.close()


    def process_station(hpFile, dd, root, pv):
        st_head = get_st_head(hpFile, None, root, pv)

        if dd.cb:
            cf_data = read_cf(dd.cfile)
            logger.debug('hpFile:' + hpFile + '\nfile size:' + str(os.path.getsize(hpFile)) +
                         '\nst_head:' + str(st_head))
            st_head_chan = st_head.select(channel=st_head[0].stats.channel)
            # for tr in st_head:
            #     if tr.stats.channel != st_head[0].stats.channel:
            #         break
            #     st_head_chan += tr
            tr_head_it = iter(st_head_chan[:-1])
            # tr_len = st_head_chan[1].stats.starttime - st_head_chan[0].stats.starttime
            times_it = iter([control_row.date_time for control_row in cf_data if control_row])
            try:
                cf_time = next(times_it)
                tr_head = next(tr_head_it)
            except StopIteration:
                return
            st_head_ = Stream()
            while True:
                try:
                    starttime = tr_head.stats.starttime
                    if starttime < cf_time:
                        if cf_time < starttime + 60:  # tr_len:
                            st_head_ += tr_head
                        tr_head = next(tr_head_it)
                    else:
                        cf_time = next(times_it)
                except StopIteration:
                    break
            if not st_head_:
                return
            if len(st_head_) > 1:
                st_head_[-1].stats.user_stats[fds.TABLE_SAMPLE_RATE] = \
                    st_head_[-2].stats.user_stats[fds.TABLE_SAMPLE_RATE]
            if dd.cb_bin:
                fp = open(get_ofile_base0(st_head_) + '.csv', 'w', newline='')
                write_spread_sheet(fp, st_head_)
                fp.close()
            cf_data_slice = []
            for control_row in cf_data:
                if control_row and control_row.date_time < st_head[0].stats.starttime:
                    continue
                if control_row and control_row.date_time > st_head[-1].stats.starttime:
                    continue  # break
                cf_data_slice.append(control_row)
            cf_len = len(cf_data_slice)
            max_len = 2 ** 31 - 1
            # if len(cf_data_slice) > max_len:
            #     error_mes = 'The number of included lines ' + str(cf_len) + ' is too high.\n' \
            #                   'Possible maximum is ' + str(max_len) + '.'
            #     logger.error(error_mes)
            #     messagebox.showerror('Error', error_mes)
            shot_sts = cfile_gen(hpFile, cf_data_slice, dd.shot_len)
            st = Stream()
            init_head = True
            len_secs = st_head[-1].stats.starttime - st_head[0].stats.starttime
            chans = st_head[0].stats['user_stats']['chans'].split(' ')
            lastTime = UTCDateTime(0)
            for st in shot_sts:
                if not st:
                    logger.info('next pass')
                    init_head = True
                    continue
                if st[0].stats.starttime > lastTime + 60 * 10:
                    logger.info('next pass')
                    init_head = True
                lastTime = st[-1].stats.endtime
                if len(st[0].data) > 0:
                    if init_head:
                        ofile_base = get_ofile_base0(st)
                        ofile_base_mseed = get_ofile_base0(st, mseed=True)
                    for chan in chans:
                        st_selected = st.select(channel=chan)
                        # logger.debug('st:' + str(st) + '\nst_selected:' + str(st_selected) +
                        #              '\nchan:' + str(chan))
                        if dd.cb_segy:
                            ofile_segy = ofile_base + '_' + chan + '.segy'
                            if init_head:
                                create_segy(st_selected, ofile_segy,
                                            calc_n_of_tr(st_head, cf_data, dd.shot_len))
                            append_to_segy(st_selected, ofile_segy)
                        if dd.cb_bin:
                            ch_num = get_ch_num(chan, st[0].stats['user_stats']['serial'])
                            ofile = get_ofile_base0(st) + '_ch' + str(ch_num) + '.bin'
                            if init_head and os.path.exists(ofile):
                                os.remove(ofile)
                            fp = open(ofile, 'wb')
                            for tr in st_selected:
                                np.frombuffer(tr.data.newbyteorder('<'), 'uint8').tofile(fp)
                            fp.close()
                        if dd.cb_sac:
                            ofile = get_ofile_base0(st_selected) + '_' + chan + '.sac'
                            # st_selected.write(ofile)
                            st_selected.merge().write(ofile)
                        if dd.cb_mseed:
                            st_selected.merge()
                            st_selected.write('temp.mseed')
                            ofile_mseed = ofile_base_mseed + '_' + chan + '.mseed'
                            if init_head:
                                st_selected.write(ofile_mseed)
                            else:
                                fo = open(ofile_mseed, 'ab')
                                fi = open('temp.mseed', 'rb')
                                fo.write(fi.read())
                                fo.close()
                                fi.close()
                    if dd.cb_txt:
                        ofile = get_ofile_base0(st) + '.txt'
                        fp = open(ofile, 'w')
                        fp.write('\t'.join(chans) + '\n')
                        fp.close()
                        save_to_txt_core(st, chans, ofile)
                init_head = False
                if root:
                    pv.set((st[0].stats.starttime - st_head[0].stats.starttime) * 100 / len_secs)
                    root.update()
            if os.path.exists('temp.mseed'):
                os.remove('temp.mseed')
        else:
            gen_sts = hFile_gen(hpFile)
            chans = st_head[0].stats['user_stats']['chans'].split(' ')
            init_head = True

            starttime = st_head[0].stats.starttime
            len_secs = st_head[-1].stats.starttime - starttime

            if dd.cb_bin:
                for chan in chans:
                    ofile = get_ofile_base0(st_head) + '_ch' + chan + '.bin'
                    if os.path.exists(ofile):
                        os.remove(ofile)  # remove previously generated files
                fp = open(get_ofile_base0(st_head) + '.csv', 'w', newline='')
                write_spread_sheet(fp, st_head)
                fp.close()

            for st in gen_sts:
                # for chan in chans:
                #     st_selected = st.select(channel=chan)
                #     logger.debug('st_selected for debug:\n' +
                #                  (st_selected[:8] + st_selected[-8:]).__str__(extended=True))
                for chan in chans:
                    st_selected = st.select(channel=chan)
                    # logger.debug('st_selected:' + st_selected[:8].__str__(extended=True))
                    # logger.debug('chan:' + str(chan) + ' chan from st:' + str(st[0].stats.channel))
                    if dd.cb_segy:
                        ofile = get_ofile_base0(st_head) + '_' + chan + '.segy'
                        if init_head:
                            k_tr = 1
                            max_len = 2 ** 15 - 1
                            number_of_traces = len(st_head.select(channel=chan))
                            if st_selected[0].stats.npts > max_len:
                                k_tr = st_selected[0].stats.npts // max_len + 1
                                logger.debug('st_head:' + str(st_head))
                                logger.debug('number of traces:' + str(number_of_traces))
                            create_segy(st_selected, ofile, k_tr * number_of_traces)
                        append_to_segy(st_selected, ofile)
                    if dd.cb_bin:
                        ofile = get_ofile_base0(st_head) + '_ch' + \
                                str(get_ch_num(chan, st_head[0].stats['user_stats']['serial'])) + '.bin'
                        if init_head and os.path.exists(ofile):
                            os.remove(ofile)
                        fp = open(ofile, 'ab')
                        bytes_data = b''
                        logger.debug('st_selected:' + str(st_selected))
                        for tr in st_selected:
                            bytes_data += np.frombuffer(tr.data.newbyteorder('<'), 'uint8').tobytes()
                        logger.debug('bytes len:' + str(len(bytes_data)))
                        fp.write(bytes_data)
                        fp.close()
                    if dd.cb_sac:
                        ofile = get_ofile_base0(st_selected) + '_' + chan + '.sac'
                        st_selected.merge().write(ofile)
                    if dd.cb_mseed:
                        ofile = get_ofile_base0(st_head, mseed=True) + '_' + chan + '.mseed'
                        # st_selected.merge()
                        st_selected.write('temp.mseed')
                        if init_head and os.path.exists(ofile):
                            os.remove(ofile)
                        fo = open(ofile, 'ab')
                        fi = open('temp.mseed', 'rb')
                        fo.write(fi.read())
                        fo.close()
                        fi.close()
                if dd.cb_txt:
                    ofile = get_ofile_base0(st_head) + '.txt'
                    if init_head:
                        fpt = open(ofile, 'w')
                        fpt.write('\t'.join(chans) + '\n')
                        fpt.close()
                    save_to_txt_core(st, chans, ofile)
                init_head = False
                if len_secs > 0:
                    pv.set((st[0].stats.starttime - starttime) * 100 / len_secs)
                    if root:
                        root.update()

            if os.path.exists('temp.mseed'):
                os.remove('temp.mseed')
        if os.path.exists('temp.segy'):
            os.remove('temp.segy')


    log_file = dd.dir_in + '/LOG.TXT'
    drift_map = {}
    if os.path.exists(log_file):
        drift_map = get_log_map(log_file)
    pattern = re.compile(r'^LOG_\d\d\d.txt$')
    if not dd.use_file:
        log_files = [dd.dir_in + '/' + x for x in os.listdir(dd.dir_in) if pattern.match(x)]
        for log_file in log_files:
            drift_map.update(get_new_map(log_file))
    if gps_pt:
        drift_map.update(gps_pt)
    sts_head = [get_st_head(hFile, None, root, progress_var) for hFile in hFiles]
    sts_head_chan = [st_head.select(channel=st_head[0].stats.channel) for st_head in sts_head]
    sts_group = reduce(fre, sts_head_chan[1:], [sts_head_chan[0]])
    sts_group = [tune_drift(st_group, root, progress_var, drift_map) for st_group in sts_group]
    sts_head_chan = []
    for st_group in sts_group:
        st_head_chan = Stream() + st_group[0]
        for tr_prev, tr in zip(st_group, st_group[1:]):
            if tr.stats.starttime - tr_prev.stats.starttime > 100:
                sts_head_chan.append(st_head_chan)
                st_head_chan = Stream() + tr
            else:
                st_head_chan += tr
        sts_head_chan.append(st_head_chan)
    for st_head, st_head_chan in zip(sts_head, sts_head_chan):
        chans = st_head[0].stats['user_stats']['chans'].split(' ')
        for chan in chans[1:]:
            st_head_raw = st_head.select(channel=chan)
            for tr_head_raw, tr_head_tuned in zip(st_head_raw, st_head_chan):
                tr_head_raw.stats['user_stats'][fds.TABLE_TUNED_DRIFT] = \
                    tr_head_tuned.stats['user_stats'][fds.TABLE_TUNED_DRIFT]
                tr_head_raw.stats['user_stats'][fds.TABLE_TEMPER_DRIFT] = \
                    tr_head_tuned.stats['user_stats'][fds.TABLE_TEMPER_DRIFT]
                tr_head_raw.stats['user_stats']['tuned_time'] = \
                    tr_head_tuned.stats.starttime + \
                    tr_head_tuned.stats['user_stats'][fds.TABLE_TUNED_DRIFT] / 10 ** 9
        if dd.corr:
            for tr in st_head:
                if 'tuned_time' in tr.stats.user_stats.keys():
                    tr.stats.starttime = tr.stats.user_stats['tuned_time']
    for hFile, st_head in zip(hFiles, sts_head):
        chans = st_head[0].stats['user_stats']['chans'].split(' ')
        st_new = Stream()
        htFile = os.path.splitext(hFile)[0][:-1] + 't.pickle'
        hp = open(htFile, 'ab')
        tr_gen = TrGen(hFile)
        for tr_head, tr in zip(st_head, tr_gen):
            if tr_head.stats.channel != tr.stats.channel:
                logger.error('channel inconsistency, tr_head:' + str(tr_head) + ' tr:' + str(tr))
                exit(1)
            tr.stats['user_stats'] = tr_head.stats['user_stats']
            if dd.corr:
                tr.stats.starttime = tr_head.stats.starttime
            st_new += tr
            if st_new[-1].stats.starttime.hour != st_new[0].stats.starttime.hour and \
                    len(st_new) % st_new[0].stats['user_stats']['n_of_chs'] == 0:
                logger.debug('dump st_new:' + str(st_new))
                pickle.dump(st_new, hp, protocol=4)
                progress_var.set((tr_head.stats.starttime - st_head[0].stats.starttime) * 100 /
                                 (st_head[-1].stats.starttime - st_head[0].stats.starttime))
                root.update()

                st_new = Stream()
        if st_new:
            logger.debug('dump st_new:' + str(st_new))
            pickle.dump(st_new, hp, protocol=4)
        hp.close()
        os.remove(hFile)
        os.rename(htFile, hFile)

    for hFile in hFiles:
        progress_step = 30 // len(hFiles)
        if dd.corr:
            interpolate_h(hFile, root=root, pv=progress_var)
        pv_files.set(pv_files.get() + progress_step)
        root.update()
        if dd.corr:
            hpFile = os.path.splitext(hFile)[0][:-1] + 'p.pickle'
        else:
            hpFile = hFile
        process_station(hpFile, dd, root, progress_var)
        if hpFile != hFile:
            os.remove(hpFile)
        os.remove(hFile)
        # st_head_file = os.path.splitext(hFile)[0][:-1] + '.pickle'
        # logger.debug('st_head_file:' + str(st_head_file) + ' work dir:' + str(os.getcwd()))
        # if os.path.exists(st_head_file):
        #     os.remove(st_head_file)
        pv_files.set(pv_files.get() + progress_step)
        root.update()

    pv_files.set(100)
    progress_var.set(100)
    button['state'] = 'normal'
    txt_var.set('Ready')
    button_back['state'] = 'normal'
    root.wait_variable(exit_var)

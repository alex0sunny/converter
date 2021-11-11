'''
Created on 16 дек. 2017 г.

@author: Alexey
'''
import numpy as np
import logging

from conv_util import format_exception

logger = logging.getLogger(__name__)


def read_frame_data(fp, n_of_bytes, n_of_chs, i24=True):
    try:
        fp.tell()
        buf = fp.read(n_of_bytes)
        ar = np.frombuffer(buf, 'uint8')
    except OSError as e:
        logger.error('er no:' + str(e.errno))
        logger.error('os error:' + format_exception(e))
        return np.array([])
    
    if i24:
        # reshape array
        if ar.size % 3:
            logger.warning('The size of frame body ' + str(ar.size) + ' is not multiple of 3.')
            return np.array([])
        ar = ar.reshape(-1, 3).transpose()
        #logger.debug('data chunked:' + str(ar))
        empty_data = np.zeros([ar.shape[1]], dtype='uint8').reshape(1, -1)
        ar = np.append(empty_data, ar, axis=0).transpose()
        #logger.debug('zero row appended:' + str(ar))
    
        # now convert back to byte array
        ar = ar.reshape(-1)
        #logger.debug('bytes finally prepared:' + str(ar))
    ar = np.frombuffer(ar, '<i')
    # ar = ar // 256 no need since we use maximum integer 32
    
    # demultiplex
    if ar.size % n_of_chs:
        logger.warning('The size of frame body ' + str(ar.size) +
                       ' is not multiple of n_of_chs:' + str(n_of_chs) + '.')
        return np.array([])
    ar = ar.reshape(-1, n_of_chs).transpose()
    #logger.debug('multiplexed:' + str(ar))
    
    return ar


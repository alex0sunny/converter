import logging

from numpy import dot, cross, array, array_equal, arccos, degrees
from numpy.linalg import norm
from obspy import Stream
from scipy.spatial.transform import Rotation

import fields_names as fds

logger = logging.getLogger(__name__)


def normalize(vec):
    return vec / norm(vec)


def calc_matrix(a_vec, m_vec):
    e_v = - normalize(a_vec)
    m_h = m_vec - dot(m_vec, e_v) * e_v
    e_n = - normalize(m_h)
    e_e = cross(e_v, e_n)
    return [e_n, e_e, e_v]


def check_channels(stream: Stream) -> bool:
    required = ['n', 'e', 'z']
    channels = [tr.stats.channel for tr in stream]
    short_names = [channel[-1].lower() for channel in channels]
    if not all([short_name in short_names for short_name in required]):
        logger.warning(f'no rotation since not all {required} channels are present in {channels}')
        return False
    if len(short_names) != 4:
        logger.warning('no rotation is applied since num of channels is ' +
                       str(len(channels)) + ' while 4 is expected')
        return False
    return True


def rotate_core(stream: Stream, m) -> Stream:
    data_debug = [tr.data[0] for tr in stream]
    logger.debug('debug data:' + str(data_debug))
    required = ['n', 'e', 'z']
    channels = [tr.stats.channel for tr in stream]
    short_names = [channel[-1].lower() for channel in channels]
    channel_inds = [short_names.index(short_name) for short_name in required]
    # logger.debug('before transpose:' +
    #              str(array([stream[channel_ind].data for channel_ind in channel_inds])))
    vectors = array([stream[channel_ind].data for channel_ind in channel_inds]).transpose()
    # logger.debug('vectors:' + str(vectors))
    r = Rotation.from_matrix(m)
    vectors = r.apply(vectors)
    #  logger.debug('vectors:' + str(vectors))
    chs_data = vectors.transpose()
    # logger.debug('before transpose:' +
    #              str(array([stream[channel_ind].data for channel_ind in channel_inds])))
    for i, channel_ind in enumerate(channel_inds):
        stream[channel_ind].data = chs_data[i].astype('float32')
    data_debug = [tr.data[0] for tr in stream]
    logger.debug('debug data:' + str(data_debug))
    #exit(1)
    return stream


def rotate(stream: Stream) -> Stream:
    if not check_channels(stream):
        return stream
    user_stats = stream[0].stats['user_stats']
    a_vec = [user_stats[fds.SIVY_ACC_X], user_stats[fds.SIVY_ACC_Y], user_stats[fds.SIVY_ACC_Z]]
    m_vec = [user_stats[fds.SIVY_COMP_X], user_stats[fds.SIVY_COMP_Y], user_stats[fds.SIVY_COMP_Z]]
    if any([array_equal(vec, [0, 0, 0]) for vec in [a_vec, m_vec]]):
        logger.warning(f'vectors of accelerometer {a_vec} and magnetometer {m_vec} ' +
                       'should not be zero, no rotation')
        return stream
    m = calc_matrix(a_vec, m_vec)
    e_v = m[-1]
    tilt = degrees(arccos(-e_v[-1]))
    e_n = m[0]
    azimuth = degrees(arccos(e_n[0]))
    if e_n[1] < 0:
        azimuth = 360-azimuth
    for tr in stream:
        tr.stats['user_stats']['tilt'] = tilt
        tr.stats['user_stats']['azimuth'] = azimuth
    return rotate_core(stream, m)


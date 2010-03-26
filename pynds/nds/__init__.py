"""

Wrapper for John Zweizig's NDS1/NDS2 library

"""
__author__       = "Leo Singer <leo.singer@ligo.org>"
__organization__ = ["LIGO", "California Institute of Technology"]
__copyright__    = "Copyright 2010, Leo Singer"

import nds_ext
from nds_ext import channel_type, nds_version

class daq_iterator:
    def __init__(self, daq, start, stop, stride):
        self.__daq = daq
        self.__start = start
        self.__stop = stop
        self.__stride = stride
        self.__started = False
    def __iter__(self):
        return self
    def next(self):
        if self.__started:
            self.__daq.recv_next()
        else:
            self.__daq.request_data(self.__start, self.__stop, self.__stride)
            self.__started = True
        return self.__daq.unpack()

class daq(nds_ext.daq):
    def __init__(self, host, port, version=None):
        if version is None:
            super(daq, self).__init__(host, port)
        else:
            super(daq, self).__init__(host, port, version)
        self.host = host
        self.port = port
    def seek(self, start, stop, stride):
        return daq_iterator(self, start, stop, stride)
    def request_channels(self, channels):
        for channel in channels:
            self.request_channel(channel)

del nds_ext

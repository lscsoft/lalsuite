import _nds2

class Daq:
    daq           = None
    channels      = []
    channel_map   = {}
    have_channels = False

    def __init__(self, host, port):
       self.daq = _nds2.daq_t()

       _nds2.daq_startup()
       _nds2.daq_connect(self.daq, host, port) # 'ldas-pcdev1.ligo.caltech.edu',31200)


    def get_channel_list(self):
        if not self.have_channels:
            self.channels = [_nds2.daq_channel_t() for i in range(65536)]
            _nds2.daq_recv_channels(self.daq, self.channels)
            
            self.channel_map = {}
            for c in self.channels:
                self.channel_map[c.name] = c

            self.have_channels = True

        return self.channels


    def get_data(self, channels, gps_start, gps_end):
        return _nds2.get_data(self.daq, channels, gps_start, gps_end)

    def get_data_old(self, channel, gps_start, gps_end):
        self.get_channel_list()

        _nds2.daq_request_channel_from_chanlist(self.daq, self.channel_map[channel])

        data = _nds2.daq_request_data(self.daq, gps_start, gps_end, gps_end - gps_start)

        return data


    def close(self):
        if self.daq:
            _nds2.daq_disconnect(self.daq)
            _nds2.daq_recv_shutdown(self.daq)
        self.daq = None

    def __del__(self):
        self.close()


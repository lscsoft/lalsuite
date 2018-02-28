#!/usr/bin/python
"""Assists pylal_mvsc_player in writing html files."""

__author__ = 'Tristan Miller <tmiller@caltech.edu>'
__prog__ = 'mvsc_htmlwriter'

##############################################################################
# import modules
import pylab,os,re

##############################################################################

class mvsc_html:
    """Contains all necessary data to be put in html file, including links to
    plots and thumbnails, and info.

    The most essential function is write_html, which writes an html file from
    a template, inserting the info and images in the right places."""

    def __init__(self,opts,fname,templatepath):
        self.opts = opts
        self.filename = fname
        self.htmlpath = opts.output_path + '/'
        self.templatepath = templatepath

        #I would make datapath and patpath relative paths rather than
        #absolute paths, but I don't know how.
        self.datapath = os.path.abspath(opts.data_path)+'/'
        self.patpath = os.path.abspath(opts.pat_path)+'/'
        
        self.letters = None
        self.figs = {}
        self.gpstime = None
        self.top_events = None
        self.op_point = None
        self.trigger_info = None
        self.treesplits = None
        self.cols = None
        self.zeronum = None
        self.efficiency = None
        self.comments = ''

    #functions to set parameters
    
    def add_figure(self,type,fig=None,dpi=None,dpi_thumb=50):
        """Adds figure of given type to be displayed in html file.

        Figure types: 'FOM','ROC','Frac_vs_SNR','MVSC_vs_effSNR','missed_inj'
        and possibly more to be added"""

        if fig is None:
            fig = pylab.gcf()
        if dpi is None:
            dpi = pylab.rcParams["savefig.dpi"]
            
        #save picture into file
        figpath = self.htmlpath+'Images/'
        figname = self.filename+'_'+type+'.png'
        figname_th = self.filename+'_'+type+'_th.png'
        fig.savefig(figpath+figname,dpi=dpi)
        fig.savefig(figpath+figname_th, dpi=dpi_thumb)
        
        self.figs[type] = [figname,figname_th]

    def set_gpstime(self,gpstime):
        self.gpstime = gpstime
    
    def set_top_events(self,events,cols):
        self.top_events = events
        self.cols = cols

    def set_op_point(self,afar,mvsc_cutoff,effsnr_cutoff):
        self.op_point = [afar,mvsc_cutoff,effsnr_cutoff]

    def set_trigger_info(self,ts_tr,inj_tr,ts_va,inj_va):
        self.trigger_info = [ts_tr,inj_tr,ts_va,inj_va]

    def set_treesplits(self,treesplits):
        self.treesplits = treesplits

    def set_data_sets(self,patletters):
        self.letters = patletters

    def set_zeronum(self,zeronum):
        self.zeronum = zeronum

    def set_efficiency(self,lower,upper):
        self.efficiency = [lower,upper]

    def set_comments(self,comments):
        self.comments = comments

##############################################################################
    
    def write_file(self):
        """Writes an html file putting all information in a template."""
        
        #open template
        tempfile = open(self.templatepath)
        template = tempfile.read()
        tempfile.close()

        #Looks for the following pattern in the template:
        p = re.compile(r'\[InsertCode:\s+(\S+)\s+AltText:\s+([^\]]*)\]')

        #Replaces pattern with appropriate figure or info
        match = p.search(template)
        while match:
            #First, try to find matching InsertCode
            #If the necessary information exists, insert it
            #If not, insert the AltText

            # insert figures
            if match.group(1)[:4] == 'fig_':
                if self.figs.has_key(match.group(1)[4:]):
                    figpaths = self.figs[match.group(1)[4:]]
                    htmlcode = '<a href="Images/' + figpaths[0] + \
                              '"><img src="Images/' + figpaths[1] + '"></a>'
                else:
                    htmlcode = match.group(2)

            # insert gps time
            elif match.group(1) == 'gpstime':
                if self.gpstime:
                    htmlcode = self.gpstime
                else:
                    htmlcode = match.group(2)

            # insert options
            elif match.group(1) == 'opts':
                try:
                    if self.opts.balance_data:
                        bstr = ' -b'
                    else:
                        bstr = ''

                    if self.opts.seed:
                        rstr = ' -R ' + self.opts.seed
                    else:
                        rstr = ''
                    
                    astr = ' -a '
                    for i in range(2):
                        astr += self.letters[i]

                    if self.opts.open_box:
                        zstr = ' --open-box'
                    elif self.opts.zero_lag:
                        zstr = ' --zero-lag'
                    elif self.opts.hardware:
                        zstr = ' --hardware'
                    else:
                        zstr = ''

                    if self.opts.s:
                        sstr = ' -s ' + self.opts.s
                    else:
                        sstr = ''

                    htmlcode = '-n '+self.opts.n+' -l '+self.opts.l+' -c '+ \
                               self.opts.c+' -g '+self.opts.g + sstr + bstr+ \
                               rstr + astr + zstr
                except:
                    htmlcode = match.group(2)

            # insert file links
            elif match.group(1) == 'files':
                try:
                    training_path = self.patpath+ \
                       self.opts.run_name+\
                       '/'+self.opts.stations+ 'set'+ self.letters[0] + \
                       'Known.pat'
                    validation_path = self.patpath+ \
                       self.opts.run_name+\
                       '/'+self.opts.stations+ 'set'+ self.letters[1] + \
                       'Known.pat'
                    testing_path = self.patpath+ \
                       self.opts.run_name+\
                       '/'+self.opts.stations+ 'sets'+ self.letters[2] + \
                       '.pat'
                    
                    tree_path = self.datapath+self.filename
                    test_path = tree_path + '_test.dat'
                    info_path = tree_path + '_info'
                    tree_path += '.spr'
                    
                    htmlcode = '<a href="'+training_path+'">Training Set</a>'+\
                        ', <a href="'+validation_path+'">Validation Set</a>'+\
                        ', <a href="'+testing_path+'">Testing Set</a><br>'+\
                        '<a href="'+info_path+'">SprBaggerDecisionTreeApp '+\
                        'printout</a>, <a href="'+tree_path+ \
                        '">Generated trees</a>, <a href="'+test_path+ \
                        '">Test results</a>'
                    
                    if self.opts.open_box | self.opts.zero_lag:
                        zeropath = self.patpath+ \
                            self.opts.run_name+\
                            '/'+self.opts.stations+ 'setZeroLag_'

                        zerotest = self.datapath+self.filename+'_'
                        
                        if self.opts.open_box:
                            zeropath += 'fulldata.pat'
                            zerotest += 'fulldata.dat'
                            zstr = 'full data'
                        else:
                            zeropath += 'playground.pat'
                            zerotest += 'playground.dat'
                            zstr = 'playground'
                            
                        htmlcode += '<br><a href="'+zeropath+ \
                            '">Zero Lag Set ('+zstr+')</a>'
                        htmlcode += ', <a href="'+zerotest+ \
                            '">Zero Lag test results</a>'
                        
                except:
                    htmlcode = match.group(2)

            # insert trigger info
            elif match.group(1) == 'trigger_info':
                if self.trigger_info:
                    htmlcode = 'Number of timeslides in training set: ' + \
                        str(self.trigger_info[0]) + \
                        '<br>Number of injections in training set: ' + \
                        str(self.trigger_info[1]) + \
                        '<br>Number of timeslides in validation set: ' + \
                        str(self.trigger_info[2]) + \
                        '<br>Number of injections in validation set: ' + \
                        str(self.trigger_info[3])
                    
                    if self.zeronum:
                        htmlcode += '<br>Number of zero lag triggers: ' + \
                                    str(self.zeronum)
                else:
                    htmlcode = match.group(2)

            # insert operating point info
            elif match.group(1) == 'op_point':
                if self.op_point:
                    htmlcode = 'False alarm rate per trigger: ' + \
                        str(self.op_point[0]) + '<br>MVSC cutoff value: ' + \
                        str(self.op_point[1]) + \
                        '<br>Combined Effective SNR squared cutoff value: ' + \
                        str(self.op_point[2])

                    if self.efficiency:
                        htmlcode += '<br>Resulting efficiency: ' + \
                            str(self.efficiency[0]) + ' < p < ' + \
                            str(self.efficiency[1]) + ' (at 68% CL)'
                else:
                    htmlcode = match.group(2)

            # insert tree splits info
            elif match.group(1) == 'treesplits':
                if self.treesplits:
                    htmlcode = ''
                    for i in range(len(self.treesplits)):
                        htmlcode += '<tr>'
                        for j in range(3):
                            htmlcode += '<td>'+str(self.treesplits[i][j])+\
                                        '</td>'
                        htmlcode += '</tr>'
                            
                else:
                    htmlcode = match.group(2)

            # insert list of top events
            elif match.group(1) == 'top_events':
              if self.top_events:
                htmlcode = ''
                for i in range(len(self.top_events)):
                  htmlcode += '<tr><td>'+str(i+1)+'</td>'

                  colinfo = ['get_end()','snr','chisq','eff_distance']
                  for j in range(len(colinfo)):
                    htmlcode += '<td>'
                    htmlcode += str(self.top_events[i][ \
                                self.cols['t1'+colinfo[j]]])
                                            
                    k = 2
                    while self.cols.has_key('t'+str(k)+colinfo[j]):
                      htmlcode += '<br>'
                      htmlcode += str(self.top_events[i][ \
                                    self.cols['t'+str(k)+colinfo[j]]])
                      k += 1

                    htmlcode += '</td>'
                            
                  htmlcode += '<td>'+\
                        str(self.top_events[i][self.cols['t1t2mchirp']]) +\
                        '</td><td>' + \
                        str(self.top_events[i][self.cols['FAN']]) + \
                        '</td><td>' + \
                        str(self.top_events[i][self.cols['Bagger']]) + \
                        '</td></tr>'
              else:
                htmlcode = match.group(2)
            elif match.group(1) == 'filename':
                try:
                    htmlcode = 'Html file for ' + self.filename
                except:
                    htmlcode = match.group(2)
            elif match.group(1) == 'comments':
                if self.comments:
                    htmlcode = 'Comments: ' + self.comments
                else:
                    htmlcode = match.group(2)
            else:
                htmlcode = 'Error! "'+match.group(1)+ \
                           '" is an invalid InsertCode!'
            
            template = template[:match.start()]+htmlcode+template[match.end():]

            m = p.findall(template)
            match = p.search(template)

        #Writes the html file
        fname = self.htmlpath + self.filename + '.html'
        os.system('touch ' + fname)

        f = open(fname,'w')
        f.write(template)
        f.close()

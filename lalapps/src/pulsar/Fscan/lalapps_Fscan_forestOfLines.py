from pathlib import Path
import os
import re
import sys
from datetime import datetime
from datetime import timedelta
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.dates import DateFormatter
import numpy
import smtplib


def create_scatter_plot(date_list, value_list, output_path, prob_cutoff):
    fig, ax = matplotlib.pyplot.subplots()
    ax.plot_date(date_list, value_list, 'o')
    ax.fmt_xdata = DateFormatter('%b %d %Y')
    fig.autofmt_xdate()
    titlestr = 'Date vs Line Count Above Cutoff for FAP = {0:.2e}'.format(prob_cutoff)
    ax.set_title(titlestr)
    ax.set_ylabel('Number of Lines Above Cutoff',fontsize = 10)
    matplotlib.pyplot.savefig(output_path)


def create_histogram(value_list, output_path, prob_cutoff):
    fig, ax = matplotlib.pyplot.subplots()
    bin_count = 10
    n, bins, patches = ax.hist(value_list, bin_count)
    titlestr = 'Line Count Above Cutoff for FAP = {0:.2e}'.format(prob_cutoff)
    ax.set_title(titlestr)
    ax.set_xlabel('Number of Lines Above Cutoff',fontsize = 10)
    ax.set_ylabel('Number of Occurrences',fontsize = 10)
    matplotlib.pyplot.savefig(output_path)


def fscanForestOfLines(filename, num_days, prob_cutoff):
    # Reads the current path.
    cur_path = Path.cwd()

    # Returns if the current path does match STRAIN.
    if not re.match('.*STRAIN', cur_path.name):
        sys.stderr.write('not running in the expected path\n')
        return

    # Returns if the current frequenies are not as expected.
    pattern = re.compile('^spec_0.00_100.00_.*')
    if not pattern.match(filename):
        sys.stderr.write('not running on the expected frequencies\n')
        return

    # Store the frequency band:
    fscan_band = filename[5:16]

    # Returns if there is no file in the current path that has a name matching
    # the regular expression. There should be only one
    # matching file in the current path (not checked). It should contain the
    # spectrum data for the 0 - 100 Hz band. The first line in the file is a
    # header. The remainder of the file is a pair of columns (freq, value).
    spectrum_file_name_list = [f for f in os.listdir('.') if re.match(r'^spec_0.00_100.00_(H1|L1)_[0-9]*_[0-9]*.txt$', f)]
    if not spectrum_file_name_list:
        sys.stderr.write('spectrum data file not found in the current path\n')
        return

    # Reads in the second column (value) from the spectrum data file after
    # discarding the header.
    spectrum_file_name = spectrum_file_name_list[0]
    with open(str(spectrum_file_name), 'r') as file_in:
        data = file_in.readlines()
    data = [float(line.split()[1]) for line in data[1:]]

    # len(prelim) is the number of sffts used to create the spectrum.
    prelim = list(cur_path.glob('logs/MakeSFTs*.out'))
    if not prelim:
        sys.stderr.write('prelim not found\n')
        return

    deg_of_freedom = len(prelim) * 2

    # The number of lines is the number of values in the spectrum greater than
    # or equal to the cutoff. The cutoff is determined statistically from the
    # number of sffts used to create the spectrum.
    cutoff = float(stats.chi2.isf(prob_cutoff, deg_of_freedom)) / float(deg_of_freedom)
    current_line_count = len(list(filter(lambda num : num >= cutoff, data)))

    # Stores the line count in a file in the current path.
    # line_count_file_name = 'line_count_0.00_100.00_' + str(prob_cutoff) + '.txt'
    line_count_file_name = 'line_count_0.00_100.00.txt'
    with open(line_count_file_name, 'w') as file_out:
        file_out.write(str(current_line_count))


    # start_date is the date referenced by the current path.
    fscan_path = cur_path.parents[0].name
    fscan_path_yyyy_mm_dd = fscan_path[0:17]
    #start_date = datetime.date(datetime.strptime(cur_path.parents[0].name, 'fscans_%Y_%m_%d_%H_%M_%S_%Z_%a'))
    start_date = datetime.date(datetime.strptime(fscan_path_yyyy_mm_dd, 'fscans_%Y_%m_%d'))

    # date_list is a list containing num_days number of consecutive days prior
    # to and including start_date.
    date_list = list(reversed([start_date - timedelta(days=x) for x in range(num_days)]))

    # Reads the number of line counts for each date in date_list.
    # Dates that do not have a line count file are reported and skipped.
    line_count_date_list = []
    line_count_value_list = []
    for date in date_list:
        # Finds the path for the given date. Can not use the exact path because the time part of it varies.
        path_list = list(cur_path.parents[1].glob(datetime.strftime(date, 'fscans_%Y_%m_%d_*/*STRAIN*')))
        if path_list:
            # path is the path to the file containing the line counts for the given date.
            path = Path(path_list[0], line_count_file_name)
            try:
                with open(str(path), "r") as file_in:
                    line_count = int(file_in.readline())
                    line_count_value_list.append(line_count)
                    line_count_date_list.append(date)
            except Exception as e:
                sys.stderr.write(str(e) + '\n')
        else:
            sys.stderr.write(datetime.strftime(date, 'no path found for %Y_%m_%d.\n'))

    # Creates a scatter plot of the line counts over the last num_days.
    #scatter_plot_file_name = ('scatter_' +
    #datetime.strftime(date_list[0], '%Y_%m_%d') + '_to_' +
    #datetime.strftime(date_list[-1], '%Y_%m_%d') +
    #'_0.00_100.00_' + str(prob_cutoff) + '.png')
    scatter_plot_file_name = filename +'_line_count_scatter_plot.png'
    create_scatter_plot(line_count_date_list, line_count_value_list, scatter_plot_file_name, prob_cutoff)

    # Creates a histogram of the line counts over the last num_days.
    #histogram_file_name = ('histogram_' +
    #datetime.strftime(date_list[0], '%Y_%m_%d') + '_to_' +
    #datetime.strftime(date_list[-1], '%Y_%m_%d') +
    #'_0.00_100.00_' + str(prob_cutoff) + '.png')
    histogram_file_name = filename +'_line_count_histogram.png'
    create_histogram(line_count_value_list, histogram_file_name, prob_cutoff)

    if (len(line_count_date_list) < 7):
        sys.stderr.write('not enough data\n')
        return

    # Sends an email alert if the line counts have changed significantly
    # in comparison to those for the last num_days.
    mean = numpy.mean(line_count_value_list[:-1])
    stdev = numpy.std(line_count_value_list[:-1])
    threshold = mean + (2 * stdev)
    if current_line_count > threshold:
        msg = 'current line count for {0}, channel {1}, and frequency band {2}: {3}\nthreshold: {4}\n'.format(fscan_path_yyyy_mm_dd, cur_path.name, fscan_band, current_line_count, threshold)
        msg = 'Subject: {}\n\n{}'.format('Line count detected above threshold by fscanForestOfLines.py', msg)
        s = smtplib.SMTP('localhost')
        s.sendmail('egoetz@phas.ubc.ca', 'aligo-lines@sympa.ligo.org', msg)
        s.quit()


def main():
    fscanForestOfLines(filename, 30, 1e-6)


if __name__ == "__main__":
    main()


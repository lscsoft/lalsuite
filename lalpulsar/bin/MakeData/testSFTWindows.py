import sys
import csv
import subprocess
import numpy as np

maxd = 2e-3

exitcode = 0

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    plot = True
except:
    plot = False

# generate window data
window_data = subprocess.check_output(['./genSFTWindows']).decode('utf-8').splitlines()

# parse window data from CSV
window_csv = csv.reader(window_data)
window_csv_itr = iter(window_csv)
header = next(window_csv)
win = dict()
for name in header:
    win[name] = list()
for row in window_csv_itr:
    for name, value in zip(header, row):
        win[name].append(float(value))
for name in win:
    win[name] = np.array(win[name])
t = np.arange(len(win[name])) / 256

# compare windows
for k in ('Tukey', 'Hann'):
    win_ref_name, win_ref = [(name, win[name]) for name in win if ('Create' + k) in name][0]
    if plot:
        headtail = 512
        ranges = [range(0, len(t)), range(0, headtail), range(len(t) - headtail, len(t))]
        legend_loc = ['lower center', 'lower right', 'lower left']
        fig = plt.figure(tight_layout=True, figsize=(16, 8))
        gs = gridspec.GridSpec(3, 2)
        axs = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]
        for i, ax in enumerate(axs[0:3]):
            ax.plot(t[ranges[i]], win_ref[ranges[i]], label=win_ref_name, color='black', linewidth=2, zorder=-10)
            ax.set_xticks(np.linspace(min(t[ranges[i]]), max(t[ranges[i]] + t[1]), 5))
            ax.ticklabel_format(axis='x', useOffset=False)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax.set_xlabel('Time / seconds')
            ax.set_ylabel('Window value')
    for n, (win_cmp_name, win_cmp) in enumerate([(name, win[name]) for name in win if (k in name and name != win_ref_name)]):
        if plot:
            for i, ax in enumerate(axs[0:3]):
                ax.plot(t[ranges[i]], win_cmp[ranges[i]], label=win_cmp_name, linewidth=1)
                ax.legend(loc=legend_loc[i])
                ax.ticklabel_format(axis='x', useOffset=False)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        d = win_cmp - win_ref
        md = max(abs(d))
        print(f'max(|{win_cmp_name} - {win_ref_name}|) = {md:0.2e}, max = {maxd:0.2e}')
        if md > maxd:
            exitcode = 1
        if plot:
            for i, ax in enumerate(axs[3:]):
                ax.plot(t[ranges[i + 1]], d[ranges[i + 1]], color='black')
                ax.ticklabel_format(axis='x', useOffset=False)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                ax.set_xticks(np.linspace(min(t[ranges[i + 1]]), max(t[ranges[i + 1]] + t[1]), 5))
                ax.ticklabel_format(axis='x', useOffset=False)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                ax.set_xlabel('Time / seconds')
                ax.set_ylabel('Window difference')
            plt.savefig(f'SFTWindows{k}{n}.png')

sys.exit(exitcode)

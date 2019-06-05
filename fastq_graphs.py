from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter, MultipleLocator
import matplotlib.patches as mpatches
from math import sqrt, ceil
import functions


class Layout(object):
    def __init__(self, structure, template, xticks, yticks):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks


class FastqPlots(object):
    @staticmethod
    def find_best_matrix(n_sample):
        """
        Trying to alway get higher than wider
        :param n_sample: Number of samples to use to create the matrix
        :return:
        """
        width = int(sqrt(n_sample))
        is_even = False
        if n_sample % 2 == 0:
            is_even = True

        if is_even:
            while True:
                if n_sample % width == 0:
                    break
                else:
                    width += 1
            height = int(n_sample / width)
        else:
            height = ceil(n_sample / width)

        return width, height

    @staticmethod
    def make_layout(maxval):
        """Make the physical layout of the MinION flowcell.
        based on https://bioinformatics.stackexchange.com/a/749/681
        returned as a numpy array
        """
        if maxval > 512:
            return Layout(
                structure=np.concatenate([np.array([list(range(10 * i + 1, i * 10 + 11))
                                                    for i in range(25)]) + j
                                          for j in range(0, 3000, 250)],
                                         axis=1),
                template=np.zeros((25, 120)),
                xticks=range(1, 121),
                yticks=range(1, 26))
        else:
            layoutlist = []
            for i, j in zip(
                    [33, 481, 417, 353, 289, 225, 161, 97],
                    [8, 456, 392, 328, 264, 200, 136, 72]):
                for n in range(4):
                    layoutlist.append(list(range(i + n * 8, (i + n * 8) + 8, 1)) +
                                      list(range(j + n * 8, (j + n * 8) - 8, -1)))
            return Layout(
                structure=np.array(layoutlist).transpose(),
                template=np.zeros((16, 32)),
                xticks=range(1, 33),
                yticks=range(1, 17))

    @staticmethod
    def jointplot_w_hue(data, x, y, hue=None, colormap=None,
                        figsize=None, fig=None, scatter_kws=None):
        """
        https://gist.github.com/ruxi/ff0e9255d74a3c187667627214e1f5fa
        :param data: pandas dataframe
        :param x:
        :param y:
        :param hue:
        :param colormap:
        :param figsize:
        :param fig:
        :param scatter_kws:
        :return:
        """

        sns.set_style('darkgrid')

        # defaults
        if colormap is None:
            colormap = sns.color_palette()  # ['blue','orange']
        if figsize is None:
            figsize = (5, 5)
        if fig is None:
            fig = plt.figure(figsize=figsize)
        if scatter_kws is None:
            scatter_kws = dict(alpha=0.4, lw=1)

        # derived variables
        if hue is None:
            return "use normal sns.jointplot"
        hue_groups = data[hue].unique()

        subdata = dict()
        colors = dict()

        active_colormap = colormap[0: len(hue_groups)]
        legend_mapping = []
        for hue_grp, color in zip(hue_groups, active_colormap):
            legend_entry = mpatches.Patch(color=color, label=hue_grp)
            legend_mapping.append(legend_entry)

            subdata[hue_grp] = data[data[hue] == hue_grp]
            colors[hue_grp] = color

        # canvas setup
        grid = gridspec.GridSpec(2, 2,
                                 width_ratios=[4, 1],
                                 height_ratios=[1, 4],
                                 hspace=0, wspace=0)

        ax_main = plt.subplot(grid[1, 0])
        ax_xhist = plt.subplot(grid[0, 0], sharex=ax_main)
        ax_yhist = plt.subplot(grid[1, 1], sharey=ax_main)

        # Set main plot x axis scale to log
        ax_main.set_xscale('log')
        # Set x-axis limits
        min_len = min(data.ix[:, 0])
        max_len = max(data.ix[:, 0])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))
        ax_main.set_xlim((min_value, max_value))

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 30)
        # Set y-axis limits
        min_phred = min(data.ix[:, 1])
        max_phred = max(data.ix[:, 1])
        # phred_range = max_phred - min_phred
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Set the y limits for the marginal plots
        # ax_xhist.set_ylim((min_value, max_value))
        # ax_yhist.set_ylim(min_phred, max_phred)

        ##########
        # Plotting
        ##########

        # histplot x-axis - Size distribution
        for hue_grp in hue_groups:
            sns.distplot(subdata[hue_grp][x], color=colors[hue_grp],
                         ax=ax_xhist, bins=len_logbins)

        # histplot y-axis - Phred score
        for hue_grp in hue_groups:
            sns.distplot(subdata[hue_grp][y], color=colors[hue_grp],
                         ax=ax_yhist, vertical=True, bins=phred_bins)

        # main scatterplot
        # note: must be after the histplots else ax_yhist messes up
        for hue_grp in hue_groups:
            sns.regplot(data=subdata[hue_grp], fit_reg=False,
                        x=x, y=y, ax=ax_main, color=colors[hue_grp],
                        scatter_kws=scatter_kws)

        # despine
        for myax in [ax_yhist, ax_xhist]:
            sns.despine(ax=myax, bottom=False, top=True, left=False, right=True, trim=False)
            plt.setp(myax.get_xticklabels(), visible=False)
            plt.setp(myax.get_yticklabels(), visible=False)

        # topright
        ax_legend = plt.subplot(grid[0, 1])  # , sharey=ax_main)
        ax_legend.set_facecolor('white')
        plt.setp(ax_legend.get_xticklabels(), visible=False)
        plt.setp(ax_legend.get_yticklabels(), visible=False)

        # Hide label and grid for histogram plots
        ax_xhist.set_xlabel('')
        ax_yhist.set_ylabel('')
        ax_legend.grid(False)  # hide grid
        ax_xhist.grid(False)
        ax_yhist.grid(False)

        ax_legend.legend(handles=legend_mapping)
        plt.close()

        return dict(fig=fig, gridspec=grid)

    @staticmethod
    def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs):
        """
        # https://stackoverflow.com/questions/41577705/how-does-2d-kernel-density-estimation-in-python-sklearn-work
        Build 2D kernel density estimate (KDE).
        """

        from sklearn.neighbors import KernelDensity

        # create grid of sample locations (default: 100x100)
        xx, yy = np.mgrid[x.min():x.max():xbins,
                 y.min():y.max():ybins]

        xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
        xy_train = np.vstack([y, x]).T

        kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
        kde_skl.fit(xy_train)

        # score_samples() returns the log-likelihood of the samples
        z = np.exp(kde_skl.score_samples(xy_sample))
        return xx, yy, np.reshape(z, xx.shape)

    @staticmethod
    def plot_total_reads_vs_time(d, out):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :param out: path to output png file
        :return:
        
        TODO -> use numpy to handle the plot data, on row per sample?
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()
        t_pass = list()  # time
        t_fail = list()

        for seq_id, seq in d.items():
            t = seq.time_string
            if seq.flag == 'pass':
                t_pass.append(t)
            else:
                t_fail.append(t)

        # Find the smallest datetime value
        t_zero_pass = None
        t_zero_fail = None
        if t_pass:
            t_zero_pass = min(t_pass)

        if t_fail:
            t_zero_fail = min(t_fail)

        if t_pass and t_fail:
            t_zero = min(t_zero_pass, t_zero_fail)
        elif t_pass:
            t_zero = t_zero_pass
        elif t_fail:
            t_zero = t_zero_fail
        else:
            raise Exception('No data!')

        # Prepare datetime value for plotting
        # Convert time object in hours from beginning of run
        y_pass = None
        y_fail = None
        if t_pass:
            t_pass[:] = [x - t_zero for x in t_pass]  # Subtract t_zero for the all time points
            t_pass.sort()  # Sort
            t_pass[:] = [x.days * 24 + x.seconds / 3600 for x in t_pass]  # Convert to hours (float)
            y_pass = range(1, len(t_pass) + 1, 1)  # Create range. 1 time point equals 1 read

        if t_fail:
            t_fail[:] = [x - t_zero for x in t_fail]
            t_fail.sort()
            t_fail[:] = [x.days * 24 + x.seconds / 3600 for x in t_fail]
            y_fail = range(1, len(t_fail) + 1, 1)

        # Create plot
        if t_pass and t_fail:
            ax.plot(t_pass, y_pass, color='blue')
            ax.plot(t_fail, y_fail, color='red')
            ax.legend(['Pass', 'Fail'])
        elif t_pass:
            ax.plot(t_pass, y_pass, color='blue')
            ax.legend(['Pass'])
        elif t_fail:
            ax.plot(t_fail, y_fail, color='red')
            ax.legend(['Fail'])

        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Total read yield')
        plt.tight_layout()
        fig.savefig(out + "/total_reads_vs_time.png")

    @staticmethod
    def plot_reads_per_sample_vs_time(d, out):
        """
        Plot yield per sample. Just the pass reads
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        # fig, ax = plt.subplots()
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                my_sample_dict[seq.name][seq_id] = seq.time_string

        # Order the dictionary by keys
        od = OrderedDict(sorted(my_sample_dict.items()))

        # Make the plot
        legend_names = list()
        for name, seq_ids in od.items():
            legend_names.append(name)
            ts_pass = list()
            for seq, time_tag in seq_ids.items():
                ts_pass.append(time_tag)

            ts_zero = min(ts_pass)
            ts_pass[:] = [x - ts_zero for x in ts_pass]  # Subtract t_zero for the all time points
            ts_pass.sort()  # Sort
            ts_pass[:] = [x.days * 24 + x.seconds / 3600 for x in ts_pass]  # Convert to hours (float)
            ys_pass = range(1, len(ts_pass) + 1, 1)  # Create range. 1 time point equals 1 read

            # ax.plot(ts_pass, ys_pass)
            ax.plot(ts_pass, ys_pass,
                    label="%s (%s)" % (name, "{:,}".format(max(ys_pass))))
            # ax.legend(legend_names)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Pass reads per sample')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        # comma-separated numbers to the y axis
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()  #
        fig.savefig(out + "/reads_per_sample_vs_time.png")

    @staticmethod
    def plot_reads_per_sample_pie(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        # Fetch required information
        my_sample_dict = defaultdict(list)
        for seq_id, seq in d.items():
            my_sample_dict[seq_id] = [seq.name, seq.flag]

        df = pd.DataFrame.from_dict(my_sample_dict, orient='index', columns=['name', 'flag'])
        df_all = df.groupby(['name']).count()
        df_pass = df[df['flag'] == 'pass'].groupby(['name']).count()
        df_fail = df[df['flag'] == 'fail'].groupby(['name']).count()

        if not df_fail.empty:
            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

            # Make the plots
            titles = ['All', 'Pass', 'Fail']
            for i, my_df in enumerate([df_all, df_pass, df_fail]):
                my_df.columns = ['count']
                my_df = my_df.sort_values(['count'], ascending=False)  # sort dataframe for better looking pie chart
                data = list(my_df['count'])
                data_sum = sum(data)
                labels = list(my_df.index)
                for j, l in enumerate(labels):
                    labels[j] = "{:s} ({:,} reads, {:.1f}%)".format(l, data[j], round(data[j] / data_sum * 100, 1))

                axs[i].pie(data, labels=labels, wedgeprops={'linewidth': 2, 'edgecolor': 'w'})
                axs[i].set_title(titles[i])
        else:
            fig, ax = plt.subplots(figsize=(10, 4))

            df_pass.columns = ['count']
            df_pass = df_pass.sort_values(['count'], ascending=False)  # sort dataframe for better looking pie chart
            data = list(df_pass['count'])
            data_sum = sum(data)
            labels = list(df_pass.index)
            for j, l in enumerate(labels):
                labels[j] = "{:s} ({:,} reads, {:.1f}%)".format(l, data[j], round(data[j] / data_sum * 100, 1))

            ax.pie(data, labels=labels, wedgeprops={'linewidth': 2, 'edgecolor': 'w'})
            ax.set_title('Pass')

        # Add label to axes
        plt.subplots_adjust(hspace=1)
        fig.suptitle('Distribution of reads among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/reads_per_sample_pie.png")

    @staticmethod
    def plot_bp_per_sample_pie(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        # Fetch required information
        my_sample_dict = defaultdict(list)
        for seq_id, seq in d.items():
            my_sample_dict[seq_id] = [seq.name, seq.flag, seq.length]

        df = pd.DataFrame.from_dict(my_sample_dict, orient='index', columns=['name', 'flag', 'length'])
        df_all = df.groupby(['name']).sum()
        df_pass = df[df['flag'] == 'pass'].groupby(['name']).sum()
        df_fail = df[df['flag'] == 'fail'].groupby(['name']).sum()

        if not df_fail.empty:
            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

            # Make the plots
            titles = ['All', 'Pass', 'Fail']
            for i, my_df in enumerate([df_all, df_pass, df_fail]):
                my_df.columns = ['count']
                my_df = my_df.sort_values(['count'], ascending=False)  # sort dataframe for better looking pie chart
                data = list(my_df['count'])
                data_sum = sum(data)
                labels = list(my_df.index)
                for j, l in enumerate(labels):
                    labels[j] = "{:s} ({:,} bp, {:.1f}%)".format(l, data[j], round(data[j] / data_sum * 100, 1))

                axs[i].pie(data, labels=labels, wedgeprops={'linewidth': 2, 'edgecolor': 'w'})
                axs[i].set_title(titles[i])
        else:
            fig, ax = plt.subplots(figsize=(10, 4))

            df_pass.columns = ['count']
            df_pass = df_pass.sort_values(['count'], ascending=False)  # sort dataframe for better looking pie chart
            data = list(df_pass['count'])
            data_sum = sum(data)
            labels = list(df_pass.index)
            for j, l in enumerate(labels):
                labels[j] = "{:s} ({:,} bp, {:.1f}%)".format(l, data[j], round(data[j] / data_sum * 100, 1))

            ax.pie(data, labels=labels, wedgeprops={'linewidth': 2, 'edgecolor': 'w'})
            ax.set_title('Pass')

        # Add label to axes
        plt.subplots_adjust(hspace=1)
        fig.suptitle('Distribution of bp among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/bp_per_sample_pie.png")

    @staticmethod
    def plot_bp_per_sample_vs_time(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                my_sample_dict[seq.name][seq_id] = (seq.time_string, seq.length)  #tuple

        # Order the dictionary by keys
        od = OrderedDict(sorted(my_sample_dict.items()))

        # Make the plot
        for name, seq_ids in od.items():
            ts_pass = list()
            for seq, data_tuple in seq_ids.items():
                ts_pass.append(data_tuple)

            # Prepare x values (time)
            ts_zero = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list
            ts_pass1 = [tuple(((x - ts_zero), y)) for x, y in ts_pass]  # Subtract t_zero for the all time points
            ts_pass1.sort(key=lambda x: x[0])  # Sort according to first element in tuples
            ts_pass2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x,y in ts_pass1]  # Convert to hours (float)
            x_values = [x for x, y in ts_pass2]  # Only get the fist value of the ordered tuples

            # Prepare y values in a cumulative way
            c = 0
            y_values = list()
            for x, y in ts_pass2:
                y = y + c
                y_values.append(y)
                c = y

            # Plot values per sample
            ax.plot(x_values, y_values,
                    label="%s (%s)" % (name, "{:,}".format(max(y_values))))

        # ax.legend(legend_names, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # New
        # Add axes labels
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Yield per sample in base pair\n("pass" only)')
        # ax.ticklabel_format(useOffset=False)  # Disable the offset on the x-axis
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        # Save figure to file
        fig.savefig(out + "/bp_per_sample_vs_time.png")

    @staticmethod
    def plot_total_bp_vs_time(d, out):
        """
        Sequence length vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        # Fetch required information
        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                ts_pass.append(tuple((seq.time_string, seq.length)))
            else:
                ts_fail.append(tuple((seq.time_string, seq.length)))

        ts_zero_pass = list()
        ts_zero_fail = list()
        if ts_pass:
            ts_zero_pass = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_fail:
            ts_zero_fail = min(ts_fail, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_pass and ts_fail:
            ts_zero = min(ts_zero_pass, ts_zero_fail)
        elif ts_pass:
            ts_zero = ts_zero_pass
        else:  # elif ts_fail:
            ts_zero = ts_zero_fail

        x_pass_values = None
        y_pass_values = None
        x_fail_values = None
        y_fail_values = None
        if ts_pass:
            ts_pass1 = [tuple(((x - ts_zero), y)) for x, y in ts_pass]  # Subtract t_zero for the all time points
            ts_pass1.sort(key=lambda x: x[0])  # Sort according to first element in tuple
            ts_pass2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_pass1]  # Convert to hours (float)
            x_pass_values = [x for x, y in ts_pass2]
            c = 0
            y_pass_values = list()
            for x, y in ts_pass2:
                y = y + c
                y_pass_values.append(y)
                c = y
        if ts_fail:
            ts_fail1 = [tuple(((x - ts_zero), y)) for x, y in ts_fail]
            ts_fail1.sort(key=lambda x: x[0])
            ts_fail2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_fail1]
            x_fail_values = [x for x, y in ts_fail2]
            c = 0
            y_fail_values = list()
            for x, y in ts_fail2:
                y = y + c
                y_fail_values.append(y)
                c = y

        # Print plot
        if ts_pass and ts_fail:
            ax.plot(x_pass_values, y_pass_values, color='blue')
            ax.plot(x_fail_values, y_fail_values, color='red')
            ax.legend(['Pass', 'Fail'])
        elif ts_pass:
            ax.plot(x_pass_values, y_pass_values, color='blue')
            ax.legend(['Pass'])
        else:  # elif ts_fail:
            ax.plot(x_fail_values, y_fail_values, color='red')
            ax.legend(['Fail'])
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Total yield in base pair')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/total_bp_vs_time.png")

    @staticmethod
    def plot_quality_vs_time(d, out):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots(figsize=(10, 6))

        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                #         average_phred = round(average_phred_full, 1)
                ts_pass.append(tuple((seq.time_string, round(seq.average_phred, 1))))
            else:
                ts_fail.append(tuple((seq.time_string, round(seq.average_phred, 1))))

        ts_zero_pass = None
        ts_zero_fail = None
        if ts_pass:
            ts_zero_pass = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_fail:
            ts_zero_fail = min(ts_fail, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_pass and ts_fail:
            ts_zero = min(ts_zero_pass, ts_zero_fail)
        elif ts_pass:
            ts_zero = ts_zero_pass
        else:  # elif ts_fail:
            ts_zero = ts_zero_fail

        ts_pass3 = list()
        ts_fail3 = list()
        if ts_pass:
            ts_pass1 = [tuple(((x - ts_zero), y)) for x, y in ts_pass]  # Subtract t_zero for the all time points
            ts_pass1.sort(key=lambda x: x[0])  # Sort according to first element in tuple
            ts_pass2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_pass1]  # Convert to hours (float)
            ts_pass3 = [tuple((int(np.round(x)), y)) for x, y in ts_pass2]  # Round hours

            df_pass = pd.DataFrame(list(ts_pass3), columns=['Sequencing time interval (h)', 'Phred score'])  # Convert to dataframe
            df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        if ts_fail:
            ts_fail1 = [tuple(((x - ts_zero), y)) for x, y in ts_fail]
            ts_fail1.sort(key=lambda x: x[0])
            ts_fail2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_fail1]
            ts_fail3 = [tuple((int(np.round(x)), y)) for x, y in ts_fail2]

            df_fail = pd.DataFrame(list(ts_fail3), columns=['Sequencing time interval (h)', 'Phred score'])
            df_fail['Flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

        # Account if there is no fail data or no pass data
        if ts_fail3 and ts_pass3:
            frames = [df_pass, df_fail]
            data = pd.concat(frames)  # Merge dataframes
        elif ts_pass3:
            data = df_pass
        else:  # elif ts_fail3:
            data = df_fail

        # Account if there is no fail data or no pass data
        if ts_fail3 and ts_pass3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, hue='Flag', split=True, inner=None)
            g.figure.suptitle('Sequence quality over time')
        elif ts_pass3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, inner=None)
            g.figure.suptitle('Sequence quality over time (pass only)')
        else:  # elif ts_fail3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, inner=None)
            g.figure.suptitle('Sequence quality over time (fail only)')

        # Major ticks every 4 hours
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        # https://matplotlib.org/2.0.2/examples/ticks_and_spines/tick-locators.html
        def my_formater(val, pos):
            val_str = '{}-{}'.format(int(val), int(val + 1))
            return val_str

        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4))

        if ts_fail:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/quality_vs_time.png")

    @staticmethod
    def plot_phred_score_distribution(d, out):
        """
        Frequency of phred scores
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [seq.average_phred, seq.flag]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'pass']
        df_fail = df.loc[df['flag'] == 'fail']

        average_qual_pass = round(functions.compute_average_quality(df_pass['Phred score'].tolist(),
                                                                    df_pass.shape[0]), 1)

        ax.hist(df_pass['Phred score'], histtype='stepfilled', color='blue', alpha=0.6,
                label='Pass (Avg: {})'.format(average_qual_pass))

        if not df_fail.empty:
            average_qual_fail = round(functions.compute_average_quality(df_fail['Phred score'].tolist(),
                                                                        df_fail.shape[0]), 1)
            ax.hist(df_fail['Phred score'], histtype='stepfilled', color='red', alpha=0.6,
                    label='Fail (Avg: {})'.format(average_qual_fail))

        plt.legend()

        ax.set(xlabel='Phred score', ylabel='Frequency', title='Phred score distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/phred_score_distribution.png")

    @staticmethod
    def plot_length_distribution(d, out):
        """
        Frequency of sizes. Bins auto-sized based on length distribution. Log scale x-axis.
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [seq.length, seq.flag]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length', 'flag'])
        df_pass = df.loc[df['flag'] == 'pass']
        df_fail = df.loc[df['flag'] == 'fail']

        # Set bin sized for histogram
        min_len = min(df['Length'])
        max_len = max(df['Length'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        len_logbins = np.logspace(min_exp, max_exp, 50)

        average_len_pass = int(df_pass['Length'].mean())
        ax.hist(df_pass['Length'], histtype='stepfilled', color='blue', alpha=0.6,
                label="Pass (Avg: {} bp)".format(average_len_pass), bins=len_logbins)

        if not df_fail.empty:
            average_len_fail = int(df_fail['Length'].mean())
            ax.hist(df_fail['Length'], histtype='stepfilled', color='red', alpha=0.6,
                    label="Fail (Avg: {} bp)".format(average_len_fail), bins=len_logbins)

        plt.legend()
        plt.xscale('log')
        ax.set(xlabel='Read length (bp)', ylabel='Frequency', title='Read length distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/length_distribution.png")

    @staticmethod
    def test_plot(d, out):
        """ TAKES 50 MIN to run! The KDE part takes way too long."""
        # from scipy import stat

        qs_pass = list()
        qs_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                qs_pass.append(tuple((seq.length, seq.average_phred)))
            else:
                qs_fail.append(tuple((seq.length, seq.average_phred)))

        df_pass = pd.DataFrame(list(qs_pass), columns=['Length (bp)', 'Phred score'])
        df_pass['flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        df_fail = pd.DataFrame(list(qs_fail), columns=['Length (bp)', 'Phred score'])
        df_fail['flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

        df_concatenated = pd.concat([df_pass, df_fail])

        # Find min and max length values
        min_len = min(df_concatenated.ix[:, 0])
        max_len = max(df_concatenated.ix[:, 0])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))
        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 50)

        # Set y-axis limits
        min_phred = min(df_concatenated.ix[:, 1])
        max_phred = max(df_concatenated.ix[:, 1])
        # phred_range = max_phred - min_phred
        phred_bins = np.linspace(min_phred, max_phred, 50)

        # Set plot
        g = sns.JointGrid(x='Length (bp)', y='Phred score', data=df_pass, space=0,
                          xlim=(min_value, max_value), ylim=(min_phred, max_phred))

        # x = df_pass['Length (bp)']
        # y = df_pass['Phred score']
        # kde = stats.gaussian_kde([x, y])
        # xx, yy = np.mgrid[min(x):max(x):(max(x)-min(x))/100, min(y):max(y):(max(y)-min(y))/100]
        # density = kde(np.c_[xx.flat, yy.flat].T).reshape(xx.shape)

        # Draw plot in top marginal area (x) for pass
        # sns.kdeplot(df_pass['Length (bp)'], ax=g.ax_marg_x, legend=False, shade=True)
        g.ax_marg_x.hist(df_pass['Length (bp)'], color='blue', alpha=0.6, bins=len_logbins)
        # Draw plot in right marginal area (y) for pass
        # sns.kdeplot(df_pass['Phred score'], ax=g.ax_marg_y, legend=False, shade=True, vertical=True)
        g.ax_marg_y.hist(df_pass['Phred score'], color='blue', alpha=0.6, bins=phred_bins, orientation="horizontal")
        # Draw plot in joint area for pass
        sns.kdeplot(df_pass['Length (bp)'], df_pass['Phred score'], ax=g.ax_joint,
                    n_levels=25, cmap="Blues", shade=True, shade_lowest=False, legend=True)
        # g.ax_joint(sns.kdeplot, n_levels=30, cmap="Blues", shade=True, shade_lowest=False)
        # g.ax_joint.contourf(xx, yy, density, 10, cmap="Blues")

        # g.x = df_fail['Length (bp)']
        # g.y = df_fail['Length (bp)']
        #
        # Draw plot in top marginal area (x) for fail
        # sns.kdeplot(df_pass['Length (bp)'], ax=g.ax_marg_x, legend=False, shade=True)
        g.ax_marg_x.hist(df_fail['Length (bp)'], color='red', alpha=0.6, bins=len_logbins)
        # Draw plot in right marginal area (y) for fail
        # sns.kdeplot(df_pass['Phred score'], ax=g.ax_marg_y, legend=False, shade=True, vertical=True)
        g.ax_marg_y.hist(df_fail['Phred score'], color='red', alpha=0.6, bins=phred_bins, orientation="horizontal")
        # Draw plot in joint area for fail
        sns.kdeplot(df_fail['Length (bp)'], df_fail['Phred score'], ax=g.ax_joint,
                    n_levels=25, cmap="Reds", shade=True, shade_lowest=False)

        # Set xscale to log
        g.ax_joint.set_xscale('log')
        g.ax_marg_x.set_xscale('log')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Add legend
        g.ax_marg_x.legend(('pass', 'fail'), loc='upper right')

        # Save figure to file
        g.savefig(out + "/quality_vs_length_kde.png")

    @staticmethod
    def plot_quality_vs_length_kde(d, out):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        sns.set(style="ticks")

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [seq.length, seq.average_phred, seq.flag]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length (bp)', 'Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'pass']
        df_fail = df.loc[df['flag'] == 'fail']

        # Set x-axis limits
        min_len = min(df.ix[:, 0])
        max_len = max(df.ix[:, 0])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))

        # Set bin sized for histogram
        # len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set y-axis limits
        min_phred = min(df.ix[:, 1])
        max_phred = max(df.ix[:, 1])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Do the Kernel Density Estimation (KDE)
        x = df_pass['Length (bp)']
        y = df_pass['Phred score']
        xx, yy, density = FastqPlots.kde2D(x, y, 1)

        # Create grid object
        g = sns.JointGrid(x='Length (bp)', y='Phred score', data=df_pass, space=0)

        cset1 = g.ax_joint.contourf(xx, yy, density, levels=25, cmap="Blues", alpha=0.6)
        # don't shade last contour
        # https://github.com/mwaskom/seaborn/blob/master/seaborn/distributions.py
        cset1.collections[0].set_alpha(0)

        g.ax_marg_x.hist(x, histtype='stepfilled', color='blue', alpha=0.6, bins=len_logbins)
        g.ax_marg_y.hist(y, histtype='stepfilled', color='blue', alpha=0.6, bins=phred_bins, orientation="horizontal")

        # Set main plot x axis scale to log
        g.ax_joint.set_xscale('log')
        g.ax_marg_x.set_xscale('log')
        g.ax_joint.set_xlim((min_value, max_value))

        ####
        # Do the same for the fail reads
        ####

        if not df_fail.empty:
            g.x = df_fail['Length (bp)']
            g.y = df_fail['Phred score']

            # Do the Kernel Density Estimation (KDE)
            x = df_fail['Length (bp)']
            y = df_fail['Phred score']
            xx, yy, density = FastqPlots.kde2D(x, y, 1)

            cset2 = g.ax_joint.contourf(xx, yy, density, levels=25, cmap="Reds", alpha=0.6)

            # don't shade last contour
            # https://github.com/mwaskom/seaborn/blob/master/seaborn/distributions.py
            cset2.collections[0].set_alpha(0)

            g.ax_marg_x.hist(x, histtype='stepfilled', color='red', alpha=0.6, bins=len_logbins)
            g.ax_marg_y.hist(y, histtype='stepfilled', color='red', alpha=0.6, bins=phred_bins,
                             orientation="horizontal")

        # Add legend to the joint plot area
        # https://matplotlib.org/tutorials/intermediate/legend_guide.html
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='Pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='Fail')
        if not df_fail.empty:
            g.ax_joint.legend(handles=[blue_patch, red_patch], loc='best')
        else:
            g.ax_joint.legend(handles=[blue_patch], loc='best')

        # Add legend to top margin_plot area
        # g.ax_marg_x.legend(('pass', 'fail'), loc='upper right')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Save figure to file
        g.savefig(out + "/quality_vs_length_kde.png")

    @staticmethod
    def plot_quality_vs_length_hex(d, out):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        sns.set(style="ticks")

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [seq.length, seq.average_phred, seq.flag]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length (bp)', 'Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'pass']
        df_fail = df.loc[df['flag'] == 'fail']

        # Set x-axis limits
        min_len = min(df['Length (bp)'])
        max_len = max(df['Length (bp)'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))

        # Set bin sized for histogram
        # len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set y-axis limits
        min_phred = min(df['Phred score'])
        max_phred = max(df['Phred score'])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Do the Kernel Density Estimation (KDE)
        x = df_pass['Length (bp)']
        y = df_pass['Phred score']

        # Create grid object
        g = sns.JointGrid(x='Length (bp)', y='Phred score', data=df_pass, space=0)

        g.ax_joint.hexbin(x, y, gridsize=50, cmap="Blues", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
        g.ax_joint.axis([min_value, max_value, min_phred, max_phred])
        g.ax_marg_x.hist(x, histtype='stepfilled', color='blue', alpha=0.6, bins=len_logbins)
        g.ax_marg_y.hist(y, histtype='stepfilled', color='blue', alpha=0.6, bins=phred_bins, orientation="horizontal")

        # Set main plot x axis scale to log
        g.ax_joint.set_xscale('log')
        g.ax_marg_x.set_xscale('log')
        g.ax_joint.set_xlim((min_value, max_value))

        ####
        # Do the same for the fail reads
        ####

        if not df_fail.empty:
            g.x = df_fail['Length (bp)']
            g.y = df_fail['Phred score']

            # Do the Kernel Density Estimation (KDE)
            x = df_fail['Length (bp)']
            y = df_fail['Phred score']

            g.ax_joint.hexbin(x, y, gridsize=50, cmap="Reds", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
            g.ax_marg_x.hist(x, histtype='stepfilled', color='red', alpha=0.6, bins=len_logbins)
            g.ax_marg_y.hist(y, histtype='stepfilled', color='red', alpha=0.6, bins=phred_bins,
                             orientation="horizontal")

        # Add legend to the joint plot area
        # https://matplotlib.org/tutorials/intermediate/legend_guide.html
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='Pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='Fail')
        if not df_fail.empty:
            g.ax_joint.legend(handles=[blue_patch, red_patch], loc='upper left')
        else:
            g.ax_joint.legend(handles=[blue_patch], loc='upper left')

        # Add legend to top margin_plot area
        # g.ax_marg_x.legend(('pass', 'fail'), loc='upper right')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Save figure to file
        g.savefig(out + "/quality_vs_length_hex.png")

    @staticmethod
    def plot_quality_vs_length_scatter(d, out):
        qs_pass = list()
        qs_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                qs_pass.append(tuple((seq.length, seq.average_phred)))
            else:
                qs_fail.append(tuple((seq.length, seq.average_phred)))

        df_pass = pd.DataFrame(list(qs_pass), columns=['Length (bp)', 'Phred Score'])
        df_pass['flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        df_fail = pd.DataFrame(list(qs_fail), columns=['Length (bp)', 'Phred Score'])
        df_fail['flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

        df_concatenated = pd.concat([df_pass, df_fail])

        fig = plt.figure(figsize=(10, 6))

        if qs_fail:
            FastqPlots.jointplot_w_hue(data=df_concatenated, x='Length (bp)', y='Phred Score',
                                 hue='flag', figsize=(10, 6), fig=fig, colormap=['blue', 'red'],
                                 scatter_kws={'s': 1, 'alpha': 0.1})
        else:
            FastqPlots.jointplot_w_hue(data=df_pass, x='Length (bp)', y='Phred Score',
                                 hue='flag', figsize=(10, 6), fig=fig, colormap=['blue'],
                                 scatter_kws={'s': 1, 'alpha': 0.1})

        fig.savefig(out + "/quality_vs_length_scatter.png")

    @staticmethod
    def plot_test_old(d, out):
        """
        seaborn jointplot (length vs quality). More manual.
        :param d: Dictionary
        :param out: path to output png file
        :return:
        """

        from matplotlib import gridspec
        from scipy.stats import gaussian_kde

        qs_pass = list()
        qs_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                qs_pass.append(tuple((seq.length, seq.average_phred)))
            else:
                qs_fail.append(tuple((seq.length, seq.average_phred)))

        df_pass = pd.DataFrame(list(qs_pass), columns=['Length (bp)', 'Phred Score'])
        df_pass['flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        df_fail = pd.DataFrame(list(qs_fail), columns=['Length (bp)', 'Phred Score'])
        df_fail['flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

        df_concatenated = pd.concat([df_pass, df_fail])

        min_len = pd.DataFrame.min(df_concatenated['Length (bp)'])
        min_exp = np.log10(min_len)
        min_value = float(10 ** (min_exp - 0.1))
        if min_value <= 0:
            min_value = 1
        max_len = pd.DataFrame.max(df_concatenated['Length (bp)'])
        max_exp = np.log10(max_len)
        max_value = float(10 ** (max_exp + 0.1))

        min_phred = df_concatenated['Phred Score'].min()
        max_phred = df_concatenated['Phred Score'].max()

        binwidth = int(np.round(10 * np.log10(max_len - min_len)))
        logbins = np.logspace(np.log10(min_len), np.log10(max_len), binwidth)

        # Initialize the figure
        grid = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 4],
                               hspace=0, wspace=0)

        ax_main = plt.subplot(grid[1, 0])
        ax_xhist = plt.subplot(grid[0, 0], sharex=ax_main)
        ax_yhist = plt.subplot(grid[1, 1])  # , sharey=ax_main)

        gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 4],
                               hspace=0, wspace=0)

        # Create scatter plot
        fig = plt.figure(figsize=(10, 6))  # In inches
        ax = plt.subplot(gs[1, 0])
        pass_ax = ax.scatter(df_pass['Length (bp)'], df_pass['Phred Score'],
                         color='blue', alpha=0.1, s=1, label='pass')

        fail_ax = ax.scatter(df_fail['Length (bp)'], df_fail['Phred Score'],
                         color='red', alpha=0.1, s=1, label='fail')
        pass_ax.xlim((0, max_value))
        pass_ax.ylim((min_phred, max_phred))
        fail_ax.xlim((0, max_value))
        fail_ax.ylim((min_phred, max_phred))

        # Create Y-marginal (right) -> Phred score
        axr = plt.subplot(gs[1, 1], sharey=fail_ax, frameon=False,
                          xticks=[])  # xlim=(0, 1), ylim = (ymin, ymax) xticks=[], yticks=[]
        axr.hist(df_pass['Phred Score'], color='#5673E0', orientation='horizontal', density=True)

        # Create X-marginal (top) -> length
        axt = plt.subplot(gs[0, 0], sharex=pass_ax, frameon=False,
                          yticks=[])  # xticks = [], , ) #xlim = (xmin, xmax), ylim=(0, 1)
        axt.set_xscale('log')
        axt.set_yscale('log')
        axt.hist(df_pass['Length (bp)'], color='#5673E0', density=True)

        legend_ax = plt.subplot(gs[0, 1], frameon=False)  # top right
        legend_ax.legend = ((pass_ax, fail_ax), ('pass', 'fail'))

        # Bring the marginals closer to the scatter plot
        fig.tight_layout()

        # sns.set_style("ticks")
        # g = sns.JointGrid(x='Length (bp)', y='Phred Score', data=df_concatenated,
        #                   xlim=[min_value, max_value], ylim=[min_phred, max_phred],
        #                   space=0)
        #
        # ax = g.ax_joint
        # ax.cla()  # clear the 2-D plot
        # ax.set_xscale('log')
        # g.ax_marg_x.set_xscale('log')
        #
        # # g.plot_joint(sns.kdeplot, shade=True, n_levels=100)
        # # g.plot_joint(sns.regplot, scatter_kws={"color":"darkred","alpha":0.1,"s":1}, fit_reg=False)
        # # g.plot_joint(sns.lmplot, x='Length (bp)', y='Phred Score', data=df_concatenated,  hue='flag')
        # # g.plot_joint(sns.lmplot, x='Length (bp)', y='Phred Score', data=df_concatenated, hue='flag', fit_reg=False)
        # g.plot_joint(plt.scatter,)
        # g.plot_marginals(sns.distplot, kde=False)

        # g.fig.set_figwidth(8)
        # g.fig.set_figheight(4)

        # Save
        # g.savefig(out + "/test.png")
        fig.savefig(out + "/test.png")

    @staticmethod
    def plot_reads_vs_bp_per_sample(d, out):
        # Fetch required information
        my_sample_dict = defaultdict()  # to get the lengths (bp)
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                if not my_sample_dict[seq.name]:
                    my_sample_dict[seq.name] = [seq.length]
                else:
                    my_sample_dict[seq.name].append(seq.length)

        # Order the dictionary by keys
        od = OrderedDict(sorted(my_sample_dict.items()))

        # Create pandas dataframe
        df = pd.DataFrame(columns=['Sample', 'bp', 'reads'])
        for name, size_list in od.items():
            df = df.append({'Sample': name, 'bp': sum(size_list), 'reads': len(size_list)}, ignore_index=True)

        fig, ax1 = plt.subplots(figsize=(10, 6))  # In inches

        ind = np.arange(len(df['Sample']))
        width = 0.35
        p1 = ax1.bar(ind, df['bp'], width, color='#4B9BFF', bottom=0, edgecolor='black')

        ax2 = ax1.twinx()
        p2 = ax2.bar(ind+width, df['reads'], width, color='#FFB46E', bottom=0, edgecolor='black')

        ax1.set_title('Total Size Versus Total Reads Per Sample')
        ax1.set_xticks(ind + width / 2)
        ax1.set_xticklabels(tuple(df['Sample']), rotation=45, ha='right')

        ax1.grid(False)
        ax2.grid(False)

        ax1.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        ax2.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

        ax1.legend((p1[0], p2[0]), ('bp', 'reads'), bbox_to_anchor=(1.1, 1), loc=2)
        ax1.yaxis.set_units('Total bp')
        ax2.yaxis.set_units('Total reads')
        ax1.autoscale_view()

        plt.tight_layout()
        fig.savefig(out + "/reads_vs_bp_per_sample.png")

        #########################################
        # ax = df.plot(kind='bar', secondary_y='reads', title='bp versus reads',
        #              x='Sample', y=['bp', 'reads'], mark_right=False)
        #
        # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        # ax.set_ylabel("Total bp")
        # plt.grid(False)
        # ax.legend(bbox_to_anchor=(1.3, 1), loc=2)
        # plt.legend(bbox_to_anchor=(1.3, 0.92), loc=2)
        #
        # plt.tight_layout()
        #
        # fig = ax.get_figure()
        # fig.savefig(out + "/reads_vs_bp_per_sample.png")

        ##########################################
        # df = pd.DataFrame(columns=['Sample', 'Value', 'Info'])
        # for name, size_list in my_sample_dict.items():
        #     df = df.append({'Sample': name, 'Value': sum(size_list), 'Info': 'Total bp'}, ignore_index=True)
        #     df = df.append({'Sample': name, 'Value': len(size_list), 'Info': 'Reads'}, ignore_index=True)
        #
        # g = sns.catplot(x='Sample', y='Value', hue='Info', data=df, kind='bar')
        #
        # plt.tight_layout()
        # g.savefig(out + "/reads_vs_bp_per_sample.png")

    @staticmethod
    def plot_pores_output_vs_time_total(d, out):

        time_list = list()
        for seq_id, seq in d.items():
            time_list.append(seq.time_string)

        time_list = sorted(time_list)  # order list
        time_zero = min(time_list)  # find smallest datetime value
        time_list1 = [x - time_zero for x in time_list]  # Subtract t_zero for the all time points
        time_list2 = [x.days * 1440 + x.seconds / 60 for x in time_list1]  # Convert to minutes (float)
        time_list3 = [int(np.round(x)) for x in time_list2]  # Round minutes
        # Check how many 15-minute bins are required to plot all the data
        nbins = max(time_list3) / 15 if max(time_list3) % 15 == 0 else int(max(time_list3) / 15) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(time_list3), max(time_list3), nbins)  # every 15 min

        # Generate counts for each bin
        hist, edges = np.histogram(time_list3, bins=x_bins, density=False)

        fig, ax = plt.subplots()

        # Plot the data
        g = sns.scatterplot(data=hist, x_bins=edges, legend=False, size=3, alpha=0.5, linewidth=0)

        # Adjust format of numbers for y axis: "1000000" -> "1,000,000"
        g.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

        # Change x axis labels chunk-of-15-min to hours
        def numfmt(m, pos):  # your custom formatter function: divide by 100.0
            h = '{}'.format(m / 4)
            return h
        ax.xaxis.set_major_formatter(FuncFormatter(numfmt))

        # Major ticks every 4 hours
        def my_formater(val, pos):
            val_str = '{}'.format(int(val/4))
            return val_str

        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        # ticks = range(0, ceil(max(time_list3) / 60) + 4, 4)
        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        # ax.xaxis.set_major_locator(MaxNLocator(len(ticks), integer=True))
        ax.xaxis.set_major_locator(MultipleLocator(4 * 4))  # 4 block of 15 min per hour. Want every 4 hours

        # Add label to axes
        plt.title('Pores output over time')
        plt.ylabel('Reads per 15 minutes')
        plt.xlabel('Sequencing time (hours)')

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig = g.get_figure()  # Get figure from FacetGrid
        fig.savefig(out + "/pores_output_vs_time.png")

    @staticmethod
    def plot_pores_output_vs_time_all(d, out):

        import matplotlib.lines as mlines

        fig, ax = plt.subplots()

        time_list_all = list()
        time_list_pass = list()
        time_list_fail = list()

        for seq_id, seq in d.items():
            time_string = seq.time_string
            time_list_all.append(time_string)
            if seq.flag == 'pass':
                time_list_pass.append(time_string)
            else:
                time_list_fail.append(time_string)

        time_zero = min(time_list_all)  # find smallest datetime value

        # Plot all
        # time_list_all = sorted(time_list_all)  # order list
        time_list_all[:] = [x - time_zero for x in time_list_all]  # Subtract t_zero for the all time points
        time_list_all[:] = [x.days * 1440 + x.seconds / 60 for x in time_list_all]  # Convert to minutes (float)
        time_list_all[:] = [int(np.round(x)) for x in time_list_all]  # Round minutes
        # Check how many 15-minute bins are required to plot all the data
        nbins = max(time_list_all) / 15 if max(time_list_all) % 15 == 0 else int(max(time_list_all) / 15) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(time_list_all), max(time_list_all), nbins)  # every 15 min

        # If not fail, just draw the pass. Else, draw total, fail and pass
        if time_list_fail:  # fail reads might be missing if plotting filtered reads for example.
            # Plot pass - Assume always pass reads present
            time_list_pass[:] = [x - time_zero for x in time_list_pass]
            time_list_pass[:] = [x.days * 1440 + x.seconds / 60 for x in time_list_pass]
            time_list_pass[:] = [int(np.round(x)) for x in time_list_pass]
            # nbins = max(time_list_pass) / 15 if max(time_list_pass) % 15 == 0 else int(max(time_list_pass) / 15) + 1
            # x_bins = np.linspace(min(time_list_pass), max(time_list_pass), nbins)
            hist, edges = np.histogram(time_list_pass, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                            cmap='Blues')  # , marker='^'

            # Plot fail
            time_list_fail[:] = [x - time_zero for x in time_list_fail]
            time_list_fail[:] = [x.days * 1440 + x.seconds / 60 for x in time_list_fail]
            time_list_fail[:] = [int(np.round(x)) for x in time_list_fail]
            # nbins = max(time_list_fail) / 15 if max(time_list_fail) % 15 == 0 else int(max(time_list_fail) / 15) + 1
            # x_bins = np.linspace(min(time_list_fail), max(time_list_fail), nbins)
            hist, edges = np.histogram(time_list_fail, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                            cmap='Reds')  # , marker='x'

        # Generate counts for each bin
        hist, edges = np.histogram(time_list_all, bins=x_bins, density=False)
        # Plot the data
        sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                        cmap='Greens')  # , marker='o'

        # Adjust format of numbers for y axis: "1000000" -> "1,000,000"
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

        # Change x axis labels chunk-of-15-min to hours
        def numfmt(m, pos):
            h = '{}'.format(m / 4)
            return h

        ax.xaxis.set_major_formatter(FuncFormatter(numfmt))

        # Major ticks every 4 hours
        def my_formater(val, pos):
            val_str = '{}'.format(int(val / 4))
            return val_str

        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4 * 4))  # 4 block of 15 min per hour. Want every 4 hours

        # Add legend to the plot area
        # https://stackoverflow.com/questions/47391702/matplotlib-making-a-colored-markers-legend-from-scratch

        all_marker = mlines.Line2D([], [], color='green', alpha=0.6, label='All', marker='o',
                                   markersize=5, linestyle='None')
        pass_marker = mlines.Line2D([], [], color='blue', alpha=0.6, label='Pass', marker='o',
                                    markersize=5, linestyle='None')
        fail_marker = mlines.Line2D([], [], color='red', alpha=0.6, label='Fail', marker='o',
                                    markersize=5, linestyle='None')

        if time_list_fail:
            ax.legend(handles=[all_marker, pass_marker, fail_marker], loc='upper right')
        else:
            green_circle = mlines.Line2D([], [], color='green', alpha=0.6, label='Pass', marker='o',
                                         markersize=5, linestyle='None')
            ax.legend(handles=[all_marker], loc='upper right')

        # Add label to axes
        plt.title('Pores output over time')
        plt.ylabel('Reads per 15 minutes')
        plt.xlabel('Sequencing time (hours)')

        plt.tight_layout()  # Get rid of extra margins around the plot
        # fig = g.get_figure()  # Get figure from FacetGrid
        fig.savefig(out + "/pores_output_vs_time_all.png")

    @staticmethod
    def plot_channel_output_all(d, out):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :param out: path to output png file
        :return:
        """

        channel_dict = defaultdict()
        for seq_id, seq in d.items():
            channel_number = int(seq.channel)
            # Pass and fail apart
            if channel_number not in channel_dict:
                channel_dict[channel_number] = [0, 0]
            if seq.flag == 'pass':
                channel_dict[channel_number][0] += 1
            else:  # seq.flag == b'fail':
                channel_dict[channel_number][1] += 1

        # convert to Pandas dataframe
        df = pd.DataFrame.from_dict(channel_dict, orient='index', columns=['Pass', 'Fail'])
        df_all = pd.DataFrame()
        df_all['All'] = df['Pass'] + df['Fail']
        df_pass = df[['Pass']]  # The double square brackets keep the column name
        df_fail = df[['Fail']]

        # Plot
        if not df_fail['Fail'].sum() == 0:
            fig, axs = plt.subplots(nrows=3, figsize=(6, 12))

            for i, my_tuple in enumerate([(df_all, 'All', 'Greens'),
                                          (df_pass, 'Pass', 'Blues'),
                                          (df_fail, 'Fail', 'Reds')]):
                my_df = my_tuple[0]
                flag = my_tuple[1]
                cmap = my_tuple[2]

                maxval = max(my_df.index)  # maximum channel value
                layout = FastqPlots.make_layout(maxval=maxval)
                value_cts = pd.Series(my_df[flag])
                for entry in value_cts.keys():
                    layout.template[np.where(layout.structure == entry)] = value_cts[entry]
                sns.heatmap(data=pd.DataFrame(layout.template, index=layout.yticks, columns=layout.xticks),
                            xticklabels="auto", yticklabels="auto",
                            square=True,
                            cbar_kws={"orientation": "horizontal"},
                            cmap=cmap,
                            linewidths=0.20,
                            ax=axs[i])
                axs[i].set_title("{} reads output per channel".format(flag))
        else:
            fig, ax = plt.subplots()

            maxval = max(df_pass.index)  # maximum channel value
            layout = FastqPlots.make_layout(maxval=maxval)
            value_cts = pd.Series(df_pass['Pass'])
            for entry in value_cts.keys():
                layout.template[np.where(layout.structure == entry)] = value_cts[entry]
            sns.heatmap(data=pd.DataFrame(layout.template, index=layout.yticks, columns=layout.xticks),
                        xticklabels="auto", yticklabels="auto",
                        square=True,
                        cbar_kws={"orientation": "horizontal"},
                        cmap='Blues',
                        linewidths=0.20)
            ax.set_title("Pass reads output per channel")

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/channel_output_all.png")

    @staticmethod
    def plot_gc_vs_time(d, out):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots(figsize=(10, 6))

        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                #         average_phred = round(average_phred_full, 1)
                ts_pass.append(tuple((seq.time_string, round(seq.gc, 1))))
            else:
                ts_fail.append(tuple((seq.time_string, round(seq.gc, 1))))

        ts_zero_pass = None
        ts_zero_fail = None
        if ts_pass:
            ts_zero_pass = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_fail:
            ts_zero_fail = min(ts_fail, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        if ts_pass and ts_fail:
            ts_zero = min(ts_zero_pass, ts_zero_fail)
        elif ts_pass:
            ts_zero = ts_zero_pass
        else:  # elif ts_fail:
            ts_zero = ts_zero_fail

        ts_pass3 = list()
        ts_fail3 = list()
        if ts_pass:
            ts_pass1 = [tuple(((x - ts_zero), y)) for x, y in ts_pass]  # Subtract t_zero for the all time points
            ts_pass1.sort(key=lambda x: x[0])  # Sort according to first element in tuple
            ts_pass2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_pass1]  # Convert to hours (float)
            ts_pass3 = [tuple((int(np.round(x)), y)) for x, y in ts_pass2]  # Round hours

            # Convert to dataframe
            df_pass = pd.DataFrame(list(ts_pass3), columns=['Sequencing time interval (h)', '%GC'])
            df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        if ts_fail:
            ts_fail1 = [tuple(((x - ts_zero), y)) for x, y in ts_fail]
            ts_fail1.sort(key=lambda x: x[0])
            ts_fail2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_fail1]
            ts_fail3 = [tuple((int(np.round(x)), y)) for x, y in ts_fail2]

            df_fail = pd.DataFrame(list(ts_fail3), columns=['Sequencing time interval (h)', '%GC'])
            df_fail['Flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

        # Account if there is no fail data or no pass data
        if ts_fail3 and ts_pass3:
            frames = [df_pass, df_fail]
            data = pd.concat(frames)  # Merge dataframes
        elif ts_pass3:
            data = df_pass
        else:  # elif ts_fail3:
            data = df_fail

        # Account if there is no fail data or no pass data
        if ts_fail3 and ts_pass3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='%GC', data=data, hue='Flag',
                               split=True, inner=None)
            g.figure.suptitle('Sequence quality over time')
        elif ts_pass3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='%GC', data=data, inner=None)
            g.figure.suptitle('Sequence quality over time (pass only)')
        else:  # elif ts_fail3:
            g = sns.violinplot(x='Sequencing time interval (h)', y='%GC', data=data, inner=None)
            g.figure.suptitle('Sequence quality over time (fail only)')

        # Major ticks every 4 hours
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        # https://matplotlib.org/2.0.2/examples/ticks_and_spines/tick-locators.html
        def my_formater(val, pos):
            val_str = '{}-{}'.format(int(val), int(val + 1))
            return val_str

        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4))

        if ts_fail:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/gc_vs_time.png")

    @staticmethod
    def plot_gc_vs_length_hex(d, out):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        sns.set(style="ticks")

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [seq.length, seq.gc, seq.flag]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length (bp)', '%GC', 'flag'])
        df_pass = df.loc[df['flag'] == 'pass']
        df_fail = df.loc[df['flag'] == 'fail']

        # Set x-axis limits
        min_len = min(df['Length (bp)'])
        max_len = max(df['Length (bp)'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set y-axis limits
        min_phred = min(df['%GC'])
        max_phred = max(df['%GC'])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Create grid object
        g = sns.JointGrid(x='Length (bp)', y='%GC', data=df, space=0)

        # Plot Fail fist
        if not df_fail.empty:
            x = df_fail['Length (bp)']
            y = df_fail['%GC']
            g.x = df_fail['Length (bp)']
            g.y = df_fail['%GC']

            g.ax_joint.hexbin(x, y, gridsize=50, cmap="Reds", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
            g.ax_marg_x.hist(x, histtype='stepfilled', color='red', alpha=0.6, bins=len_logbins)
            g.ax_marg_y.hist(y, histtype='stepfilled', color='red', alpha=0.6, bins=phred_bins,
                             orientation="horizontal")

        # Plot Pass second
        x = df_pass['Length (bp)']
        y = df_pass['%GC']
        g.x = df_pass['Length (bp)']
        g.y = df_pass['%GC']

        g.ax_joint.hexbin(x, y, gridsize=50, cmap="Blues", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
        g.ax_joint.axis([min_value, max_value, min_phred, max_phred])
        g.ax_marg_x.hist(x, histtype='stepfilled', color='blue', alpha=0.6, bins=len_logbins)
        g.ax_marg_y.hist(y, histtype='stepfilled', color='blue', alpha=0.6, bins=phred_bins, orientation="horizontal")

        # Set main plot x axis scale to log
        g.ax_joint.set_xscale('log')
        g.ax_marg_x.set_xscale('log')
        g.ax_joint.set_xlim((min_value, max_value))

        # Add legend to the joint plot area
        # https://matplotlib.org/tutorials/intermediate/legend_guide.html
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='Pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='Fail')
        if not df_fail.empty:
            g.ax_joint.legend(handles=[blue_patch, red_patch], loc='best')
        else:
            g.ax_joint.legend(handles=[blue_patch], loc='best')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Save figure to file
        g.savefig(out + "/gc_vs_length_hex.png")

    @staticmethod
    def plot_pores_gc_output_vs_time_all(d, out):

        fig, ax = plt.subplots()

        time_list_all = list()
        time_list_pass = list()
        time_list_fail = list()

        for seq_id, seq in d.items():
            time_string = seq.time_string
            time_list_all.append(tuple((time_string, seq.gc)))
            if seq.flag == 'pass':
                time_list_pass.append(tuple((time_string, seq.gc)))
            else:
                time_list_fail.append(tuple((time_string, seq.gc)))

        time_zero = min(time_list_all, key=lambda x: x[0])[0]  # looking for min of 1st elements of list of tuples

        # Compute x_bins
        time_list_all.sort(key=lambda x: x[0])  # order list
        time_list_all[:] = [tuple((x - time_zero, y)) for x, y in time_list_all]  # Subtract t_zero
        time_list_all[:] = [tuple((x.days * 24 + x.seconds / 3600, y)) for x, y in time_list_all]  # Convert to minutes
        # time_list_all[:] = [tuple((int(np.round(x)), y)) for x, y in time_list_all]  # Round minutes
        x = [x for x, y in time_list_all]
        y = [y for x, y in time_list_all]

        # How many bins to plot data
        nbins = int(max(x)) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(x), max(x), nbins)

        # Plot pass - Assume always pass reads present
        time_list_pass.sort(key=lambda x: x[0])  # order list
        time_list_pass[:] = [tuple((x - time_zero, y)) for x, y in time_list_pass]  # Subtract t_zero
        time_list_pass[:] = [tuple((x.days * 24 + x.seconds / 3600, y)) for x, y in time_list_pass]  # Convert to min
        x = [x for x, y in time_list_pass]
        y = [y for x, y in time_list_pass]

        sns.regplot(x=x, y=y, x_bins=x_bins, fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='Pass', color='blue')

        # sns.scatterplot(x=x, y=y, x_bins=x_bins, alpha=0.8, cmap='Blues')

        if time_list_fail:
            time_list_fail.sort(key=lambda x: x[0])  # order list
            time_list_fail[:] = [tuple((x - time_zero, y)) for x, y in time_list_fail]  # Subtract t_zero
            time_list_fail[:] = [tuple((x.days * 24 + x.seconds / 3660, y)) for x, y in time_list_fail]
            x = [x for x, y in time_list_fail]
            y = [y for x, y in time_list_fail]

            # Plot the data
            sns.regplot(x=x, y=y, x_bins=x_bins, fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                        label='Fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('%GC over time')
        plt.ylabel('%GC')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_gc_output_vs_time_all.png")

    @staticmethod
    def plot_pores_length_output_vs_time_all(d, out):
        # Fetch and prepare data from dictionary
        my_dict = defaultdict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = (seq.time_string, seq.length, seq.flag)

        # convert dictionary to pandas dataframe
        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['time_string', 'length', 'flag'])

        # convert datatime to elapsed hours
        time_zero = min(df['time_string'])  # looking for min of 1st elements of list of tuples
        df['time_string'] = df['time_string'] - time_zero
        df['time_string'] = df['time_string'].dt.total_seconds() / 3600

        # Compute x_bins
        nbins = int(max(df['time_string'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['time_string']), max(df['time_string']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['flag'] == 'pass']
        sns.regplot(x=pass_df['time_string'], y=pass_df['length'], x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='Pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['flag'] == 'fail']
        if not fail_df.empty:
            sns.regplot(x=fail_df['time_string'], y=fail_df['length'], x_bins=x_bins,
                        fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                        label='Fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('Read length over time')
        plt.ylabel('Average reads length per 15 min.')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_length_output_vs_time_all.png")

    @staticmethod
    def plot_pores_qual_output_vs_time_all(d, out):
        # Fetch and prepare data from dictionary
        my_dict = defaultdict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = (seq.time_string, seq.average_phred, seq.flag)

        # convert dictionary to pandas dataframe
        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['time_string', 'average_phred', 'flag'])

        # convert datatime to elapsed hours
        time_zero = min(df['time_string'])  # looking for min of 1st elements of list of tuples
        df['time_string'] = df['time_string'] - time_zero
        df['time_string'] = df['time_string'].dt.total_seconds() / 3600

        # Compute x_bins
        nbins = int(max(df['time_string'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['time_string']), max(df['time_string']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['flag'] == 'pass']
        sns.regplot(x=pass_df['time_string'], y=pass_df['average_phred'], x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='Pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['flag'] == 'fail']
        if not fail_df.empty:
            sns.regplot(x=fail_df['time_string'], y=fail_df['average_phred'], x_bins=x_bins,
                        fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                        label='Fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('Read quality over time')
        plt.ylabel('Average reads quality per 15 min.')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_qual_output_vs_time_all.png")

    @staticmethod
    def plot_pores_gc_output_vs_time_per_sample(d, out):
        my_dict = defaultdict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = (seq.time_string, seq.gc, seq.flag, seq.name)

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['time_string', '%GC', 'flag', 'name'])

        # convert datatime to elapsed hours
        time_zero = min(df['time_string'])  # looking for min of 1st elements of list of tuples
        df['time_string'] = df['time_string'] - time_zero
        df['time_string'] = df['time_string'].dt.total_seconds() / 3600

        # Compute x_bins
        nbins = int(max(df['time_string'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['time_string']), max(df['time_string']), nbins)  # Create the bin boundaries

        sample_list = sorted((df['name'].unique()))
        n_sample = len(sample_list)
        width, height = FastqPlots.find_best_matrix(n_sample)
        # print(n_sample, width, height)  # debug

        # Make grid for all samples
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
        fig, ax = plt.subplots(height, width, sharex='col', sharey='row', figsize=(height*5, width*5))
        sample_index = 0
        for i in range(height):
            for j in range(width):
                if sample_index >= len(sample_list):
                    ax[i, j].axis('off')  # don't draw the plot is no more sample for the 'too big' matrix
                else:
                    sample_name = sample_list[sample_index]
                    tmp_df = df[df['name'].str.match(sample_name)]
                    # Pass
                    pass_df = tmp_df[tmp_df['flag'].str.match('pass')]
                    sns.regplot(x=pass_df['time_string'], y=pass_df['%GC'], x_bins=x_bins, ax=ax[i, j],
                                fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                                label='Pass', color='blue')  # capsize=4, capthick=1
                    # Fail
                    fail_df = tmp_df[tmp_df['flag'].str.match('fail')]
                    if not fail_df.empty:
                        sns.regplot(x=fail_df['time_string'], y=fail_df['%GC'], x_bins=x_bins, ax=ax[i, j],
                                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                                    label='Fail', color='red')

                    # Set major ticks every 4 h
                    ax[i, j].xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

                    # Add sample name to graph
                    ax[i, j].set_title(sample_name)
                    ax[i, j].set_xlabel(None)
                    ax[i, j].set_ylabel(None)

                # Move to next sample
                sample_index += 1

        # Add label to axes
        fig.suptitle('%GC over time per sample', fontsize=24)

        # Set common x and y labels
        fig.text(0.5, 0.01, 'Sequencing time (hours)', horizontalalignment='center', verticalalignment='center')
        fig.text(0.01, 0.5, '%GC', horizontalalignment='center', verticalalignment='center',
                 rotation='vertical')

        # Create legend without duplicates
        # https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.figlegend(by_label.values(), by_label.keys())

        plt.tight_layout(rect=[0.02, 0.02, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/pores_gc_output_vs_time_per_sample.png")

    @staticmethod
    def plot_gc_vs_qual_vs_time_3D(d, out):
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        my_dict = defaultdict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = (seq.gc, seq.average_phred, seq.flag, seq.time_string)

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['%GC', 'Phred score', 'flag', 'time_string'])

        # convert datatime to elapsed hours
        time_zero = min(df['time_string'])  # looking for min of 1st elements of list of tuples
        df['time_string'] = df['time_string'] - time_zero
        df['time_string'] = df['time_string'].dt.total_seconds() / 3600

        df_pass = df[df['flag'] == 'pass']
        df_fail = df[df['flag'] == 'fail']

        # g = sns.regplot(x=df_pass['%GC'], y=df_pass['Phred score'], scatter=True,
        #                 scatter_kws={'s': 0.5, 'alpha': 0.01}, label='Pass', color='blue')
        # if not df_fail.empty:
        #     sns.regplot(x=df_fail['%GC'], y=df_fail['Phred score'], scatter=True,
        #                 scatter_kws={'s': 0.5, 'alpha': 0.01}, label='Fail', color='red')

        ax.scatter(df_pass['time_string'].tolist(), df_pass['%GC'].tolist(), df_pass['Phred score'].tolist(),
                   c='blue', s=0.5, alpha=0.01)
        if not df_fail.empty:
            ax.scatter(df_fail['time_string'].tolist(), df_fail['%GC'].tolist(), df_fail['Phred score'].tolist(),
                       c='red', s=0.5, alpha=0.01)

        # ax.legend()
        ax.set_xlabel('Sequencing time (h)')
        ax.set_ylabel('%GC')
        ax.set_zlabel('Phred score')
        ax.set_title('Correlation between %GC and Phred score over time')
        ax.view_init(60, 35)
        # plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/gc_vs_qual_vs_time_3D.png")
        # plt.show()

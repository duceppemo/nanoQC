from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter, MultipleLocator
import matplotlib.patches as mpatches
from dateutil.parser import parse
import functions


class Layout(object):
    def __init__(self, structure, template, xticks, yticks):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks


class SummaryPlots(object):
    @staticmethod
    def make_layout(maxval):
        """
        Make the physical layout of the MinION flowcell.
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
    def get_percentage(pct, allvals):
        absolute = int(pct / 100. * np.sum(allvals))
        return "({:.1f}%, {:d})".format(pct, absolute)

    @staticmethod
    def plot_pores_length_output_vs_time_all_summary(d, out):
        # Fetch and prepare data from dictionary
        my_dict = defaultdict()
        for seq_id, seq in d.items():
            try:
                my_dict[seq_id] = (float(seq.time_stamp), int(seq.length), seq.flag)
            except:
                t =1

        # convert dictionary to pandas dataframe
        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['time_string', 'length', 'flag'])

        # convert datatime to elapsed hours
        df['time_string'] = df['time_string'] / 3600

        # Compute x_bins
        nbins = int(max(df['time_string'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['time_string']), max(df['time_string']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['flag'] == b'TRUE']
        sns.regplot(x=pass_df['time_string'], y=pass_df['length'], x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='Pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['flag'] == b'FALSE']
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
    def plot_pores_qual_output_vs_time_all_summary(d, out):
        # Fetch and prepare data from dictionary
        my_dict = defaultdict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = (float(seq.time_stamp), float(seq.average_phred), seq.flag)

        # convert dictionary to pandas dataframe
        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['time_string', 'average_phred', 'flag'])

        # convert datatime to elapsed hours
        df['time_string'] = df['time_string'] / 3600

        # Compute x_bins
        nbins = int(max(df['time_string'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['time_string']), max(df['time_string']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['flag'] == b'TRUE']
        sns.regplot(x=pass_df['time_string'], y=pass_df['average_phred'], x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='Pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['flag'] == b'FALSE']
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
    def plot_total_reads_vs_time_summary(d, out):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :param out: path to output png file
        :return:
        TODO -> use numpy to handle the plot data, on row per sample?
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()
        t_pass = list()  # time
        t_fail = list()

        for seq_id, seq in d.items():
            t = seq.time_stamp
            if seq.flag == b'TRUE':
                t_pass.append(t)
            else:
                t_fail.append(t)

        # Prepare datetime value for plotting
        # Convert time object in hours from beginning of run
        y_pass = None
        y_fail = None
        if t_pass:
            t_pass[:] = [float(x) / 3600 for x in t_pass]  # Convert to hours (float)
            t_pass.sort()  # Sort
            y_pass = range(1, len(t_pass) + 1, 1)  # Create range. 1 time point equals 1 read

        if t_fail:
            t_fail[:] = [float(x) / 3600 for x in t_fail]
            t_fail.sort()
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
    def plot_reads_per_sample_vs_time_summary(d, out):
        """
        Plot yield per sample. Just the pass reads
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        # fig, ax = plt.subplots()
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == b'TRUE':
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                my_sample_dict[seq.name][seq_id] = seq.time_stamp

        # Order the dictionary by keys
        od = OrderedDict(sorted(my_sample_dict.items()))

        # Make the plot
        legend_names = list()
        for name, seq_ids in od.items():
            name = name.decode('ascii')
            legend_names.append(name)
            ts_pass = list()
            for seq, time_tag in seq_ids.items():
                ts_pass.append(time_tag)

            ts_pass[:] = [float(x) / 3600 for x in ts_pass]  # Convert to hours (float)
            ts_pass.sort()  # Sort
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
    def plot_reads_per_sample_pie_summary(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

        # Fetch required information
        my_sample_dict = defaultdict(list)
        for seq_id, seq in d.items():
            my_sample_dict[seq_id] = [seq.name.decode('ascii'),
                                      seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_sample_dict, orient='index', columns=['name', 'flag'])
        df_all = df.groupby(['name']).count()
        df_pass = df[df['flag'] == 'TRUE'].groupby(['name']).count()
        df_fail = df[df['flag'] == 'FALSE'].groupby(['name']).count()

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

        # Add label to axes
        plt.subplots_adjust(hspace=1)
        fig.suptitle('Distribution of reads among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/reads_per_sample_pie.png")

    @staticmethod
    def plot_bp_per_sample_vs_time_summary(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == b'TRUE':
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                my_sample_dict[seq.name][seq_id] = (seq.time_stamp, seq.length)  # tuple

        # Order the dictionary by keys
        od = OrderedDict(sorted(my_sample_dict.items()))

        # Make the plot
        for name, seq_ids in od.items():
            name = name.decode('ascii')
            ts_pass = list()
            for seq, data_tuple in seq_ids.items():
                ts_pass.append(data_tuple)

            # Prepare x values (time)
            ts_pass[:] = [tuple(((float(x) / 3600), int(y))) for x, y in ts_pass]  # Convert to hours (float)
            ts_pass.sort(key=lambda x: x[0])  # Sort according to first element in tuples
            x_values = [x for x, y in ts_pass]  # Only get the fist value of the ordered tuples

            # Prepare y values in a cumulative way
            c = 0
            y_values = list()
            for x, y in ts_pass:
                y = y + c
                y_values.append(y)
                c = y

            # Plot values per sample
            ax.plot(x_values, y_values,
                    label="%s (%s)" % (name, "{:,}".format(max(y_values))))

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # New
        # Add axes labels
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Yield per sample in base pair\n("pass" only)')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        # Save figure to file
        fig.savefig(out + "/bp_per_sample_vs_time.png")

    @staticmethod
    def plot_bp_per_sample_pie_summary(d, out):
        """
        Read length per sample vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

        # Fetch required information
        my_sample_dict = defaultdict(list)
        for seq_id, seq in d.items():
            my_sample_dict[seq_id] = [seq.name.decode('ascii'),
                                      int(seq.length.decode('ascii')),
                                      seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_sample_dict, orient='index', columns=['name', 'length', 'flag'])
        df_all = df.groupby(['name']).sum()
        df_pass = df[df['flag'] == 'TRUE'].groupby(['name']).sum()
        df_fail = df[df['flag'] == 'FALSE'].groupby(['name']).sum()

        # Make the plots
        titles = ['All', 'Pass', 'Fail']
        for i, my_df in enumerate([df_all, df_pass, df_fail]):
            my_df = my_df.sort_values(['length'], ascending=False)  # sort dataframe for better looking pie chart
            data = list(my_df['length'])
            data_sum = sum(data)
            labels = list(my_df.index)
            # label2 = [NanoQC.get_percentage(pct, data) for pct in data]
            for j, l in enumerate(labels):
                labels[j] = "{:s} ({:,} bp, {:.1f}%)".format(l, data[j], round(data[j] / data_sum * 100, 1))

            axs[i].pie(data, labels=labels, wedgeprops={'linewidth': 2, 'edgecolor': 'w'})
            axs[i].set_title(titles[i])

        # Add label to axes
        plt.subplots_adjust(hspace=0.5)
        fig.suptitle('Distribution of base pairs among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/bp_per_sample_pie.png")

    @staticmethod
    def plot_total_bp_vs_time_summary(d, out):
        """
        Sequence length vs time
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()

        # Fetch required information
        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == b'TRUE':
                ts_pass.append(tuple((seq.time_stamp, seq.length)))
            else:
                ts_fail.append(tuple((seq.time_stamp, seq.length)))

        x_pass_values = None
        y_pass_values = None
        x_fail_values = None
        y_fail_values = None
        if ts_pass:
            ts_pass[:] = [tuple(((float(x) / 3600), int(y))) for x, y in ts_pass]  # Convert to hours (float)
            ts_pass.sort(key=lambda x: x[0])  # Sort according to first element in tuple
            x_pass_values = [x for x, y in ts_pass]
            c = 0
            y_pass_values = list()
            for x, y in ts_pass:
                y = y + c
                y_pass_values.append(y)
                c = y
        if ts_fail:
            ts_fail[:] = [tuple(((float(x) / 3600), int(y))) for x, y in ts_fail]
            ts_fail.sort(key=lambda x: x[0])
            x_fail_values = [x for x, y in ts_fail]
            c = 0
            y_fail_values = list()
            for x, y in ts_fail:
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
    def plot_quality_vs_time_summary(d, out):
            """
            Quality vs time (bins of 1h). Violin plot
            :param d: Dictionary
            :param out: path to output png file
            :return:
            name, length, flag, average_phred, gc, time_string
            """

            fig, ax = plt.subplots(figsize=(10, 4))


            ################################################
            # Fetch required information
            my_sample_dict = defaultdict()
            for seq_id, seq in d.items():
                if seq.name not in my_sample_dict:
                    my_sample_dict[seq.name] = defaultdict()
                my_sample_dict[seq.name][seq_id] = (seq.time_stamp,  # hours
                                                    round(float(seq.average_phred), 2),
                                                    seq.flag)  # time in h

            # Order the dictionary by keys
            od = OrderedDict(sorted(my_sample_dict.items()))

            # Make the plot
            ts_pass = list()
            ts_fail = list()
            for name, sub_dict in od.items():
                for seq_id, data_tuple in sub_dict.items():
                    if data_tuple[2] == b'TRUE':
                        ts_pass.append(data_tuple)
                    else:
                        ts_fail.append(data_tuple)

            # Prepare x values (time)
            ts_pass[:] = [tuple(((float(x) / 3600), int(y))) for x, y, z in ts_pass]  # Convert to hours (float)
            ts_pass.sort(key=lambda x: x[0])  # Sort according to first element in tuples
            ts_pass[:] = [tuple((int(np.round(x)), y)) for x, y in ts_pass]  # Round hours
            # x_values_pass = [x for x, y in ts_pass]  # Only get the fist value of the ordered tuples

            ts_fail[:] = [tuple(((float(x) / 3600), int(y))) for x, y, z in ts_fail]  # Convert to hours (float)
            ts_fail.sort(key=lambda x: x[0])  # Sort according to first element in tuples
            ts_fail[:] = [tuple((int(np.round(x)), y)) for x, y in ts_fail]  # Round hours
            # x_values_fail = [x for x, y in ts_fail]  # Only get the fist value of the ordered tuples

            ################################################

            if ts_pass:
                df_pass = pd.DataFrame(list(ts_pass),
                                       columns=['Sequencing time interval (h)', 'Phred score'])  # Convert to dataframe
                df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

            if ts_fail:
                df_fail = pd.DataFrame(list(ts_fail), columns=['Sequencing time interval (h)', 'Phred score'])
                df_fail['Flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value

            # Account if there is no fail data or no pass data
            if ts_fail and ts_pass:
                frames = [df_pass, df_fail]
                data = pd.concat(frames)  # Merge dataframes
            elif ts_pass:
                data = df_pass
            else:  # elif ts_fail3:
                data = df_fail

            # Account if there is no fail data or no pass data
            if ts_fail and ts_pass:
                g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, hue='Flag', split=True,
                                   inner=None)
                g.figure.suptitle('Sequence quality over time')
            elif ts_pass:
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
    def plot_phred_score_distribution_summary(d, out):
        """
        Frequency of phred scores
        :param d: Dictionary
        :param out: path to output png file
        :return:
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [float(seq.average_phred), seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'TRUE']
        df_fail = df.loc[df['flag'] == 'FALSE']

        average_qual_pass = functions.compute_average_quality(df_pass['Phred score'].tolist(), df_pass.shape[0])
        average_qual_fail = functions.compute_average_quality(df_fail['Phred score'].tolist(), df_fail.shape[0])

        ax.hist(df_pass['Phred score'], histtype='stepfilled', color='blue', alpha=0.6,
                label='Pass (Avg: {})'.format(round(average_qual_pass, 2)))

        if not df_fail.empty:
            ax.hist(df_fail['Phred score'], histtype='stepfilled', color='red', alpha=0.6,
                    label='Fail (Avg: {})'.format(round(average_qual_fail, 2)))

        plt.legend()

        ax.set(xlabel='Phred score', ylabel='Frequency', title='Phred score distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/phred_score_distribution.png")

    @staticmethod
    def plot_length_distribution_summary(d, out):
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
            my_dict[seq_id] = [float(seq.length), seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length', 'flag'])
        df_pass = df.loc[df['flag'] == 'TRUE']
        df_fail = df.loc[df['flag'] == 'FALSE']

        # Set bin sized for histogram
        min_len = min(df['Length'])
        max_len = max(df['Length'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        len_logbins = np.logspace(min_exp, max_exp, 50)

        average_len_pass = int(df_pass['Length'].mean())
        average_len_fail = int(df_fail['Length'].mean())

        ax.hist(df_pass['Length'], histtype='stepfilled', color='blue', alpha=0.6,
                label="Pass (Avg: {} bp)".format(average_len_pass), bins=len_logbins)

        if not df_fail.empty:
            ax.hist(df_fail['Length'], histtype='stepfilled', color='red', alpha=0.6,
                    label="Fail (Avg: {} bp)".format(average_len_fail), bins=len_logbins)

        plt.legend()
        plt.xscale('log')
        ax.set(xlabel='Read length (bp)', ylabel='Frequency', title='Read length distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/length_distribution.png")

    @staticmethod
    def plot_quality_vs_length_hex_summary(d, out):
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
            my_dict[seq_id] = [int(seq.length), float(seq.average_phred), seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length (bp)', 'Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'TRUE']
        df_fail = df.loc[df['flag'] == 'FALSE']

        # Set x-axis limits
        min_len = min(df['Length (bp)'])
        max_len = max(df['Length (bp)'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 50)

        # Set y-axis limits
        min_phred = min(df['Phred score'])
        max_phred = max(df['Phred score'])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 30)

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
    def plot_reads_vs_bp_per_sample_summary(d, out):
        # Fetch required information
        my_sample_dict = defaultdict()  # to get the lengths (bp)
        for seq_id, seq in d.items():
            if seq.flag == b'TRUE':
                length = int(seq.length)
                if length == 0:
                    continue
                name = seq.name.decode('ascii')
                if name not in my_sample_dict:
                    my_sample_dict[name] = defaultdict()
                if not my_sample_dict[name]:
                    my_sample_dict[name] = [length]
                else:
                    my_sample_dict[name].append(length)

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

    @staticmethod
    def plot_pores_output_vs_time_summary(d, out):

        time_list = list()
        for seq_id, seq in d.items():
            time_list.append(seq.time_stamp)

        time_list[:] = [float(x) / 60 for x in time_list]  # Convert to minutes (float)
        time_list[:] = [int(np.round(x)) for x in time_list]  # Round minutes
        # Check how many 15-minute bins are required to plot all the data
        nbins = max(time_list) / 15 if max(time_list) % 15 == 0 else int(max(time_list) / 15) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(time_list), max(time_list), nbins)  # every 15 min

        # Generate counts for each bin
        hist, edges = np.histogram(time_list, bins=x_bins, density=False)

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
    def plot_pores_output_vs_time_all_summary(d, out):

        import matplotlib.lines as mlines

        fig, ax = plt.subplots()

        time_list_all = list()
        time_list_pass = list()
        time_list_fail = list()

        for seq_id, seq in d.items():
            time_stamp = seq.time_stamp
            time_list_all.append(time_stamp)
            if seq.flag == b'TRUE':
                time_list_pass.append(time_stamp)
            else:
                time_list_fail.append(time_stamp)

        time_zero = min(time_list_all)  # find smallest datetime value

        # Plot all
        time_list_all[:] = [float(x) / 60 for x in time_list_all]  # Convert to minutes (float)
        time_list_all[:] = [int(np.round(x)) for x in time_list_all]  # Round minutes
        # Check how many 15-minute bins are required to plot all the data
        nbins = max(time_list_all) / 15 if max(time_list_all) % 15 == 0 else int(max(time_list_all) / 15) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(time_list_all), max(time_list_all), nbins)  # every 15 min

        # If not fail, just draw the pass. Else, draw total, fail and pass
        if time_list_fail:  # fail reads might be missing if plotting filtered reads for example.
            # Plot pass - Assume always pass reads present
            time_list_pass[:] = [float(x) / 60 for x in time_list_pass]
            time_list_pass[:] = [int(np.round(x)) for x in time_list_pass]
            hist, edges = np.histogram(time_list_pass, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0, cmap='Blues')

            # Plot fail
            time_list_fail[:] = [float(x) / 60 for x in time_list_fail]
            time_list_fail[:] = [int(np.round(x)) for x in time_list_fail]
            hist, edges = np.histogram(time_list_fail, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0, cmap='Reds')

        # Plot all
        hist, edges = np.histogram(time_list_all, bins=x_bins, density=False)
        sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0, cmap='Greens')

        # Adjust format of numbers for y axis: "1000000" -> "1,000,000"
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

        # Change x axis labels chunk-of-15-min to hours
        def numfmt(m, pos):
            h = '{}'.format(m / 4)
            return h
        ax.xaxis.set_major_formatter(FuncFormatter(numfmt))

        # Major ticks every 4 hours
        def my_formater(val, pos):
            val_str = '{}'.format(int(val/4))
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
    def plot_channel_output_total(d, out):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :param out: path to output png file
        :return:
        """

        channel_dict = defaultdict()
        for seq_id, seq in d.items():
            length = int(seq.length)
            if length == 0:
                continue
            channel_number = int(seq.channel)
            # Pass and fail appart
            # if channel_number not in channel_dict:
            #     channel_dict[channel_number] = [0, 0]
            # if seq.flag == b'TRUE':
            #     channel_dict[channel_number][0] += 1
            # else:  # seq.flag == b'FALSE':
            #     channel_dict[channel_number][1] += 1

            # Pass and fail together
            if channel_number not in channel_dict:
                channel_dict[channel_number] = 0
            channel_dict[channel_number] += 1

        # convert to Pandas dataframe
        df = pd.DataFrame.from_dict(channel_dict, orient='index', columns=['Reads'])
        # Sort df on value
        # df = df.sort_values(by='Reads', ascending=True)

        # Plot

        maxval = max(df.index)
        layout = SummaryPlots.make_layout(maxval=maxval)
        value_cts = pd.Series(df['Reads'])
        for entry in value_cts.keys():
            layout.template[np.where(layout.structure == entry)] = value_cts[entry]
        plt.figure()
        ax = sns.heatmap(data=pd.DataFrame(layout.template, index=layout.yticks, columns=layout.xticks),
                         xticklabels="auto", yticklabels="auto",
                         square=True,
                         cbar_kws={"orientation": "horizontal"},
                         cmap='Greens',
                         linewidths=0.20)
        plt.tight_layout()  # Get rid of extra margins around the plot
        fig = ax.get_figure()  # Get figure from FacetGrid
        fig.savefig(out + "/channel_output_total.png")

    @staticmethod
    def plot_channel_output_pass_fail(d, out):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :param out: path to output png file
        :return:
        """

        channel_dict = defaultdict()
        for seq_id, seq in d.items():
            length = int(seq.length)
            if length == 0:
                continue
            channel_number = int(seq.channel)
            # Pass and fail apart
            if channel_number not in channel_dict:
                channel_dict[channel_number] = [0, 0]
            if seq.flag == b'TRUE':
                channel_dict[channel_number][0] += 1
            else:  # seq.flag == b'FALSE':
                channel_dict[channel_number][1] += 1

        # convert to Pandas dataframe
        df = pd.DataFrame.from_dict(channel_dict, orient='index', columns=['Pass', 'Fail'])
        df_pass = df[['Pass']]  # The double square brakets keep the column name
        df_fail = df[['Fail']]

        # Plot
        fig, axs = plt.subplots(nrows=2, figsize=(6, 8))

        for i, my_tuple in enumerate([(df_pass, 'Pass', 'Blues'), (df_fail, 'Fail', 'Reds')]):
            my_df = my_tuple[0]
            flag = my_tuple[1]
            cmap = my_tuple[2]
            maxval = max(my_df.index)  # maximum channel value
            layout = SummaryPlots.make_layout(maxval=maxval)
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
        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/channel_output_pass_fail.png")

    @staticmethod
    def plot_channel_output_all_summary(d, out):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :param out: path to output png file
        :return:
        """

        channel_dict = defaultdict()
        for seq_id, seq in d.items():
            length = int(seq.length)
            if length == 0:
                continue
            channel_number = int(seq.channel)
            # Pass and fail apart
            if channel_number not in channel_dict:
                channel_dict[channel_number] = [0, 0]
            if seq.flag == b'TRUE':
                channel_dict[channel_number][0] += 1
            else:  # seq.flag == b'FALSE':
                channel_dict[channel_number][1] += 1

        # convert to Pandas dataframe
        df = pd.DataFrame.from_dict(channel_dict, orient='index', columns=['Pass', 'Fail'])
        df_all = pd.DataFrame()
        df_all['All'] = df['Pass'] + df['Fail']
        df_pass = df[['Pass']]  # The double square brakets keep the column name
        df_fail = df[['Fail']]

        # Plot
        fig, axs = plt.subplots(nrows=3, figsize=(6, 12))

        for i, my_tuple in enumerate([(df_all, 'All', 'Greens'),
                                      (df_pass, 'Pass', 'Blues'),
                                      (df_fail, 'Fail', 'Reds')]):
            my_df = my_tuple[0]
            flag = my_tuple[1]
            cmap = my_tuple[2]

            maxval = max(my_df.index)  # maximum channel value
            layout = SummaryPlots.make_layout(maxval=maxval)
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
        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/channel_output_all.png")

    # @staticmethod
    # def plot_quality_vs_time_summary(d, out):
    #     """
    #     Quality vs time (bins of 1h). Violin plot
    #     :param d: Dictionary
    #     :param out: path to output png file
    #     :return:
    #     name, length, flag, average_phred, gc, time_string
    #     """
    #
    #     fig, ax = plt.subplots(figsize=(10, 6))
    #
    #     ts_pass = list()
    #     ts_fail = list()
    #     for seq_id, seq in d.items():
    #         if seq.flag == b'TRUE':
    #             ts_pass.append(tuple((seq.time_stamp, round(float(seq.average_phred), 1))))
    #         else:
    #             ts_fail.append(tuple((seq.time_stamp, round(float(seq.average_phred), 1))))
    #
    #     if ts_pass:
    #         ts_pass[:] = [tuple(((float(x) / 3600), y)) for x, y in ts_pass]  # Convert to hours (float)
    #         ts_pass.sort(key=lambda x: x[0])  # Sort according to first element in tuple
    #         ts_pass[:] = [tuple((int(np.round(x)), y)) for x, y in ts_pass]  # Round hours
    #
    #         df_pass = pd.DataFrame(list(ts_pass), columns=['Sequencing time interval (h)', 'Phred score'])  # Convert to dataframe
    #         df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value
    #
    #     if ts_fail:
    #         ts_fail[:] = [tuple(((float(x) / 3600), y)) for x, y in ts_fail]
    #         ts_fail.sort(key=lambda x: x[0])
    #         ts_fail[:] = [tuple((int(np.round(x)), y)) for x, y in ts_fail]
    #
    #         df_fail = pd.DataFrame(list(ts_fail), columns=['Sequencing time interval (h)', 'Phred score'])
    #         df_fail['Flag'] = pd.Series('fail', index=df_fail.index)  # Add a 'Flag' column to the end with 'fail' value
    #
    #     # Account if there is no fail data or no pass data
    #     if ts_fail and ts_pass:
    #         frames = [df_pass, df_fail]
    #         data = pd.concat(frames)  # Merge dataframes
    #     elif ts_pass:
    #         data = df_pass
    #     else:  # elif ts_fail3:
    #         data = df_fail
    #
    #     # Account if there is no fail data or no pass data
    #     if ts_fail and ts_pass:
    #         g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, hue='Flag', split=True, inner=None)
    #         g.figure.suptitle('Sequence quality over time')
    #     elif ts_pass:
    #         g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, inner=None)
    #         g.figure.suptitle('Sequence quality over time (pass only)')
    #     else:  # elif ts_fail:
    #         g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, inner=None)
    #         g.figure.suptitle('Sequence quality over time (fail only)')
    #
    #     # Major ticks every 4 hours
    #     # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
    #     # https://matplotlib.org/2.0.2/examples/ticks_and_spines/tick-locators.html
    #     def my_formater(val, pos):
    #         val_str = '{}-{}'.format(int(val), int(val + 1))
    #         return val_str
    #
    #     ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
    #     ax.xaxis.set_major_locator(MultipleLocator(4))
    #
    #     if ts_fail:
    #         plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #
    #     plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
    #     fig.savefig(out + "/channel_output_all.png")

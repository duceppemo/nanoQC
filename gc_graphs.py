from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter, MultipleLocator
import matplotlib.patches as mpatches
from math import sqrt, ceil


class Layout(object):
    def __init__(self, structure, template, xticks, yticks):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks


class GcPlots(object):
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
    def plot_gc_vs_time(df1, out):
        # Get data
        df = df1.loc[:, ('Time', 'Flag', 'GC')]

        # Round %GC to one decimal
        df['GC'] = df['GC'].round(1)

        # convert datatime to elapsed minutes
        time_zero = df.loc[:, 'Time'].min()  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600 / 4  # every hour
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer

        # Make plot
        fig, ax = plt.subplots(figsize=(10, 6))

        # Account if there is no fail data or no pass data
        # g = sns.violinplot(data=df, x='Time', y='GC', hue='Flag',
        #                    split=True, inner=None)
        g = sns.boxplot(data=df, x='Time', y='GC', hue='Flag', palette=['blue', 'red'], showfliers=False)
        g.figure.suptitle('%GC over time')

        # Major ticks every 4 hours
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        # https://matplotlib.org/2.0.2/examples/ticks_and_spines/tick-locators.html
        def my_formater(val, pos):
            val_str = '{}-{}'.format(int(val), int(val + 1))
            return val_str

        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4))

        g.legend_.set_title(None)  # Remove legend title
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/gc_vs_time.png")
        plt.close()

    @staticmethod
    def plot_gc_vs_length_hex(df1, out):

        sns.set(style="ticks")

        df = df1.loc[:, ('GC', 'Flag', 'Length')]
        df_pass = df.loc[df['Flag'] == 'pass']
        df_fail = df.loc[df['Flag'] == 'fail']

        # Set x-axis limits
        min_len = min(df['Length'])
        max_len = max(df['Length'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        min_value = float(10 ** (min_exp - 0.1))
        max_value = float(10 ** (max_exp + 0.1))

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set y-axis limits
        min_phred = min(df['GC'])
        max_phred = max(df['GC'])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Create grid object
        g = sns.JointGrid(x='Length', y='GC', data=df, space=0)

        # Plot Fail fist
        if not df_fail.empty:
            x = df_fail['Length']
            y = df_fail['GC']
            g.x = df_fail['Length']
            g.y = df_fail['GC']

            g.ax_joint.hexbin(x, y, gridsize=50, cmap="Reds", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
            g.ax_marg_x.hist(x, histtype='stepfilled', color='red', alpha=0.6, bins=len_logbins)
            g.ax_marg_y.hist(y, histtype='stepfilled', color='red', alpha=0.6, bins=phred_bins,
                             orientation="horizontal")

        # Plot Pass second
        x = df_pass['Length']
        y = df_pass['GC']
        g.x = df_pass['Length']
        g.y = df_pass['GC']

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
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='fail')
        if not df_fail.empty:
            g.ax_joint.legend(handles=[blue_patch, red_patch], loc='best')
        else:
            g.ax_joint.legend(handles=[blue_patch], loc='best')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Save figure to file
        g.savefig(out + "/gc_vs_length_hex.png")
        plt.close()

    @staticmethod
    def plot_pores_gc_output_vs_time_all(df1, out):
        df = df1.loc[:, ('Time', 'GC', 'Flag')]

        # convert datatime to elapsed minutes
        time_zero = df.loc[:, 'Time'].min()  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3660
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer

        df_pass = df.loc[df['Flag'] == 'pass']
        df_fail = df.loc[df['Flag'] == 'fail']

        # How many bins to plot data
        nbins = int(max(df['Time'])) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min(df['Time']), max(df['Time']), nbins)

        # Make plot
        fig, ax = plt.subplots()

        sns.regplot(data=df_pass, x='Time', y='GC', x_bins=x_bins, fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='pass', color='blue')

        if not df_fail.empty:
            sns.regplot(data=df_fail, x='Time', y='GC', x_bins=x_bins, fit_reg=False,
                        scatter_kws={'alpha': 0.6, 's': 30},
                        label='fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('%GC over time')
        plt.ylabel('%GC')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_gc_output_vs_time_all.png")
        plt.close()

    @staticmethod
    def plot_pores_gc_output_vs_time_per_sample(df1, out):
        df = df1.loc[:, ('Name', 'Time', 'Flag', 'GC')]

        # convert datetime to elapsed hours
        time_zero = min(df.loc[:, 'Time'])  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600

        # Compute x_bins
        nbins = int(max(df['Time'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['Time']), max(df['Time']), nbins)  # Create the bin boundaries

        sample_list = sorted((df['Name'].unique()))
        n_sample = len(sample_list)
        width, height = GcPlots.find_best_matrix(n_sample)
        # print(n_sample, width, height)  # debug

        # Make grid for all samples
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
        fig, ax = plt.subplots(height, width, sharex='col', sharey='row', figsize=(height*5, width*5), squeeze=False)
        sample_index = 0
        for i in range(height):
            for j in range(width):
                if sample_index >= len(sample_list):
                    ax[i, j].axis('off')  # don't draw the plot is no more sample for the 'too big' matrix
                else:
                    sample_name = sample_list[sample_index]
                    tmp_df = df[df['Name'].str.match(sample_name)]
                    # Pass
                    pass_df = tmp_df[tmp_df['Flag'].str.match('pass')]
                    sns.regplot(x=pass_df['Time'], y=pass_df['GC'], x_bins=x_bins, ax=ax[i, j],
                                fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                                label='pass', color='blue')  # capsize=4, capthick=1
                    # Fail
                    fail_df = tmp_df[tmp_df['Flag'].str.match('fail')]
                    if not fail_df.empty:
                        sns.regplot(x=fail_df['Time'], y=fail_df['GC'], x_bins=x_bins, ax=ax[i, j],
                                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                                    label='fail', color='red')

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
        plt.close()

    @staticmethod
    def plot_gc_vs_qual_vs_time_3D(d, out):
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

        my_dict = defaultdict()
        for seq_id, seq in d.items():
            if not seq.time_string:
                return
            my_dict[seq_id] = (seq.gc, seq.average_phred, seq.flag, seq.time_string)

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['%GC', 'Phred score', 'flag', 'time_string'])

        # convert datatime to elapsed hours
        time_zero = min(df['time_string'])  # looking for min of 1st elements of list of tuples
        df['time_string'] = df['time_string'] - time_zero
        df['time_string'] = df['time_string'].dt.total_seconds() / 3600

        df_pass = df[df['flag'] == 'pass']
        df_fail = df[df['flag'] == 'fail']

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # g = sns.regplot(x=df_pass['%GC'], y=df_pass['Phred score'], scatter=True,
        #                 scatter_kws={'s': 0.5, 'alpha': 0.01}, label='pass', color='blue')
        # if not df_fail.empty:
        #     sns.regplot(x=df_fail['%GC'], y=df_fail['Phred score'], scatter=True,
        #                 scatter_kws={'s': 0.5, 'alpha': 0.01}, label='fail', color='red')

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
        plt.close()

    @staticmethod
    def plot_gc_distribution_per_sample(df1, out):
        # Only keep pass reads
        df = df1.loc[:, ('Name', 'Flag', 'GC')]

        # Make plot
        fig, ax = plt.subplots(figsize=(10, 6))

        g = sns.boxplot(data=df, x='Name', y='GC', hue='Flag', palette=['blue', 'red'], showfliers=False)
        g.figure.suptitle('%GC distribution per sample')

        g.legend_.set_title(None)  # Remove legend title
        ax.tick_params(axis='x', rotation=45)  # Rotate tick text on 45 degree

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/gc_distribution_per_sample.png")
        plt.close()

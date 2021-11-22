from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter, MultipleLocator, StrMethodFormatter
import matplotlib.patches as mpatches
from math import sqrt, ceil
import plotly
import plotly.express as px
from datetime import datetime
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
    def plot_total_reads_vs_time_plotly(df1, out):
        # Extract columns of interest from master dataframe
        df = df1.loc[:, ('Time', 'Flag')]

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Time and Flag, count how many reads for each unit of time and transpose the dataframe
        df = df.groupby(['Time', 'Flag'], as_index=False).size().pivot('Time', 'Flag', 'size')

        # Cumulative sum of counts for y values
        df['pass'] = df['pass'].cumsum()
        df['fail'] = df['fail'].cumsum()

        # Transpose back the new dataframe with the cumulative sum values
        df = pd.melt(df.reset_index(), id_vars=['Time'], value_vars=['pass', 'fail'],
                     var_name='Flag', value_name='Count')

        # Make the plot using seaborn
        fig = px.line(df, x='Time', y='Count', color='Flag')
        # Save to file
        plotly.offline.plot(fig, filename=out + "/total_reads_vs_time.html")
        plt.close()

    @staticmethod
    def plot_total_reads_vs_time(df1, out):
        # Extract columns of interest from master dataframe
        df = df1.loc[:, ('Time', 'Flag')]

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Time and Flag, count how many reads for each unit of time and transpose the dataframe
        df = df.groupby(['Time', 'Flag'], as_index=False).size().pivot('Time', 'Flag', 'size')

        # Cumulative sum of counts for y values
        df['pass'] = df['pass'].cumsum()
        df['fail'] = df['fail'].cumsum()

        # Transpose back the new dataframe with the cumulative sum values
        df = pd.melt(df.reset_index(), id_vars=['Time'], value_vars=['pass', 'fail'],
                     var_name='Flag', value_name='Count')

        # Make the plot using seaborn
        fig, ax = plt.subplots()
        sns.set_palette(sns.color_palette(['blue', 'red']))
        g = sns.lineplot(data=df, x='Time', y='Count', hue='Flag')
        # Add a comma for the thousands
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        # Set the axes and figure titles
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Total read yield')
        # Remove legend title:
        g.legend_.set_title(None)

        # Remove extra white space around the figure
        plt.tight_layout()
        # Save to file
        fig.savefig(out + "/total_reads_vs_time.png")
        plt.close()

    @staticmethod
    def plot_reads_per_sample_vs_time(df1, out):
        # Only keep pass reads
        df = df1.loc[(df1['Flag'] == 'pass'), ['Name', 'Flag', 'Time']]

        # Drop the Flag column
        df = df.loc[:, ['Name', 'Time']].reset_index(drop=True)

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Sample and Time, count how many reads for each unit of time
        df = df.groupby(['Name', 'Time'], as_index=False).agg(Count=('Time', 'size'))
        # Add ccumulative sum for each time point
        df['CumSum'] = df.groupby(['Name', 'Time'], as_index=False).sum().groupby('Name')['Count'].cumsum()
        # print('\n{}'.format(df))  # debug

        # Make Plot
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        g = sns.lineplot(data=df, x='Time', y='CumSum', hue='Name')
        # Add a comma for the thousands
        # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        # Set the axes and figure titles
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Pass reads per sample')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        # Remove legend title:
        g.legend_.set_title(None)
        # Remove extra white space around the figure
        plt.tight_layout()
        # Save to file
        fig.savefig(out + "/reads_per_sample_vs_time.png")
        plt.close()

    @staticmethod
    def plot_reads_per_sample_pie(df1, out):
        # Fetch required information
        df = df1.loc[:, ('Name', 'Flag')]

        df_all = df.groupby(['Name']).count()
        df_pass = df[df['Flag'] == 'pass'].groupby(['Name']).count()
        df_fail = df[df['Flag'] == 'fail'].groupby(['Name']).count()

        if not df_fail.empty:
            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

            # Make the plots
            titles = ['all', 'pass', 'fail']
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
            ax.set_title('pass')

        # Add label to axes
        plt.subplots_adjust(hspace=1)
        fig.suptitle('Distribution of reads among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/reads_per_sample_pie.png")
        plt.close()

    @staticmethod
    def plot_bp_per_sample_pie(df1, out):
        # Fetch required information
        df = df1.loc[:, ('Name', 'Flag', 'Length')]

        df_all = df.groupby(['Name']).sum()
        df_pass = df[df['Flag'] == 'pass'].groupby(['Name']).sum()
        df_fail = df[df['Flag'] == 'fail'].groupby(['Name']).sum()

        mpl.style.use('default')

        if not df_fail.empty:
            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 14))

            # Make the plots
            titles = ['all', 'pass', 'fail']
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
            ax.set_title('pass')

        # Add label to axes
        plt.subplots_adjust(hspace=1)
        fig.suptitle('Distribution of bp among samples', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.5)  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/bp_per_sample_pie.png")
        plt.close()

    @staticmethod
    def plot_bp_per_sample_vs_time(df1, out):
        # Only keep pass reads
        df = df1.loc[(df1['Flag'] == 'pass'), ['Name', 'Flag', 'Time', 'Length']]

        # Drop the Flag column
        df = df.loc[:, ['Name', 'Time', 'Length']].reset_index(drop=True)

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Sample and Time, count how many reads for each unit of time
        df = df.groupby(['Name', 'Time'], as_index=False).agg(Sum=('Length', 'sum'))
        # Add cumulative sum for each time point
        df['CumSum'] = df.groupby(['Name', 'Time'], as_index=False).sum().groupby('Name')['Sum'].cumsum()
        # print('\n{}'.format(df))  # debug

        # Make Plot
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        g = sns.lineplot(data=df, x='Time', y='CumSum', hue='Name')
        # Add a comma for the thousands
        # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        # Set the axes and figure titles
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Yield per sample in base pair\n("pass" only)')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        # Remove legend title:
        g.legend_.set_title(None)
        # Remove extra white space around the figure
        plt.tight_layout()
        # Save to file
        fig.savefig(out + "/bp_per_sample_vs_time.png")
        plt.close()

    @staticmethod
    def plot_total_bp_vs_time(df1, out):
        # Extract columns of interest from master dataframe
        df = df1.loc[:, ('Time', 'Flag', 'Length')]

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Time and Flag, count how many reads for each unit of time and transpose the dataframe
        df = df.groupby(['Time', 'Flag'], as_index=False).sum()
        df['cumsum'] = df.groupby(['Time', 'Flag'], as_index=False)['Length'].cumsum()
        df = df.pivot(index='Time', columns='Flag', values='cumsum')

        # Cumulative sum of counts for y values
        df['pass'] = df['pass'].cumsum()
        df['fail'] = df['fail'].cumsum()

        # Transpose back the new dataframe with the cumulative sum values
        df = pd.melt(df.reset_index(), id_vars=['Time'], value_vars=['pass', 'fail'],
                     var_name='Flag', value_name='cumsum')

        # Make the plot usinf seaborn
        fig, ax = plt.subplots()
        sns.set_palette(sns.color_palette(['blue', 'red']))
        g = sns.lineplot(data=df, x='Time', y='cumsum', hue='Flag')
        # Add a comma for the thousands
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        # Set the axes and figure titles
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Total read yield')
        # Remove legend title:
        g.legend_.set_title(None)

        # Remove extra white space around the figure
        plt.tight_layout()
        # Save to file
        fig.savefig(out + "/total_bp_vs_time.png")
        plt.close()

    @staticmethod
    def plot_quality_vs_time(df1, out):
        # Fetch data of interest
        df = df1.loc[:, ('Time', 'Flag', 'Qual')]

        # Find the smallest datetime value
        time_zero = df.loc[:, 'Time'].min()
        # Subtract time zero from all timedate values to get elapsed time
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600 / 4  # convert to every 4 hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600 / 4  # convert to every 4 hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Make plot
        fig, ax = plt.subplots(figsize=(10, 6))

        # Account if there is no fail data or no pass data
        # g = sns.violinplot(x='Time', y='Qual', data=df, hue='Flag', split=True, inner=None)
        g = sns.boxplot(data=df, x='Time', y='Qual', hue='Flag', palette=['blue', 'red'], showfliers=False)
        g.figure.suptitle('Sequence quality over time')
        # Remove legend title:
        g.legend_.set_title(None)

        # Major ticks every 4 hours
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        # https://matplotlib.org/2.0.2/examples/ticks_and_spines/tick-locators.html
        def my_formater(val, pos):
            val_str = '{}-{}'.format(int(val), int(val + 1))
            return val_str

        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4))

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/quality_vs_time.png")
        plt.close()

    @staticmethod
    def plot_phred_score_distribution(df1, out):
        """
        Frequency of phred scores
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        df = df1.loc[:, ('Qual', 'Flag')]
        df_pass = df.loc[df['Flag'] == 'pass']
        df_fail = df.loc[df['Flag'] == 'fail']

        average_qual_pass = round(functions.compute_average_quality(df_pass['Qual'].tolist(),
                                                                    df_pass.shape[0]), 1)

        ax.hist(df_pass['Qual'], histtype='stepfilled', color='blue', alpha=0.6,
                label='pass (avg: {})'.format(average_qual_pass))

        if not df_fail.empty:
            average_qual_fail = round(functions.compute_average_quality(df_fail['Qual'].tolist(),
                                                                        df_fail.shape[0]), 1)
            ax.hist(df_fail['Qual'], histtype='stepfilled', color='red', alpha=0.6,
                    label='fail (avg: {})'.format(average_qual_fail))

        plt.legend()

        ax.set(xlabel='Phred score', ylabel='Frequency', title='Phred score distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/phred_score_distribution.png")
        plt.close()

    @staticmethod
    def plot_length_distribution(df1, out):
        """
        Frequency of sizes. Bins auto-sized based on length distribution. Log scale x-axis.
        :param d: Dictionary
        :param out: path to output png file
        :return:
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        df = df1.loc[:, ('Length', 'Flag')]
        df_pass = df.loc[df['Flag'] == 'pass']
        df_fail = df.loc[df['Flag'] == 'fail']

        # Set bin sized for histogram
        min_len = min(df['Length'])
        max_len = max(df['Length'])
        min_exp = np.log10(min_len)
        max_exp = np.log10(max_len)
        len_logbins = np.logspace(min_exp, max_exp, 50)

        average_len_pass = int(df_pass['Length'].mean())
        ax.hist(df_pass['Length'], histtype='stepfilled', color='blue', alpha=0.6,
                label="pass (avg: {} bp)".format(average_len_pass), bins=len_logbins)

        if not df_fail.empty:
            average_len_fail = int(df_fail['Length'].mean())
            ax.hist(df_fail['Length'], histtype='stepfilled', color='red', alpha=0.6,
                    label="fail (avg: {} bp)".format(average_len_fail), bins=len_logbins)

        plt.legend()
        plt.xscale('log')
        ax.set(xlabel='Read length (bp)', ylabel='Frequency', title='Read length distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(out + "/length_distribution.png")
        plt.close()

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
        plt.close()

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
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='fail')
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
        plt.close()

    @staticmethod
    def plot_quality_vs_length_hex(df1, out):
        sns.set(style="ticks")

        df = df1.loc[:, ('Flag', 'Length', 'Qual')]
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
        # len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set bin sized for histogram
        len_logbins = np.logspace(min_exp, max_exp, 25)

        # Set y-axis limits
        min_phred = min(df['Qual'])
        max_phred = max(df['Qual'])

        # Set bin sized for histogram
        phred_bins = np.linspace(min_phred, max_phred, 15)

        # Do the Kernel Density Estimation (KDE)
        x = df_pass['Length']
        y = df_pass['Qual']

        # Create grid object
        g = sns.JointGrid(x='Length', y='Qual', data=df_pass, space=0)

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
            g.x = df_fail['Length']
            g.y = df_fail['Qual']

            # Do the Kernel Density Estimation (KDE)
            x = df_fail['Length']
            y = df_fail['Qual']

            g.ax_joint.hexbin(x, y, gridsize=50, cmap="Reds", xscale='log', alpha=0.6, mincnt=1, edgecolor='none')
            g.ax_marg_x.hist(x, histtype='stepfilled', color='red', alpha=0.6, bins=len_logbins)
            g.ax_marg_y.hist(y, histtype='stepfilled', color='red', alpha=0.6, bins=phred_bins,
                             orientation="horizontal")

        # Add legend to the joint plot area
        # https://matplotlib.org/tutorials/intermediate/legend_guide.html
        blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='pass')
        red_patch = mpatches.Patch(color='red', alpha=0.6, label='fail')
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
        plt.close()

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
        plt.close()

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

        fig.savefig(out + "/test.png")
        plt.close()

    @staticmethod
    def plot_reads_vs_bp_per_sample(df1, out):
        # Only keep pass reads
        df = df1.loc[(df1['Flag'] == 'pass'), ['Name', 'Flag', 'Length']]

        # Drop the Flag column
        df = df.loc[:, ['Name', 'Length']].reset_index(drop=True)

        # Count number of reads and total bp per sample and add columns to dataframe
        df = df.groupby('Name', as_index=False).agg(Count=('Length', 'size'), TotalLength=('Length', 'sum'))
        # print('\n{}'.format(df))  # Debug

        # Make plot
        fig, ax1 = plt.subplots(figsize=(10, 6))  # In inches

        ind = np.arange(len(df['Name']))
        width = 0.35
        p1 = ax1.bar(ind, df['TotalLength'], width, color='#4B9BFF', bottom=0, edgecolor='black')

        ax2 = ax1.twinx()
        p2 = ax2.bar(ind+width, df['Count'], width, color='#FFB46E', bottom=0, edgecolor='black')

        ax1.set_title('Total Size Versus Total Reads Per Sample')
        ax1.set_xticks(ind + width / 2)
        ax1.set_xticklabels(tuple(df['Name']), rotation=45, ha='right')

        ax1.grid(False)
        ax2.grid(False)

        ax1.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        ax2.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

        ax1.legend((p1[0], p2[0]), ('TotalLength', 'Count'), bbox_to_anchor=(1.1, 1), loc=2)
        ax1.yaxis.set_units('Total bp')
        ax2.yaxis.set_units('Total reads')
        ax1.autoscale_view()

        plt.tight_layout()
        fig.savefig(out + "/reads_vs_bp_per_sample.png")
        plt.close()

    @staticmethod
    def plot_pores_output_vs_time_all(df1, out):
        import matplotlib.lines as mlines
        df = df1.loc[:, ('Flag', 'Time')]

        # convert datatime to elapsed minutes
        time_zero = df.loc[:, 'Time'].min()  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 60
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 60  # convert to hours
        df.loc[:, 'Time'] = df.loc[:, 'Time'].astype(int)  # convert to integer

        df_pass = df[df.loc[:, 'Flag'] == 'pass']
        df_fail = df[df.loc[:, 'Flag'] == 'fail']

        # Check how many 15-minute bins are required to plot all the data
        max_time = max(df.loc[:, 'Time'])
        min_time = min(df.loc[:, 'Time'])
        nbins = int(max_time / 15) if max_time % 15 == 0 else int(max_time / 15) + 1
        # Create the bin boundaries
        x_bins = np.linspace(min_time, max_time, nbins)  # every 15 min

        # Make plot
        fig, ax = plt.subplots()

        # If not fail, just draw the pass. Else, draw total, fail and pass
        if not df_fail.empty:  # fail reads might be missing if plotting filtered reads for example.
            # Pass
            hist, edges = np.histogram(df_pass['Time'], bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, alpha=0.5, linewidth=0, color='blue')  # cmap='Blues')

            # Fail
            hist, edges = np.histogram(df_fail['Time'], bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, alpha=0.5, linewidth=0, color='red')  # cmap='Reds')

        # Generate counts for each bin
        # hist, edges = np.histogram(time_list_all, bins=x_bins, density=False)
        hist, edges = np.histogram(df['Time'], bins=x_bins, density=False)
        # Plot the data
        sns.scatterplot(data=hist, x_bins=edges, alpha=0.5, linewidth=0, color='green')  # cmap='Greens')

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

        all_marker = mlines.Line2D([], [], color='green', alpha=0.6, label='all', marker='o',
                                   markersize=5, linestyle='None')
        pass_marker = mlines.Line2D([], [], color='blue', alpha=0.6, label='pass', marker='o',
                                    markersize=5, linestyle='None')
        fail_marker = mlines.Line2D([], [], color='red', alpha=0.6, label='fail', marker='o',
                                    markersize=5, linestyle='None')

        # if time_list_fail:
        if not df_fail.empty:
            ax.legend(handles=[all_marker, pass_marker, fail_marker], loc='upper right')
        else:
            green_circle = mlines.Line2D([], [], color='green', alpha=0.6, label='pass', marker='o',
                                         markersize=5, linestyle='None')
            ax.legend(handles=[all_marker], loc='upper right')

        # Add label to axes
        plt.title('Pores output over time')
        plt.ylabel('Reads per 15 minutes')
        plt.xlabel('Sequencing time (hours)')

        plt.tight_layout()  # Get rid of extra margins around the plot
        # fig = g.get_figure()  # Get figure from FacetGrid
        fig.savefig(out + "/pores_output_vs_time_all.png")
        plt.close()

    @staticmethod
    def plot_channel_output_all(df1, out):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749
        """

        # Sum pass and fail reads per channel
        df = df1.loc[:, ('Channel', 'Flag')]
        df = df.groupby(['Channel', 'Flag'], as_index=False).size().pivot('Channel', 'Flag', 'size')

        df_all = pd.DataFrame()
        df_all['All'] = df['pass'] + df['fail']
        df_pass = df[['pass']]  # The double square brackets keep the column name
        df_fail = df[['fail']]

        # Plot
        if not df_fail['fail'].sum() == 0:
            fig, axs = plt.subplots(nrows=3, figsize=(6, 12))

            for i, my_tuple in enumerate([(df_all, 'All', 'Greens'),
                                          (df_pass, 'pass', 'Blues'),
                                          (df_fail, 'fail', 'Reds')]):
                my_df = my_tuple[0]
                flag = my_tuple[1]
                cmap = my_tuple[2]

                maxval = max(my_df.index.astype(int))  # maximum channel value
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
            value_cts = pd.Series(df_pass['pass'])
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
        plt.close()

    @staticmethod
    def plot_pores_length_output_vs_time_all(df1, out):
        # Fetch and prepare data from dictionary
        df = df1.loc[:, ('Time', 'Length', 'Flag')]

        # convert datatime to elapsed hours
        time_zero = min(df.loc[:, 'Time'])  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours

        # Compute x_bins
        nbins = int(max(df['Time'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['Time']), max(df['Time']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['Flag'] == 'pass']
        sns.regplot(x=pass_df['Time'], y=pass_df['Length'], x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['Flag'] == 'fail']
        if not fail_df.empty:
            sns.regplot(x=fail_df['Time'], y=fail_df['Length'], x_bins=x_bins,
                        fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                        label='fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('Read length over time')
        plt.ylabel('Average reads length per 15 min.')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_length_output_vs_time_all.png")
        plt.close()

    @staticmethod
    def plot_pores_qual_output_vs_time_all(df1, out):
        # Fetch and prepare data from dictionary
        df = df1.loc[:, ('Time', 'Qual', 'Flag')]

        # convert datatime to elapsed hours
        time_zero = min(df.loc[:, 'Time'])  # looking for min of 1st elements of list of tuples
        df.loc[:, 'Time'] = df.loc[:, 'Time'] - time_zero
        try:
            df.loc[:, 'Time'] = df.loc[:, 'Time'].dt.total_seconds() / 3600  # convert to hours
        except AttributeError:
            df.loc[:, 'Time'] = df.loc[:, 'Time'] / 3600  # convert to hours

        # Compute x_bins
        nbins = int(max(df['Time'])) + 1  # How many bins to plot data
        x_bins = np.linspace(min(df['Time']), max(df['Time']), nbins)  # Create the bin boundaries

        # Plot
        fig, ax = plt.subplots()
        # Pass
        pass_df = df[df['Flag'] == 'pass']
        sns.regplot(data=pass_df, x='Time', y='Qual', x_bins=x_bins,
                    fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                    label='pass', color='blue')  # capsize=4, capthick=1
        # Fail
        fail_df = df[df['Flag'] == 'fail']
        if not fail_df.empty:
            sns.regplot(data=fail_df, x='Time', y='Qual', x_bins=x_bins,
                        fit_reg=False, scatter_kws={'alpha': 0.6, 's': 30},
                        label='fail', color='red')

        # Set major ticks every 4 h
        ax.xaxis.set_major_locator(MultipleLocator(4))  # Want every 4 hours

        # Add label to axes
        plt.title('Read quality over time')
        plt.ylabel('Average reads quality per 15 min.')
        plt.xlabel('Sequencing time (hours)')
        plt.legend()

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig.savefig(out + "/pores_qual_output_vs_time_all.png")
        plt.close()

    @staticmethod
    def plot_size_distribution_per_sample(df1, out):
        # Only keep pass reads
        df = df1.loc[:, ('Name', 'Flag', 'Length')]

        # Make plot
        fig, ax = plt.subplots(figsize=(10, 6))

        g = sns.boxplot(data=df, x='Name', y='Length', hue='Flag', palette=['blue', 'red'], showfliers=False)
        g.figure.suptitle('Length distribution per sample')
        # Remove legend title:
        g.legend_.set_title(None)

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/length_distribution_per_sample.png")
        plt.close()

    @staticmethod
    def plot_bp_reads_stacked_histo(df1, out):
        df = df1.loc[:, ('Name', 'Flag', 'Length')]
        df.reset_index(inplace=True, drop=True)  # Drop the index column from the dataframe

        # Group data by Sample and Time, count how many reads for each unit of time
        df = df.groupby(['Name', 'Flag'], as_index=False).agg(TotalLength=('Length', 'sum'), Count=('Length', 'size'))
        # df_total = df.groupby(['Name', 'Flag'], as_index=False).agg(TotalLength=('Length', 'sum'), Count=('Length', 'size'))
        # df_fail = df_total[df_total['Flag'] == 'fail']
        # Add cumulative sum for each time point
        df['LengthCumSum'] = df.groupby(['Name', 'Flag'], as_index=False).sum().groupby('Name')['TotalLength'].cumsum()

        # Make plot
        fig, ax = plt.subplots(figsize=(10, 6))

        g = sns.barplot(data=df, x='Name', y='Count', hue='Flag', ci=None)
        # g = sns.barplot(data=df, x='Name', y='LengthCumSum', hue='Flag', ci=None)

        # # Top
        # g = sns.barplot(data=df_total, x='Name', y='Count', color='blue', ci=None)
        #
        # # Bottom
        # g = sns.barplot(data=df_fail, x='Name', y='Count', color='red', estimator=sum, ci=None)

        g.figure.suptitle('Total read and length (bp) per sample')
        # Remove legend title
        g.legend_.set_title(None)

        # Print samples at 45 degrees
        # ind = np.arange(len(df['Name']))
        # width = 0.35
        # ax.set_xticks(ind + width / 2)
        # ax.set_xticklabels(tuple(df['Name']), rotation=45, ha='right')

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # accounts for the "suptitile" [left, bottom, right, top]
        fig.savefig(out + "/read_vs_length_per_sample.png")
        plt.close()


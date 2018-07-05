#!/usr/local/env python3

from nested_dict import nested_dict
import os.path
import glob
import numpy as np
from time import time


__author__ = 'duceppemo'


class NanoQC(object):

    def __init__(self, args):

        """Define objects based on supplied arguments"""
        self.args = args
        self.input_folder = args.input
        self.output_folder = args.output

        # Shared data structure(s)
        # self.sample_dict = collections.defaultdict(lambda: dict)
        self.sample_dict = nested_dict()

        # Create a list of fastq files in input folder
        self.input_fastq_list = list()
        for fastq_file in glob.iglob(self.input_folder + '/**/*fastq*', recursive=True):
            if os.path.isfile(fastq_file):
                 self.input_fastq_list.append(fastq_file)

        # run the script
        self.run()

    def run(self):
        """Run everything"""
        self.parse_fastq(self.input_fastq_list, self.sample_dict)
        print("\nPlotting stats...")
        self.graph_read_number(self.sample_dict)  # Reads over time per sample. One graph per sample

    def hbytes(self, num):
        for x in ['bytes', 'KB', 'MB', 'GB']:
            if num < 1024.0:
                return "%3.1f%s" % (num, x)
            num /= 1024.0
        return "%3.1f%s" % (num, 'TB')

    def elapsed_time(self, seconds):
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string

    def parse_fastq(self, f, d):
        """
        Parse a fastq file into a dictionary
        :param f: A fastq file handle
        :param d: A dictionary to store the relevant information about each sequence
        :return: An update dictionary
        https://www.biostars.org/p/317524/
        http://www.blopig.com/blog/2016/08/processing-large-files-using-python/
        """
        from dateutil.parser import parse
        import gzip

        for f in self.input_fastq_list:
            with gzip.open(f, 'rt') if f.endswith('gz') else open(f, 'rU') as file:
                start_time = time()

                name = os.path.basename(file.name).split('.')[0].split('_')[0]  # Everything before the 1st underscore
                # name = file.name.split('/')[-2]  # Parent folder name

                # get some stats about the input file
                statinfo = os.stat(file.name)
                file_size = self.hbytes(statinfo.st_size)

                flag = 'pass'  # Default value
                if 'fail' in file.name:
                    flag = 'fail'  # Check in path for the word "fail"
                    print("Parsing sample \"%s\" from \"%s\" folder (%s)... "
                          % (name, flag, file_size), end="", flush=True)
                elif 'pass' in file.name:
                    print("Parsing sample \"%s\" from \"%s\" folder (%s)..."
                          % (name, flag, file_size), end="", flush=True)
                else:
                    print("Parsing sample \"%s\" from \"%s\" folder (%s). Assuming \"pass\" reads..."
                          % (name, flag, file_size), end="", flush=True)

                lines = list()
                for line in file:
                    lines.append(line.rstrip())
                    if len(lines) == 4:
                        header, seq, extra, qual = lines  # get each component in a variable
                        # Sequence ID
                        seqid = header.split()[0][1:]
                        # Read Time stamp
                        try:
                            time_string = header.split()[4].split('=')[1]
                        except IndexError:
                            print(seqid)
                        #Sequence length
                        length = len(seq)
                        # Average phred score
                        phred_list = list()
                        for letter in qual:
                            phred_list.append(ord(letter))
                        average_phred = int(np.round(np.mean(phred_list)) - 33)
                        # GC percentage
                        gc = int(round((sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / len(seq)) * 100))

                        # Add to dictionary
                        d[name][seqid]['datetime'] = parse(time_string)  # 2018-05-04T03:14:13Z
                        d[name][seqid]['length'] = length
                        d[name][seqid]['quality'] = average_phred
                        d[name][seqid]['flag'] = flag
                        d[name][seqid]['gc'] = gc

                        lines = []  # empty list

                end_time = time()
                interval = end_time - start_time
                print("took %s" % self.elapsed_time(interval))

    def graph_read_number(self, d):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :return: A png file with the graph
        TODO -> use numpy to handle the plot data, on row per sample?
        TODO -> what if just fail reads? Adapt code!
        """

        import matplotlib.pyplot as plt

        #
        # Total yield, all samples combined
        #

        fig, ax = plt.subplots()
        t_pass = list()  # time
        t_fail = list()
        for name in d:
            for seqid in d[name]:
                t = d[name][seqid]['datetime']
                if d[name][seqid]['flag'] == "pass":
                    t_pass.append(t)
                else:
                    t_fail.append(t)

        # Convert time object in hours from beginning of run
        t_zero_pass = min(t_pass)

        # Find the smallest datetime value
        if t_fail:
            t_zero_fail = min(t_fail)
            t_zero = min(t_zero_pass, t_zero_fail)
        else:
            t_zero = t_zero_pass

        # Prepare datetime value for plotting
        t_pass[:] = [x - t_zero for x in t_pass]  # Subtract t_zero for the all time points
        t_pass.sort()  # Sort
        t_pass[:] = [x.days * 24 + x.seconds / 3600 for x in t_pass]  # Convert to hours (float)
        y_pass = range(1, len(t_pass) + 1, 1)  # Create range. 1 time point equals 1 read

        y_fail = list()
        if t_fail:
            t_fail[:] = [x - t_zero for x in t_fail]
            t_fail.sort()
            t_fail[:] = [x.days * 24 + x.seconds / 3600 for x in t_fail]
            y_fail = range(1, len(t_fail) + 1, 1)

        # Create plot
        ax.plot(t_pass, y_pass, color='blue')
        if t_fail:
            ax.plot(t_fail, y_fail, color='red')
            ax.legend(['Pass', 'Fail'])
        else:
            ax.legend(['Pass'])
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Total read yield')
        fig.savefig(self.output_folder + "/total_reads_vs_time.png")

        #
        # Yield per sample (stack plot), just the pass reads
        #

        fig, ax = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        legend_names = list()
        for name in d:
            legend_names.append(name)
            ts_pass = list()
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(ts)
            ts_zero = min(ts_pass)
            ts_pass[:] = [x - ts_zero for x in ts_pass]  # Subtract t_zero for the all time points
            ts_pass.sort()  # Sort
            ts_pass[:] = [x.days * 24 + x.seconds / 3600 for x in ts_pass]  # Convert to hours (float)
            ys_pass = range(1, len(ts_pass) + 1, 1)  # Create range. 1 time point equals 1 read
            ax.plot(ts_pass, ys_pass)
            ax.legend(legend_names)
        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Read yield per sample')
        ax.ticklabel_format(style='sci')
        fig.savefig(self.output_folder + "/reads_per_sample_vs_time.png")

        #
        # Read length per sample vs time
        #

        fig, ax = plt.subplots()
        legend_names = list()
        for name in d:
            legend_names.append(name)
            ts_pass = list()
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                length = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(tuple((ts, length)))

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
            ax.plot(x_values, y_values)

        ax.legend(legend_names)
        # Add axes labels
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Yield per sample in base pair\n("pass" only)')
        # ax.ticklabel_format(useOffset=False)  # Disable the offset on the x-axis
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis

        # Save figure to file
        fig.savefig(self.output_folder + "/bp_per_sample_vs_time.png")

        #
        # Sequence length vs time
        #

        fig, ax = plt.subplots()
        ts_pass = list()
        ts_fail = list()
        for name in d:
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                l = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(tuple((ts, l)))
                else:
                    ts_fail.append(tuple((ts, l)))

        ts_zero = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

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
        ax.plot(x_pass_values, y_pass_values, color='blue')
        ax.plot(x_fail_values, y_fail_values, color='red')
        ax.legend(['Pass', 'Fail'])
        ax.set(xlabel='Time (h)', ylabel='Number of base pairs', title='Total yield in base pair')
        fig.savefig(self.output_folder + "/total_bp_vs_time.png")

        #
        # Quality vs time
        #

        fig, ax = plt.subplots()
        ts_pass = list()
        ts_fail = list()
        for name in d:
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                l = d[name][seqid]['quality']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(tuple((ts, l)))
                else:
                    ts_fail.append(tuple((ts, l)))

        ts_zero = min(ts_pass, key=lambda x: x[0])[0]  # looking for min of 1st elements of the tuple list

        ts_pass1 = [tuple(((x - ts_zero), y)) for x, y in ts_pass]  # Subtract t_zero for the all time points
        ts_pass1.sort(key=lambda x: x[0])  # Sort according to first element in tuple
        ts_pass2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_pass1]  # Convert to hours (float)
        x_pass_values = [x for x, y in ts_pass2]
        y_pass_values = [y for x, y in ts_pass2]

        ts_fail1 = [tuple(((x - ts_zero), y)) for x, y in ts_fail]
        ts_fail1.sort(key=lambda x: x[0])
        ts_fail2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_fail1]
        x_fail_values = [x for x, y in ts_fail2]
        y_fail_values = [y for x, y in ts_fail2]

        # Print plot
        ax.plot(x_pass_values, y_pass_values, color='blue')
        ax.plot(x_fail_values, y_fail_values, color='red', alpha=0.5)
        ax.legend(['Pass', 'Fail'])
        ax.set(xlabel='Time (h)', ylabel='Phred score', title='Read quality evolution')
        fig.savefig(self.output_folder + "/quality_vs_time.png")

        #
        # Frequency of phred scores
        #

        fig, ax = plt.subplots()
        qual_pass = list()
        qual_fail = list()
        for name in d:
            for seqid in d[name]:
                qual = d[name][seqid]['quality']
                if d[name][seqid]['flag'] == "pass":
                    qual_pass.append(qual)
                else:
                    qual_fail.append(qual)
        mean_pass_qual = np.round(np.mean(qual_pass), 1)
        mean_fail_qual = np.round(np.mean(qual_fail), 1)

        # Print plot
        ax.hist([qual_pass, qual_fail],
                bins=np.arange(min(min(qual_pass), min(qual_fail)) - 1, max(max(qual_pass), max(qual_fail)) + 1),
                color=['blue', 'red'],
                label=["pass (Avg: %s)" % mean_pass_qual, "fail (Avg: %s)" % mean_fail_qual])
        plt.legend()
        ax.set(xlabel='Phred score', ylabel='Frequency', title='Phred score distribution')
        fig.savefig(self.output_folder + "/phred_score distribution.png")

        #
        # Frequency of sizes (bins of 1kb). Log scale
        #

        size_pass = list()
        size_fail = list()
        for name in d:
            for seqid in d[name]:
                size = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    size_pass.append(size)
                else:
                    size_fail.append(size)
        min_len = min(min(size_pass), min(size_fail))
        max_len = max(max(size_pass), max(size_fail))
        mean_pass_size = np.round(np.mean(size_pass), 1)
        mean_fail_size = np.round(np.mean(size_fail), 1)

        # Print plot
        fig, ax = plt.subplots()
        binwidth = int(np.round(10 * np.log10(max_len - min_len)))
        logbins = np.logspace(np.log10(min_len), np.log10(max_len), binwidth)
        plt.hist([size_pass, size_fail],
                bins=logbins,
                color=['blue', 'red'],
                 label=["pass (Avg: %s)" % mean_pass_size, "fail (Avg: %s)" % mean_fail_size])
        plt.legend()
        plt.xscale('log')
        ax.set(xlabel='Read length (bp)', ylabel='Frequency', title='Read length distribution')
        fig.savefig(self.output_folder + "/length_distribution.png")

        #
        # Heat map (length vs quality)
        #

        # fig, ax = plt.subplots()


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Check status of submissions on Phaster\'s server'
                                        'and download results when available')
    parser.add_argument('-i', '--input', metavar='myfile.fastq.gz',
                        required=True,
                        help='Input folder with fastq file(s),gzipped or not')
    parser.add_argument('-o', '--output', metavar='/qc/',
                        required=True,
                        help='Output folder')

    # Get the arguments into an object
    arguments = parser.parse_args()

    NanoQC(arguments)

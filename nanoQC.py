#!/usr/local/env python3

import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import multiprocessing as mp
from dateutil.parser import parse
import logging
import collections


__author__ = 'duceppemo'
__version__ = '0.2.1'


# TODO -> parse the fastq file if parallel, if gives a performance gain
#         Either share dictionary across threads, or create one dictionary per thread and merge at the end???
# TODO -> Add a function to rename the samples with a conversion table (two-column tab-separated file)
# TODO -> make proper log file
# TODO -> check for dependencies (e.g. pigz)
# TODO -> unit testing
# TODO -> benchmark using gzip library to decompress versus pigz system call
# TODO -> add option to use the "sequencing_summary.txt" file as input instead of the fastq files


class NanoQC(object):

    def __init__(self, args):

        """Define objects based on supplied arguments"""
        self.args = args
        self.input_folder = args.fastq
        self.input_summary = args.summary
        self.output_folder = args.output

        # Shared data structure(s)
        self.sample_dict = dict()
        self.summary_dict = dict()
        # self.sample_dict = nested_dict()
        # self.summary_dict = nested_dict()

        # Create a list of fastq files in input folder
        self.input_fastq_list = list()

        # Threading
        self.cpu = mp.cpu_count()

        # run the script
        self.run()

    def run(self):
        """
        Run everything
        :return:
        """

        self.check_dependencies()

        # Check if correct argument combination is being used
        self.check_args()

        # Select appropriate parser based on input type
        if self.input_folder:
            self.parse_fastq(self.input_fastq_list, self.sample_dict)

            # Check if there is data
            if not self.sample_dict:
                raise Exception('No data!')
            else:
                self.make_fastq_plots()  # make the plots for fastq files
        else:  # elif self.input_summary:
            self.parse_summary(self.summary_dict)

            # Check if there is data
            if not self.summary_dict:
                raise Exception('No data!')

        # import pprint
        # pp = pprint.PrettyPrinter(indent=4)
        # pp.pprint(self.sample_dict)

    def check_dependencies(self):
        pass

    def check_args(self):
        """
        Check if correct argument combination is being used
        Check if output folder exists
        :return:
        """
        import pathlib

        if self.input_folder and self.summary_dict:
            print('Please use only one input type ("-f" or "-s")')
            parser.print_help(sys.stderr)
            sys.exit(1)
        elif not self.input_folder and not self.input_summary:
            print('Please use one of the following input types ("-f" or "-s")')
            parser.print_help(sys.stderr)
            sys.exit(1)

        if not self.output_folder:
            print('Please specify an output folder ("-o")')
            parser.print_help(sys.stderr)
            sys.exit(1)
        else:
            pathlib.Path(self.output_folder).mkdir(parents=True, exist_ok=True)  # Create if if it does not exist

        if self.input_folder:
            for root, directories, filenames in os.walk(self.input_folder):
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(('.fastq','.fastq.gz')):
                        self.input_fastq_list.append(absolute_path)

            # check if input_fastq_list is not empty
            if not self.input_fastq_list:
                raise Exception("No fastq file found in %s!" % self.input_folder)

    def hbytes(self, num):
        """
        Convert bytes to KB, MB, GB or TB
        :param num: a file size in bytes
        :return: A string representing the size of the file with proper units
        """
        for x in ['bytes', 'KB', 'MB', 'GB']:
            if num < 1024.0:
                return "%3.1f%s" % (num, x)
            num /= 1024.0
        return "%3.1f%s" % (num, 'TB')

    def elapsed_time(self, seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formated time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string

    def optimal_line_number(self, f, cpu):
        # Determine optimal number of lines in the fastq file to put in a chunk
        total_lines = 0
        with open(f, 'rb') as file:
            while True:
                buffer = file.read(1024 * 1024 * 8)  # 8 MB buffer
                if not buffer:
                    break
                total_lines += buffer.count(b'\n')

        # make sure its a multiple of 4
        optimal = int(round(total_lines / cpu, 0))
        optimal += optimal % 4

        return optimal

    # Iterator that yields start and end locations of a file chunk of default size 1MB.
    def chunkify(self, f, size=1024 * 1024):
        file_end = os.path.getsize(f)
        with open(f, 'rb', 1024 * 1024) as file:
            chunk_end = file.tell()
            while True:
                chunk_start = chunk_end
                file.seek(size, 1)
                self._EOC(file)
                chunk_end = file.tell()
                yield chunk_start, chunk_end - chunk_start
                if chunk_end >= file_end:
                    break

    # read chunk
    def read(self, f, chunk):
        with open(f, 'rb') as file:
            file.seek(chunk[0])
            return file.read(chunk[1])

    # End of chunk
    def _EOC(self, file):
        l = file.readline()  # incomplete line
        p = file.tell()
        l = file.readline()
        while l and b'runid' not in l:  # find the start of sequence
            p = file.tell()
            l = file.readline()
        file.seek(p)  # revert one line

    def parse_fastq_to_dict(self, l, my_dict, name, flag):
        header, seq, extra, qual = l  # get each component od list in a variable

        # Sequence ID
        seq_id = header.split()[0][1:]

        # Read Time stamp
        time_string = header.split()[4].split(b'=')[1]
        time_string = parse(time_string)

        # Sequence length
        length = len(seq)

        # Average phred score
        phred_list = list()
        for letter in qual:
            phred_list.append(letter)
        average_phred = sum(phred_list) / len(phred_list) - 33
        average_phred = round(average_phred, 1)

        # GC percentage
        g_count = float(seq.count(b'G'))
        c_count = float(seq.count(b'C'))
        gc = int(round((g_count + c_count) / float(length) * 100))

        # Add to dictionary
        if name not in my_dict:  # empty dictionary
            my_dict[name] = {}
        if seq_id not in my_dict[name]:
            my_dict[name][seq_id] = {}

        my_dict[name][seq_id]['length'] = length
        my_dict[name][seq_id]['quality'] = average_phred
        my_dict[name][seq_id]['flag'] = flag
        my_dict[name][seq_id]['gc'] = gc
        my_dict[name][seq_id]['datetime'] = time_string  # 2018-05-04T03:14:13Z

    def get_chunk_data(self, f, name, flag, chunk):
        my_dict = {}

        while True:
            data = self.read(f, chunk)
            lines = []
            my_data = data.split(b'\n')

            for line in my_data:
                if not line:  # end of file?
                    break
                line = line.rstrip()

                if lines and b'runid' in line:  # new entry
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            return my_dict

    def update_nested_dict(self, original_dict, new_dict):
        """
        Update the dictionary and its nested dictionary fields.
        my_dict[name][seq_id]['length'] = length
        """

        for name, dict1 in new_dict.items():
            for seq_id, dict2 in dict1.items():
                if name not in original_dict.keys():
                    original_dict[name] = dict1  # Add sample and read
                else:
                # elif seq_id not in original_dict.values():  # it's not going to be. Every read has a unique seq_id
                    original_dict[name][seq_id] = dict2

        # for k, v in new_dict.items():
        #     if (k in original_dict.keys()
        #             and isinstance(original_dict[k], dict)
        #             and isinstance(new_dict[k], collections.Mapping)):
        #         self.update_nested_dict(original_dict[k], new_dict[k])
        #     else:
        #         original_dict[k] = new_dict[k]

        return original_dict

    def parse_fastq(self, l, d):
        """
        Parse a basecalled nanopore fastq file by Albacore into a dictionary.
        Multiple fastq files shall be from the same sequencing run.
        :param l: A list of fastq files.
        :param d: An empty dictionary
        :return: A dictionary with the time stamp, length, average quality, pass/fail flag and gc %
                 for each read of each sample
        https://www.biostars.org/p/317524/
        http://www.blopig.com/blog/2016/08/processing-large-files-using-python/
        """

        import subprocess

        for f in l:
            name = os.path.basename(f).split('.')[0].split('_')[0]  # basename before 1st "_" -> sample name
            # name = file.name.split('/')[-2]  # Parent folder name

            # get some stats about the input file
            stat_info = os.stat(f)
            file_size = self.hbytes(stat_info.st_size)

            flag = 'pass'  # Default value
            if 'fail' in f:
                flag = 'fail'  # Check in path for the word "fail"
                print("Parsing sample \"%s\" from \"%s\" folder (%s)... "
                      % (name, flag, file_size), end="", flush=True)
            elif 'pass' in f:
                print("Parsing sample \"%s\" from \"%s\" folder (%s)..."
                      % (name, flag, file_size), end="", flush=True)
            else:
                print("Parsing sample \"%s\" from \"%s\" folder (%s). Assuming \"pass\" reads..."
                      % (name, flag, file_size), end="", flush=True)

            filename = os.path.basename(f)
            filename_no_gz = '.'.join(filename.split('.')[:-1])
            gz_flag = 0

            start_time = time()

            if f.endswith('.fastq.gz'):
                p1 = subprocess.Popen(['pigz', '-d', '-k', '-f', '-c', f], stdout=subprocess.PIPE)
                (outs, errs) = p1.communicate()
                with open('/tmp/' + filename_no_gz, 'b+w') as tmp_file:
                    tmp_file.write(outs)
                f = "/tmp/" + filename_no_gz
                gz_flag = 1
            if f.endswith(".fastq"):
                pass
            else:
                raise Exception('Wrong file extension. Please use ".fastq" or ".fastq.gz"')

            # start_time = time()
            pool = mp.Pool(self.cpu)

            jobs = []
            for chunk in self.chunkify(f):
                job = pool.apply_async(self.get_chunk_data, args=(f, name, flag, chunk,))
                jobs.append(job)

            output = []
            for job in jobs:
                output.append(job.get())

            pool.close()
            pool.join()

            # Merge dictionary
            # Update self.sample_dict with results from every chunk
            tmp_dict = dict()
            reads = 0
            # print("Merging processes outputs")
            for dictionary in output:
                # Count number of reads in the sample
                for key, value in dictionary.items():
                    reads += len(value)
                # d = {**d, **dictionary}  # Do the merge
                # d.update(dictionary)  # Do the merge
                self.update_nested_dict(tmp_dict, dictionary)  # merge the output from the sample

            # print("Merging to main dictionary")
            self.update_nested_dict(d, tmp_dict)  # Merge sample to main dictionary

            end_time = time()
            interval = end_time - start_time
            print("took {} ({} reads)".format(self.elapsed_time(interval), reads))

            #  Remove tmp file
            if gz_flag == 1:
                os.remove('/tmp/' + filename_no_gz)

    def parse_summary(self, d):
        """
        Parse "sequencing_summary.txt" file from Albacore
        :param d: Empty summary dictionary
        :return: Dictionary with info about
        """

        start_time = time()
        with open(self.input_summary, 'rb', 8192*1024) as file:
            fields = list()
            read_counter = 0
            next(file)  # skip first line
            for line in file:
                line = line.rstrip()
                fields.append(line.split(b"\t"))
                if not line:
                    continue

                # Only keep  fields of interest
                name = fields[19]
                seqid = fields[1]
                channel = fields[3]
                events = fields[6]
                start_time = fields[4]
                flag = fields[7]

                # Populate the dictionary
                d[name][seqid]['channel'] = channel
                d[name][seqid]['events'] = events
                d[name][seqid]['start_time'] = start_time
                d[name][seqid]['flag'] = flag

        # Print read stats
        end_time = time()
        interval = end_time - start_time
        print("took %s for %d reads" % (self.elapsed_time(interval), read_counter))

    def make_fastq_plots(self):
        # print("\nMaking plots...", end="", flush=True)
        d = self.sample_dict

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_reads_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_reads_per_sample_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_bp_per_sample_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_bp_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        self.plot_phred_score_distribution(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_length...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_length(d)
        end_time = time()
        interval = end_time - start_time
        print("Took %s" % self.elapsed_time(interval))

        # self.plot_test(d)

    def plot_total_reads_vs_time(self, d):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :return: A png file with the graph
        TODO -> use numpy to handle the plot data, on row per sample?
        """

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

        # Find the smallest datetime value
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
        if t_pass:
            t_pass[:] = [x - t_zero for x in t_pass]  # Subtract t_zero for the all time points
            t_pass.sort()  # Sort
            t_pass[:] = [x.days * 24 + x.seconds / 3600 for x in t_pass]  # Convert to hours (float)
            y_pass = range(1, len(t_pass) + 1, 1)  # Create range. 1 time point equals 1 read

        if t_fail:
            y_fail = list()
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
        fig.savefig(self.output_folder + "/total_reads_vs_time.png")

    def plot_reads_per_sample_vs_time(self, d):
        """
        Plot yield per sample. Just the pass reads
        :param d: Dictionary
        :return: png file
        """

        # fig, ax = plt.subplots()
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches
        legend_names = list()
        for name in d:
            legend_names.append(name)
            ts_pass = list()
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(ts)
            if not ts_pass:
                return

            ts_zero = min(ts_pass)
            ts_pass[:] = [x - ts_zero for x in ts_pass]  # Subtract t_zero for the all time points
            ts_pass.sort()  # Sort
            ts_pass[:] = [x.days * 24 + x.seconds / 3600 for x in ts_pass]  # Convert to hours (float)
            ys_pass = range(1, len(ts_pass) + 1, 1)  # Create range. 1 time point equals 1 read

            # ax.plot(ts_pass, ys_pass)
            ax.plot(ts_pass, ys_pass,
                    label="%s (%s)" % (name, "{:,}".format(max(ys_pass))))
            # ax.legend(legend_names)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # New

        ax.set(xlabel='Time (h)', ylabel='Number of reads', title='Read yield per sample')
        ax.ticklabel_format(style='plain')  # Disable the scientific notation on the y-axis
        # comma-separated numbers to the y axis
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()  #
        fig.savefig(self.output_folder + "/reads_per_sample_vs_time.png")

    def plot_bp_per_sample_vs_time(self, d):
        """
        Read length per sample vs time
        :param d: Dictionary
        :return: png file
        """

        fig, ax = plt.subplots(figsize=(10, 6))  # In inches
        for name in d:
            ts_pass = list()
            for seqid in d[name]:
                ts = d[name][seqid]['datetime']
                length = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    ts_pass.append(tuple((ts, length)))

            if not ts_pass:
                return
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
        fig.savefig(self.output_folder + "/bp_per_sample_vs_time.png")

    def plot_total_bp_vs_time(self, d):
        """
        Sequence length vs time
        :param d: Dictionary
        :return: png file
        """

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
        fig.savefig(self.output_folder + "/total_bp_vs_time.png")

    def plot_quality_vs_time(self, d):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :return: png file
        """

        fig, ax = plt.subplots(figsize=(10, 4))
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

            df_pass = pd.DataFrame(list(ts_pass3), columns=['Time (h)', 'Phred Score'])  # Convert to dataframe
            df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        if ts_fail:
            ts_fail1 = [tuple(((x - ts_zero), y)) for x, y in ts_fail]
            ts_fail1.sort(key=lambda x: x[0])
            ts_fail2 = [tuple(((x.days * 24 + x.seconds / 3600), y)) for x, y in ts_fail1]
            ts_fail3 = [tuple((int(np.round(x)), y)) for x, y in ts_fail2]

            df_fail = pd.DataFrame(list(ts_fail3), columns=['Time (h)', 'Phred Score'])
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
            g = sns.violinplot(x='Time (h)', y='Phred Score', data=data, hue='Flag', split=True, inner=None)
            g.figure.suptitle('Phred score evolution')
        elif ts_pass3:
            g = sns.violinplot(x='Time (h)', y='Phred Score', data=data, inner=None)
            g.figure.suptitle('Phred score evolution\n(pass)')
        else:  # elif ts_fail3:
            g = sns.violinplot(x='Time (h)', y='Phred Score', data=data, inner=None)
            g.figure.suptitle('Phred score evolution\n(fail)')

        plt.tight_layout()
        fig.savefig(self.output_folder + "/quality_vs_time.png")

    def plot_phred_score_distribution(self, d):
        """
        Frequency of phred scores
        :param d: Dictionary
        :return: png file
        """

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

        if qual_pass:
            mean_pass_qual = np.round(np.mean(qual_pass), 1)

        if qual_fail:
            mean_fail_qual = np.round(np.mean(qual_fail), 1)

        # Print plot
        if qual_pass and qual_fail:
            ax.hist([qual_pass, qual_fail],
                    bins=np.arange(min(min(qual_pass), min(qual_fail)), max(max(qual_pass), max(qual_fail))),
                    color=['blue', 'red'],
                    label=["pass (Avg: %s)" % mean_pass_qual, "fail (Avg: %s)" % mean_fail_qual])
        elif qual_pass:
            ax.hist(qual_pass,
                    bins=np.arange(min(qual_pass), max(qual_pass)),
                    color='blue',
                    label="pass (Avg: %s)" % mean_pass_qual)
        else:
            ax.hist(qual_fail,
                    bins=np.arange(min(qual_fail), max(qual_fail)),
                    color='red',
                    label="fail (Avg: %s)" % mean_fail_qual)
        plt.legend()
        ax.set(xlabel='Phred score', ylabel='Frequency', title='Phred score distribution')
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(self.output_folder + "/phred_score_distribution.png")

    def plot_length_distribution(self, d):
        """
        Frequency of sizes. Bins auto-sized based on length distribution. Log scale x-axis.
        :param d: Dictionary
        :return: png file
        """

        size_pass = list()
        size_fail = list()
        for name in d:
            for seqid in d[name]:
                size = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    size_pass.append(size)
                else:
                    size_fail.append(size)
        if size_fail:  # assume "pass" always present
            min_len = min(min(size_pass), min(size_fail))
            max_len = max(max(size_pass), max(size_fail))
            mean_pass_size = np.round(np.mean(size_pass), 1)
            mean_fail_size = np.round(np.mean(size_fail), 1)
        elif size_pass:
            min_len = min(size_pass)
            max_len = max(size_pass)
            mean_pass_size = np.round(np.mean(size_pass), 1)
        else:
            print("No reads detected!")
            return

        # Print plot
        fig, ax = plt.subplots()
        binwidth = int(np.round(10 * np.log10(max_len - min_len)))
        logbins = np.logspace(np.log10(min_len), np.log10(max_len), binwidth)
        if size_fail:
            plt.hist([size_pass, size_fail],
                     bins=logbins,
                     color=['blue', 'red'],
                     label=["pass (Avg: %s)" % mean_pass_size, "fail (Avg: %s)" % mean_fail_size])
        else:
            plt.hist(size_pass, bins=logbins, color='blue', label="pass (Avg: %s)" % mean_pass_size)
        plt.legend()
        plt.xscale('log')
        ax.set(xlabel='Read length (bp)', ylabel='Frequency', title='Read length distribution')
        plt.tight_layout()
        fig.savefig(self.output_folder + "/length_distribution.png")

    def plot_quality_vs_length(self, d):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :return: png file
        """

        qs_pass = list()
        for name in d:
            for seqid in d[name]:
                qual = d[name][seqid]['quality']
                size = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    qs_pass.append(tuple((size, qual)))

        if not qs_pass:
            return

        df_pass = pd.DataFrame(list(qs_pass), columns=['Length (bp)', 'Phred Score'])

        min_len = pd.DataFrame.min(df_pass['Length (bp)'])
        min_exp = np.log10(min_len)
        min_value = float(10 ** (min_exp - 0.1))
        if min_value <= 0:
            min_value = 1
        max_len = pd.DataFrame.max(df_pass['Length (bp)'])
        max_exp = np.log10(max_len)
        max_value = float(10 ** (max_exp + 0.1))
        g = sns.jointplot(x='Length (bp)', y='Phred Score', data=df_pass, kind='kde',
                          stat_func=None,
                          xlim=[min_value, max_value],
                          space=0,
                          n_levels=50)
        ax = g.ax_joint
        ax.set_xscale('log')
        g.ax_marg_x.set_xscale('log')
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        g.savefig(self.output_folder + "/quality_vs_length.png")

    def plot_test(self, d):

        from scipy import stats

        qs_pass = list()
        for name in d:
            for seqid in d[name]:
                qual = d[name][seqid]['quality']
                size = d[name][seqid]['length']
                if d[name][seqid]['flag'] == "pass":
                    qs_pass.append(tuple((size, qual)))

        if not qs_pass:
            return

        df_pass = pd.DataFrame(list(qs_pass), columns=['Length (bp)', 'Phred Score'])

        min_len = pd.DataFrame.min(df_pass['Length (bp)'])
        min_exp = np.log10(min_len)
        min_value = float(10 ** (min_exp - 0.1))
        if min_value <= 0:
            min_value = 1
        max_len = pd.DataFrame.max(df_pass['Length (bp)'])
        max_exp = np.log10(max_len)
        max_value = float(10 ** (max_exp + 0.1))

        min_phred = df_pass['Phred Score'].min()
        max_phred = df_pass['Phred Score'].max()

        binwidth = int(np.round(10 * np.log10(max_len - min_len)))
        logbins = np.logspace(np.log10(min_len), np.log10(max_len), binwidth)

        # Initialize the figure
        sns.set_style("ticks")
        g = sns.JointGrid(x='Length (bp)', y='Phred Score', data=df_pass,
                          xlim=[min_value, max_value],
                          ylim=[min_phred, max_phred],
                          space=0)
        ax = g.ax_joint
        ax.set_xscale('log')
        g.ax_marg_x.set_xscale('log')

        g.plot_joint(sns.kdeplot, shade=True, n_levels=100)
        # g.plot_marginals(sns.kdeplot, shade=True)
        g.plot_marginals(sns.distplot, kde=True)

        # g.fig.set_figwidth(8)
        # g.fig.set_figheight(4)

        # Save
        g.savefig(self.output_folder + "/test.png")


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Plot QC data from nanopore sequencing run')
    parser.add_argument('-f', '--fastq', metavar='/basecalled/folder/',
                        required=False,
                        help='Input folder with fastq file(s),gzipped or not')
    parser.add_argument('-s', '--summary', metavar='sequencing_summary.txt',
                        required=False,
                        help='The "sequencing_summary.txt" file produced by the Albacore basecaller')
    parser.add_argument('-o', '--output', metavar='/qc/',
                        required=True,
                        help='Output folder')

    # Get the arguments into an object
    arguments = parser.parse_args()

    NanoQC(arguments)

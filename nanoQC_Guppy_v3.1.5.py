#!/usr/local/env python3.6

import os
import sys
from argparse import ArgumentParser
import numpy as np
from time import time
import multiprocessing as mp
from collections import defaultdict
import functions
from fastq_parser import FastqParser
from fastq_graphs import FastqPlots
from summary_parser import SummaryParser
from summary_graphs import SummaryPlots


__author__ = 'duceppemo'
__version__ = '0.3.6'


# TODO -> Check if can parallel parse (chunks) the sample processed in parallel?
# TODO -> Add a function to rename the samples with a conversion table (two-column tab-separated file)
# TODO -> make proper log file
# TODO -> Add stats (# of reads, # of pass, # of fail, top 10 channels, etc.)
# TODO -> check for dependencies
# TODO -> unit testing


class NanoQC(object):

    def __init__(self, args):

        """Define objects based on supplied arguments"""
        self.args = args
        self.input_folder = args.fastq
        self.input_summary = args.summary
        self.output_folder = args.output

        # Shared data structure(s)
        self.sample_dict = defaultdict()
        self.summary_dict = defaultdict()

        # Create a list of fastq files in input folder
        self.input_fastq_list = list()

        # Time tracking
        self.total_time = list()

        # Threading
        self.cpu = mp.cpu_count()

        # run the script
        self.run()

    def run(self):
        """
        Run everything
        :return:
        """

        self.total_time.append(time())

        self.check_dependencies()

        # Check if correct argument combination is being used
        self.check_args()

        # Select appropriate parser based on input type
        if self.input_folder:
            print("Parsing fastq files...", end="", flush=True)
            start_time = time()
            read_counter = FastqParser.parse_fastq_parallel(self.input_fastq_list, self.sample_dict, self.cpu)
            end_time = time()
            interval = end_time - start_time
            print(" took %s for %d reads" % (self.elapsed_time(interval), read_counter))

            # Check if there is data
            if not self.sample_dict:
                raise Exception('No data!')
            else:
                self.make_fastq_plots(self.sample_dict)  # make the plots for fastq files
        else:  # elif self.input_summary:
            print("Parsing summary file...", end='', flush=True)
            start_time = time()
            read_counter = SummaryParser.parse_summary(self.input_summary, self.summary_dict)
            # Print read stats
            end_time = time()
            interval = end_time - start_time
            print(" took %s for %d reads" % (self.elapsed_time(interval), read_counter))

            # Check if there is data
            if not self.summary_dict:
                raise Exception('No data!')
            else:
                self.make_summary_plots(self.summary_dict)

        # import pprint
        # pp = pprint.PrettyPrinter(indent=4)
        # pp.pprint(self.sample_dict)

        self.total_time.append(time())
        print("\n Total run time: {}".format(self.elapsed_time(self.total_time[1] - self.total_time[0])))

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
        :return: A string repreresenting the size of the file with proper units
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

    def make_fastq_plots(self, d):
        out = self.output_folder
        print("\nMaking plots:")

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_output_vs_time_all(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_pie...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_bp_per_sample_pie(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_phred_score_distribution(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_length_distribution(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_qual_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_qual_output_vs_time_all(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_length_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_length_output_vs_time_all(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_pie...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_per_sample_pie(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_channel_output_all(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_total_reads_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_total_bp_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_vs_bp_per_sample(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_per_sample_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_bp_per_sample_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting pores_output_vs_time_total...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_pores_output_vs_time_total(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_quality_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_scatter...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_quality_vs_length_scatter(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_quality_vs_length_hex(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_kde...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_quality_vs_length_kde(d, out)
        # # self.test_plot(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting gc_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_vs_time(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting gc_vs_length_hex...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_vs_length_hex(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_gc_output_vs_time_all(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting gc_vs_qual_vs_time_3D...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_vs_qual_vs_time_3D(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_gc_output_vs_time_per_sample(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

    def make_summary_plots(self, d):
        out = self.output_folder

        print("\nMaking plots:")

        print('\tPlotting pores_length_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_length_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_qual_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_qual_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_phred_score_distribution_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_length_distribution_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_pie...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_per_sample_pie_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_pie...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_bp_per_sample_pie_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_total_reads_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_total_bp_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_per_sample_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_bp_per_sample_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_quality_vs_length_hex_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_vs_bp_per_sample_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting pores_output_vs_time...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_pores_output_vs_time_summary(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting channel_output_total...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_channel_output_total(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))
        #
        # print('\tPlotting channel_output_pass_fail...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_channel_output_pass_fail(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_channel_output_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_quality_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))


if __name__ == '__main__':
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

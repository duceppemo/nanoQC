#!/usr/local/env python3
import multiprocessing as mp
import os
import sys
from argparse import ArgumentParser
import numpy as np
from time import time
from collections import defaultdict
from fastq_parser import FastqParser
from fastq_graphs import FastqPlots
from summary_parser import SummaryParser
from summary_graphs import SummaryPlots
import pathlib
import pandas as pd


__author__ = 'duceppemo'
__version__ = '0.2.0'


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
        self.cpu = args.threads
        self.parallel = args.parallel
        self.single = args.individual

        # Shared data structure(s)
        self.sample_dict = defaultdict()
        self.summary_dict = defaultdict()

        # run the script
        self.run()

    def run(self):
        # Runtime tracking
        beginning_time = time()

        self.check_dependencies()

        # Check if correct argument combination is being used
        self.check_args()

        # Select appropriate parser based on input type
        # Create a list of fastq files in input folder
        input_fastq_list = NanoQC.list_fastq_files(self.input_folder)
        if self.input_folder:
            print("Parsing fastq files to dictionary...", end="", flush=True)
            start_time = time()
            sample_dict = FastqParser.parallel_process_fastq(input_fastq_list, self.cpu, self.parallel)
            end_time = time()
            interval = end_time - start_time
            print(" took {} for {} reads".format(NanoQC.elapsed_time(interval), len(sample_dict)))

            # Check if there is data
            if not sample_dict:
                raise Exception('No data!')
            else:
                df = NanoQC.dict_2_df(sample_dict)
                # df.to_csv(self.output_folder + '/pd.tsv', sep='\t', index=False)  # Debug

                if self.single:
                    print('\tPlotting pores_gc_output_vs_time_per_sample...', end="", flush=True)
                    start_time = time()
                    FastqPlots.plot_pores_gc_output_vs_time_per_sample(df, self.output_folder)
                    end_time = time()
                    interval = end_time - start_time
                    print(" took %s" % NanoQC.elapsed_time(interval))

                    name_list = list(
                        map(lambda x: os.path.basename(x).split('.')[0].replace('_pass', '').replace('_fail', ''),
                            input_fastq_list))
                    for name in name_list:
                        single_sample_dict = dict()
                        for read_id, info in sample_dict.items():
                            if info.name == name:
                                single_sample_dict[read_id] = info
                        out_folder = self.output_folder + '/' + name
                        NanoQC.make_folder(out_folder)
                        print("Making plots for {}:".format(name))
                        dfs = NanoQC.dict_2_df(single_sample_dict)

                        NanoQC.make_fastq_plots(dfs, out_folder)
                else:
                    NanoQC.make_fastq_plots(df, self.output_folder)  # make the plots for fastq files
        else:  # elif self.input_summary:
            print("Parsing summary file...", end='', flush=True)
            start_time = time()
            read_counter = SummaryParser.parse_summary(self.input_summary, self.summary_dict)
            # Print read stats
            end_time = time()
            interval = end_time - start_time
            print(" took %s for %d reads" % (NanoQC.elapsed_time(interval), read_counter))

            # Check if there is data
            if not self.summary_dict:
                raise Exception('No data!')
            else:
                self.make_summary_plots(self.summary_dict)

        ending_time = time()
        print("\n Total run time: {}".format(NanoQC.elapsed_time(ending_time - beginning_time)))

    def check_dependencies(self):
        pass

    def check_args(self):
        """
        Check if correct argument combination is being used
        Check if output folder exists
        :return:
        """
        # Input folder
        if self.input_folder and self.summary_dict:
            print('Please use only one input type ("-f" or "-s")')
            parser.print_help(sys.stderr)
            sys.exit(1)

        elif not self.input_folder and not self.input_summary:
            print('Please use one of the following input types ("-f" or "-s")')
            parser.print_help(sys.stderr)
            sys.exit(1)

        # Output folder
        if not self.output_folder:
            print('Please specify an output folder ("-o")')
            parser.print_help(sys.stderr)
            sys.exit(1)
        else:
            NanoQC.make_folder(self.output_folder)

        # Simple
        if self.single:
            if self.input_summary:
                raise Exception('Per sample reporting has not been implemented yet for the summary file, sorry...')

        # Number of CPU
        if self.cpu > max_cpu:
            self.cpu = max_cpu

        # Number of samples to process in parallel
        if self.parallel > max_cpu:
            self.parallel = max_cpu

    @staticmethod
    def list_fastq_files(input_folder):
        input_fastq_list = list()

        # Walk input folder recursively to find fastq files
        if input_folder and os.path.isdir(input_folder):
            for root, directories, filenames in os.walk(input_folder):
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                        input_fastq_list.append(absolute_path)

        # Check if input_fastq_list is not empty
        if not input_fastq_list:
            raise Exception("No fastq file found in %s!" % input_folder)

        return input_fastq_list

    @staticmethod
    def hbytes(num):
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

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string

    @staticmethod
    def make_folder(folder):
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def dict_2_df(sample_dict):
        print('Converting dictionary to dataframe...', end="", flush=True)
        start_time = time()
        # Convert dictionary to pandas data frame
        my_dict = defaultdict()
        for seq_id, seq in sample_dict.items():
            if not seq.time_string:
                continue
            my_dict[seq_id] = (seq.name, seq.length, seq.flag, seq.average_phred, seq.gc, seq.time_string, seq.channel)
        df = pd.DataFrame.from_dict(my_dict, orient='index',
                                    columns=['Name', 'Length', 'Flag', 'Qual', 'GC', 'Time', 'Channel'])
        end_time = time()
        interval = end_time - start_time
        print(" took {}".format(NanoQC.elapsed_time(interval)))

        return df

    @staticmethod
    def make_fastq_plots(df, out):
        print("\nMaking plots:")

        print('\tPlotting read_vs_length_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_bp_reads_stacked_histo(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting plot_gc_distribution_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_distribution_per_sample(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting plot_size_distribution_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_size_distribution_per_sample(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting gc_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_quality_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_per_sample_pie...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_per_sample_pie(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_gc_output_vs_time_per_sample(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_gc_output_vs_time_all(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting gc_vs_length_hex...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_gc_vs_length_hex(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_quality_vs_length_hex(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_bp_per_sample_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_per_sample_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_reads_vs_bp_per_sample(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_total_bp_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_total_reads_vs_time(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_channel_output_all(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_output_vs_time_all(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting bp_per_sample_pie...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_bp_per_sample_pie(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_phred_score_distribution(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_length_distribution(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_qual_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_qual_output_vs_time_all(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_length_output_vs_time_all...', end="", flush=True)
        start_time = time()
        FastqPlots.plot_pores_length_output_vs_time_all(df, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_scatter...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_quality_vs_length_scatter(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_kde...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_quality_vs_length_kde(d, out)
        # # self.test_plot(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))



        # print('\tPlotting gc_vs_qual_vs_time_3D...', end="", flush=True)
        # start_time = time()
        # FastqPlots.plot_gc_vs_qual_vs_time_3D(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))

    def make_summary_plots(self, d):
        out = self.output_folder

        print("\nMaking plots:")

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_quality_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_length_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_length_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_qual_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_qual_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_phred_score_distribution_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_length_distribution_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_per_sample_pie...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_per_sample_pie_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting bp_per_sample_pie...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_bp_per_sample_pie_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_total_reads_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_total_bp_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_per_sample_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_bp_per_sample_vs_time_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_quality_vs_length_hex_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_reads_vs_bp_per_sample_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_pores_output_vs_time_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))

        # print('\tPlotting pores_output_vs_time...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_pores_output_vs_time_summary(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))

        # print('\tPlotting channel_output_total...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_channel_output_total(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))
        #
        # print('\tPlotting channel_output_pass_fail...', end="", flush=True)
        # start_time = time()
        # SummaryPlots.plot_channel_output_pass_fail(d, out)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % NanoQC.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        SummaryPlots.plot_channel_output_all_summary(d, out)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % NanoQC.elapsed_time(interval))


if __name__ == '__main__':
    max_cpu = mp.cpu_count()

    parser = ArgumentParser(description='Create QC plots using nanopore sequencing or basecalling data')
    parser.add_argument('-f', '--fastq', metavar='/basecalled/folder/',
                        required=False,
                        help='Input folder with fastq file(s), gzipped or not')
    parser.add_argument('-s', '--summary', metavar='sequencing_summary.txt',
                        required=False,
                        help='The "sequencing_summary.txt" file produced by guppy_basecaller')
    parser.add_argument('-o', '--output', metavar='/qc/',
                        required=True,
                        help='Output folder')
    parser.add_argument('-t', '--threads', metavar='{}'.format(max_cpu),
                        required=False, default=max_cpu,
                        type=int,
                        help='Number of CPU. Default {}'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='4',
                        required=False, default=4,
                        type=int,
                        help='Number of samples to process in parallel.')
    parser.add_argument('-i', '--individual',
                        action='store_true',
                        help='Produce a nanoQC report for each files instead of combining them in a single report')

    # Get the arguments into an object
    arguments = parser.parse_args()

    NanoQC(arguments)

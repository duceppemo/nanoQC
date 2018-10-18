#!/usr/local/env python3

import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec
import multiprocessing as mp
import threading
from queue import Queue
from multiprocessing.pool import ThreadPool
from dateutil.parser import parse
import logging
from collections import defaultdict, OrderedDict
from itertools import zip_longest
from functools import partial
from contextlib import closing
from math import ceil
import subprocess
import gzip
from matplotlib.ticker import FuncFormatter, MultipleLocator
import matplotlib.patches as mpatches
import datetime as dt
from math import sqrt
from itertools import count, islice


__author__ = 'duceppemo'
__version__ = '0.3.3'


# TODO -> Check if can parallel parse (chunks) the sample processed in parallel?
# TODO -> Add a function to rename the samples with a conversion table (two-column tab-separated file)
# TODO -> make proper log file
# TODO -> check for dependencies
# TODO -> unit testing
# TODO -> add option to use the "sequencing_summary.txt" file as input instead of the fastq files


class FastqObjects(object):
    def __init__(self, name, length, flag, average_phred, gc, time_string, channel):
        # Create seq object with its attributes
        self.name = name
        self.length = length
        self.flag = flag
        self.average_phred = average_phred
        self.gc = gc
        self.time_string = time_string
        self.channel = channel


class SummaryObjects(object):
    def __init__(self, name, length, channel, events, average_phred, time_stamp, flag):
        # Create seq object with its attributes
        self.name = name
        self.channel = channel
        self.events = events
        self.time_stamp = time_stamp
        self.length = length
        self.flag = flag
        self.average_phred = average_phred


class Layout(object):
    def __init__(self, structure, template, xticks, yticks):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks


class NanoQC(object):

    def __init__(self, args):

        """Define objects based on supplied arguments"""
        self.args = args
        self.input_folder = args.fastq
        self.input_summary = args.summary
        self.output_folder = args.output

        # Shared data structure(s)
        # self.sample_dict = dict()
        # self.summary_dict = dict()
        # self.sample_dict = nested_dict()
        # self.summary_dict = nested_dict()
        self.sample_dict = defaultdict()
        self.summary_dict = defaultdict()

        # Create a list of fastq files in input folder
        self.input_fastq_list = list()

        #Time tracking
        self.total_time = list()

        # Threading
        self.cpu = mp.cpu_count()
        # self.my_queue = Queue(maxsize=0)

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
            # self.parse_fastq(self.input_fastq_list, self.sample_dict)
            self.parse_fastq_parallel(self.input_fastq_list, self.sample_dict)

            # Check if there is data
            if not self.sample_dict:
                raise Exception('No data!')
            else:
                self.make_fastq_plots(self.sample_dict)  # make the plots for fastq files
        else:  # elif self.input_summary:
            self.parse_summary(self.summary_dict)

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
        with open(f, 'rb', 1024 * 1024) as file_handle:  # 1MB buffer
            while True:
                buffer = file_handle.read(1024 * 1024)  # 1MB buffer
                if not buffer:
                    break
                total_lines += buffer.count(b'\n')

        # make sure its a multiple of 4
        optimal = int(round(total_lines / cpu, 0))
        optimal += optimal % 4

        return optimal

    # Iterator that yields start and end locations of a file chunk of default size 1MB.
    def chunkify(self, f):
        file_end = os.path.getsize(f)

        # Adjust size of chunks so its equal to number of CPUs
        size = 1024 * 1024 * 4 # default was 1MB
        # size = ceil(file_end / self.cpu)

        with open(f, 'rb', size) as file_handle:
            chunk_end = file_handle.tell()
            while True:
                chunk_start = chunk_end
                file_handle.seek(size, 1)
                self.find_end_of_chunk(file_handle)
                chunk_end = file_handle.tell()
                yield chunk_start, chunk_end - chunk_start
                if chunk_end >= file_end:
                    break

    def make_chunks(self, file_handle, size):
        while True:
            chunk = list(islice(file_handle, size))
            if not chunk:
                break
            yield chunk

    # read chunk
    def read_chunk(self, f, chunk_info):
        with open(f, 'rb', 1024 * 1024) as file_handle:
            file_handle.seek(chunk_info[0])
            return file_handle.read(chunk_info[1])

    # End of chunk
    def find_end_of_chunk(self, file_handle):
        line = file_handle.readline()  # incomplete line
        position = file_handle.tell()
        line = file_handle.readline()
        while line and b'runid' not in line:  # find the start of sequence
            position = file_handle.tell()
            line = file_handle.readline()
        file_handle.seek(position)  # revert one line

    def parse_fastq_to_dict(self, l, my_dict, name, flag):
        header, seq, extra, qual = l  # get each component of list in a variable

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
            # phred_list.append(ord(letter))
            phred_list.append(letter)
        average_phred = sum(phred_list) / len(phred_list) - 33

        # GC percentage
        g_count = float(seq.count(b'G'))
        c_count = float(seq.count(b'C'))
        gc = round((g_count + c_count) / float(length) * 100, 1)

        # Channel
        channel = header.split()[3].split(b'=')[1]

        seq = FastqObjects(name, length, flag, average_phred, gc, time_string, channel)

        my_dict[seq_id] = seq

    def parse_fastq_to_dict_islice(self, l, d, name, flag):
        l = map(str.strip, l)
        header, seq, extra, qual = l  # get each component of list in a variable

        # Sequence ID
        seq_id = header.split()[0][1:]

        # Read Time stamp
        time_string = header.split()[4].split('=')[1]
        time_string = parse(time_string)

        # Sequence length
        length = len(seq)

        # Average phred score
        phred_list = list()
        for letter in qual:
            phred_list.append(ord(letter))
            # phred_list.append(letter)
        average_phred = sum(phred_list) / len(phred_list) - 33

        # GC percentage
        g_count = float(seq.count('G'))
        c_count = float(seq.count('C'))
        gc = int(round((g_count + c_count) / float(length) * 100))

        seq = FastqObjects(name, length, flag, average_phred, gc, time_string)
        d[seq_id] = seq

    def parse_fastq_to_dict_islice_pool(self, l, name, flag):
        my_dict = {}
        l = map(str.strip, l)
        header, seq, extra, qual = l  # get each component of list in a variable

        # Sequence ID
        seq_id = header.split()[0][1:]

        # Read Time stamp
        time_string = header.split()[4].split('=')[1]
        time_string = parse(time_string)

        # Sequence length
        length = len(seq)

        # Average phred score
        phred_list = list()
        for letter in qual:
            phred_list.append(ord(letter))
            # phred_list.append(letter)
        average_phred = sum(phred_list) / len(phred_list) - 33

        # GC percentage
        g_count = float(seq.count('G'))
        c_count = float(seq.count('C'))
        gc = int(round((g_count + c_count) / float(length) * 100))

        seq = FastqObjects(name, length, flag, average_phred, gc, time_string)
        my_dict[seq_id] = seq

        return my_dict

    def get_chunk_data(self, f, name, flag, chunk_info):
        my_dict = {}

        while True:
            data = self.read_chunk(f, chunk_info)
            lines = []
            my_data = data.split(b'\n')
            data = None

            for line in my_data:
                if not line:  # end of file?
                    break
                line = line.rstrip()

                if lines and b'runid' in line:  # new entry
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            my_data = None

            return my_dict

    def get_chunk_data_process(self, f, name, flag, q):
        my_dict = {}
        chunk_info = q.get()

        while True:
            data = self.read_chunk(f, chunk_info)
            lines = []
            my_data = data.split(b'\n')
            data = None

            for line in my_data:
                if not line:  # end of file?
                    break
                line = line.rstrip()

                if lines and b'runid' in line:  # new entry
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            my_data = None

            return my_dict

    def get_chunk_data_queue(self, f, name, flag, q):
        my_dict = {}

        while True:
            chunk_info = q.get()
            if chunk_info is None:
                break
            # print(f, chunk_info)  # debug
            data = self.read_chunk(f, chunk_info)
            lines = []
            my_data = data.split(b'\n')
            data = None

            for line in my_data:
                if not line:  # end of file?
                    break
                line = line.rstrip()

                if lines and b'runid' in line:  # new entry
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            my_data = None

            self.sample_dict.update(my_dict)
            q.task_done()

    def get_chunk_data_map(self, f, name, flag, chunk_info):
        my_dict = {}

        while True:
            data = self.read_chunk(f, chunk_info)
            lines = []
            my_data = data.split(b'\n')
            data = None

            for line in my_data:
                if not line:  # end of file?
                    break
                line = line.rstrip()

                if lines and b'runid' in line:  # new entry
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            my_data = None

            return my_dict

    def get_chunk_data_new(self, f, chunk_info, chunk_data):
        data = self.read_chunk(f, chunk_info)
        chunk_data.append(data.split(b'\n'))

        return chunk_data

    def process_chunk(self, name, flag, chunk):
        """
        Put a fastq entry in a list and send it to self.parse_fastq_to_dict
        :param f: file path
        :param name: name of the file, typically 'barcodeXX'
        :param flag: 'pass' or 'fail'
        :param chunk: a tuple containing a multiple of 4 lines
        :return: dictioanry with relevant information
        """
        my_dict = {}

        while True:
            lines = []
            for line in chunk:
                if not line:  # end of file?
                    continue
                line = line.rstrip()

                # if lines and b'runid' in line:  # new entry -> specific to nanopore
                if len(lines) == 4:
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

            return my_dict

    def process_chunk_new(self, chunk, name, flag):
        lines = []
        my_dict = {}
        for line in chunk:
            if not line:  # end of file?
                break
            line = line.rstrip()

            if lines and b'runid' in line:  # new entry
                self.parse_fastq_to_dict(lines, my_dict, name, flag)
                lines = []
            lines.append(line)
        self.parse_fastq_to_dict(lines, my_dict, name, flag)  # For the last entry

        return my_dict

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

        for f in l:
            name = os.path.basename(f).split('.')[0].split('_')[0]  # basename before 1st "_" -> sample name
            # name = file.name.split('/')[-2]  # Parent folder name

            filename = os.path.basename(f)
            filename_no_gz = '.'.join(filename.split('.')[:-1])
            gz_flag = 0

            start_time = time()
            if f.endswith('.fastq.gz'):
                print("Decompressing {} to '/tmp'".format(os.path.basename(f)))
                p1 = subprocess.Popen(['pigz', '-d', '-k', '-f', '-c', f], stdout=subprocess.PIPE)
                (outs, errs) = p1.communicate()
                with open('/tmp/' + filename_no_gz, 'b+w', 1024 * 1024) as tmp_file:
                    tmp_file.write(outs)
                f = "/tmp/" + filename_no_gz
                gz_flag = 1
            if f.endswith(".fastq"):
                pass
            else:
                raise Exception('Wrong file extension. Please use ".fastq" or ".fastq.gz"')

            # get some stats about the input file
            if gz_flag == 1:
                stat_info = os.stat('/tmp/' + filename_no_gz)
            else:
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

            ####################################
            # Multiprocessing -> pool.apply_async

            # pool = mp.Pool(self.cpu)
            # # # lines_per_chunk = self.optimal_line_number(f, self.cpu)
            # lines_per_chunk = 4
            # #
            # jobs = []
            # with open(f, 'rb', 1024 * 1024) as file_handle:
            #     # for chunk in self.make_chunks(file_handle, 40):
            #     #     job = pool.apply_async(self.process_chunk, args=(name, flag, chunk))
            #     #     jobs.append(job)
            #     for chunk in zip_longest(*[file_handle] * lines_per_chunk, fillvalue=None):
            #         job = pool.apply_async(self.process_chunk, args=(name, flag, chunk))
            #         jobs.append(job)
            # #     for chunk_info in self.chunkify(f):
            # #         job = pool.apply_async(self.get_chunk_data, args=(f, name, flag, chunk_info))
            # #         jobs.append(job)
            # # # with open(f, 'rb', 1024 * 1024) as file_handle:
            # # #     for chunk in self.make_chunks(file_handle, lines_per_chunk):
            # # #         job = pool.apply_async(self.process_chunk, args=(name, flag, chunk))
            # # #         jobs.append(job)
            # #
            # output = []
            # for job in jobs:
            #     output.append(job.get())
            #
            # jobs.clear()
            #
            # pool.close()
            # pool.join()
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in output:
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #
            # output.clear()
            #########################################

            ####################################
            # Multiprocessing -> pool.map

            # pool = mp.Pool(self.cpu)
            # # pool = mp.Pool(4)
            #
            # chunk_info_list = []
            # for chunk_info in self.chunkify(f):
            #     chunk_info_list.append(chunk_info)
            # func = partial(self.get_chunk_data, f, name, flag)
            # results = pool.map_async(func, chunk_info_list)
            # # results = pool.starmap_async(func, chunk_info_list)
            #
            # pool.close()
            # pool.join()
            # # pool.terminate()  # Needed to do proper garbage collection?
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in results.get():
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #########################################

            #########################################
            # Process
            # q = mp.Queue()
            # jobs = []
            # for i in range(self.cpu):
            #     p = mp.Process(target=self.get_chunk_data_process, args=(f, name, flag, q,))
            #     jobs.append(p)
            #     p.start()
            #
            # for chunk_info in self.chunkify(f):
            #     q.put(chunk_info)
            #
            # for proc in jobs:
            #     proc.join()
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in return_dict.values():
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #########################################

            #########################################
            # asyncio
            import asyncio


            #########################################

            #########################################
            # islice
            # pool = mp.Pool(self.cpu)
            # # lines_per_chunk = self.optimal_line_number(f, self.cpu)
            # lines_per_chunk = 40
            #
            # jobs = []
            # with open(f, 'rb') as file_handle:
            #     for chunk in zip_longest(*[file_handle] * lines_per_chunk, fillvalue=None):
            #         job = pool.apply_async(self.process_chunk, args=(name, flag, chunk))
            #         jobs.append(job)
            #
            # output = []
            # for job in jobs:
            #     output.append(job.get())
            #
            # jobs.clear()
            #
            # pool.close()
            # pool.join()
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in output:
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #
            # output.clear()
            #########################################

            #########################################
            # pool = mp.Pool(self.cpu)
            # # lines_per_chunk = self.optimal_line_number(f, self.cpu)
            # lines_per_chunk = 40
            #
            # chunk_list = []
            # with open(f, 'rb') as file_handle:
            #     for chunk in zip_longest(*[file_handle] * lines_per_chunk, fillvalue=None):
            #         chunk_list.append(chunk)
            #
            # results = []
            #
            # func = partial(self.process_chunk, name, flag)
            # results = pool.map_async(func, chunk_list, chunksize=1)
            # # results = pool.map(func, chunk_list)
            #
            # pool.close()
            # pool.join()
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in results.get():
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #########################################

            #########################################
            # put chunk in list
            # pool = mp.Pool(self.cpu)
            #
            # chunk_data = []
            # for c in self.chunkify(f):
            #     self.get_chunk_data_new(f, c, chunk_data)
            #
            # jobs = []
            # for chunk in chunk_data:
            #     job = pool.apply_async(self.process_chunk_new, args=(chunk, name, flag,))
            #     jobs.append(job)
            #
            # output = []
            # for job in jobs:
            #     output.append(job.get())
            #
            # pool.close()
            # pool.join()
            #
            # # Merge dictionary
            # reads = 0
            # # print("Merging processes outputs")
            # for dictionary in output:
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge

            #########################################

            #########################################
            # Threading approach

            # # Setup some threads
            # threads = []
            # for i in range(self.cpu):
            #     t = threading.Thread(target=self.get_chunk_data_queue, args=(f, name, flag, self.my_queue))
            #     t.setDaemon(True)
            #     t.start()
            #     threads.append(t)
            #
            # # Feed the threads
            # for chunk_info in self.chunkify(f):
            #     self.my_queue.put(chunk_info)
            #
            # # Wait for queue to be empty
            # self.my_queue.join()
            #
            # #Stop the workers
            # for i in range(self.cpu):
            #     self.my_queue.put(None)
            # for t in threads:
            #     t.join()
            #########################################

            #########################################
            # Line by line approach

            reads = 0
            with open(f, 'rb', 1024 * 1024) as file_handle:
                lines = []
                for line in file_handle:
                    if not line:  # end of file?
                        break
                    line = line.rstrip()
                    if len(lines) == 4:
                        reads += 1
                        self.parse_fastq_to_dict(lines, d, name, flag)
                        lines = []
                    lines.append(line)
                reads += 1
                self.parse_fastq_to_dict(lines, d, name, flag)
            #########################################

            #########################################
            # linear with chunks
            # results = []
            # for chunk_info in self.chunkify(f):
            #     res = self.get_chunk_data(f, name, flag, chunk_info)
            #     results.append(res)
            #
            # # Update self.sample_dict with results from every chunk
            # reads = 0
            # for dictionary in results:
            #     # Count number of reads in the sample
            #     reads += len(dictionary)
            #     d.update(dictionary)  # Do the merge
            #########################################

            #########################################
            # Serial, 4 lines at the time
            # n_lines = 4
            # with open(f, 'r') as file_handle:
            #     for chunk in zip_longest(*[iter(file_handle)] * n_lines, fillvalue=None):
            #         self.parse_fastq_to_dict_islice(chunk, d, name, flag)
            #
            # reads = len(d.keys())
            #########################################

            end_time = time()
            interval = end_time - start_time
            print("took {} ({} reads)".format(self.elapsed_time(interval), reads))

            #  Remove tmp file
            if gz_flag == 1:
                os.remove('/tmp/' + filename_no_gz)

    def parse_file(self, f):
        name = os.path.basename(f).split('.')[0].split('_')[0]  # basename before 1st "_" -> sample name

        # get some stats about the input file
        stat_info = os.stat(f)
        file_size = self.hbytes(stat_info.st_size)

        flag = 'pass'  # Default value
        if 'fail' in f:
            flag = 'fail'  # Check in path for the word "fail"

        # Parse
        my_dict = {}
        with gzip.open(f, 'rb', 1024 * 1024) if f.endswith('gz') else open(f, 'rb', 1024 * 1024) as file_handle:
            lines = []
            for line in file_handle:
                if not line:  # end of file?
                    break
                line = line.rstrip()
                if len(lines) == 4:
                    self.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            self.parse_fastq_to_dict(lines, my_dict, name, flag)

        return my_dict

    def parse_fastq_parallel(self, l, d):
        print("Parsing fastq files...", end="", flush=True)
        start_time = time()

        # Parse the files in parallel
        pool = mp.Pool(self.cpu)

        jobs = []
        for f in l:
            job = pool.apply_async(self.parse_file, [f])
            jobs.append(job)

        results = []
        for j in jobs:
            results.append(j.get())

        pool.close()
        pool.join()
        # pool.terminate()  # Needed to do proper garbage collection?

        # Update self.sample_dict with results from every chunk
        read_counter = 0
        for dictionary in results:
            read_counter += len(dictionary.keys())
            d.update(dictionary)  # Do the merge

        end_time = time()
        interval = end_time - start_time
        print(" took %s for %d reads" % (self.elapsed_time(interval), read_counter))

    # Fastq plots

    def make_fastq_plots(self, d):

        print("\nMaking plots:")

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_reads_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_bp_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        self.plot_reads_vs_bp_per_sample(d)
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

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        self.plot_phred_score_distribution(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        self.plot_length_distribution(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting pores_output_vs_time_total...', end="", flush=True)
        # start_time = time()
        # self.plot_pores_output_vs_time_total(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        self.plot_pores_output_vs_time_all(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_scatter...', end="", flush=True)
        # start_time = time()
        # self.plot_quality_vs_length_scatter(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_length_hex(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting quality_vs_length_kde...', end="", flush=True)
        # start_time = time()
        # self.plot_quality_vs_length_kde(d)
        # # self.test_plot(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        self.plot_channel_output_all(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting gc_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_gc_vs_time(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting gc_vs_length_hex...', end="", flush=True)
        start_time = time()
        self.plot_gc_vs_length_hex(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_all...', end="", flush=True)
        start_time = time()
        self.plot_pores_gc_output_vs_time_all(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_gc_output_vs_time_per_sample...', end="", flush=True)
        start_time = time()
        self.plot_pores_gc_output_vs_time_per_sample(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

    def plot_total_reads_vs_time(self, d):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :return: A png file with the graph
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
        fig.savefig(self.output_folder + "/total_reads_vs_time.png")

    def plot_reads_per_sample_vs_time(self, d):
        """
        Plot yield per sample. Just the pass reads
        :param d: Dictionary
        :return: png file
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
        fig.savefig(self.output_folder + "/reads_per_sample_vs_time.png")

    def plot_bp_per_sample_vs_time(self, d):
        """
        Read length per sample vs time
        :param d: Dictionary
        :return: png file
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
        fig.savefig(self.output_folder + "/bp_per_sample_vs_time.png")

    def plot_total_bp_vs_time(self, d):
        """
        Sequence length vs time
        :param d: Dictionary
        :return: png file
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
        fig.savefig(self.output_folder + "/total_bp_vs_time.png")

    def plot_quality_vs_time(self, d):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :return: png file
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
        fig.savefig(self.output_folder + "/quality_vs_time.png")

    def plot_phred_score_distribution(self, d):
        """
        Frequency of phred scores
        :param d: Dictionary
        :return: png file
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots()

        qual_pass = list()
        qual_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                qual_pass.append(round(seq.average_phred, 1))
            else:
                qual_fail.append(round(seq.average_phred, 1))

        mean_pass_qual = None
        mean_fail_qual = None
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
        name, length, flag, average_phred, gc, time_string
        """

        size_pass = list()
        size_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == 'pass':
                size_pass.append(seq.length)
            else:
                size_fail.append(seq.length)

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
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(self.output_folder + "/length_distribution.png")

    def test_plot(self, d):
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
        g.savefig(self.output_folder + "/quality_vs_length_kde.png")

    def kde2D(self, x, y, bandwidth, xbins=100j, ybins=100j, **kwargs):
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

    def plot_quality_vs_length_kde(self, d):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :return: png file
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
        xx, yy, density = self.kde2D(x, y, 1)

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
            xx, yy, density = self.kde2D(x, y, 1)

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
        g.savefig(self.output_folder + "/quality_vs_length_kde.png")

    def plot_quality_vs_length_hex(self, d):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :return: png file
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
            g.ax_joint.legend(handles=[blue_patch, red_patch], loc='best')
        else:
            g.ax_joint.legend(handles=[blue_patch], loc='best')

        # Add legend to top margin_plot area
        # g.ax_marg_x.legend(('pass', 'fail'), loc='upper right')

        # Set figure size
        g.fig.set_figwidth(8)
        g.fig.set_figheight(4)

        # Save figure to file
        g.savefig(self.output_folder + "/quality_vs_length_hex.png")

    def jointplot_w_hue(self, data, x, y, hue=None, colormap=None,
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

    def plot_quality_vs_length_scatter(self, d):
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
            self.jointplot_w_hue(data=df_concatenated, x='Length (bp)', y='Phred Score',
                                 hue='flag', figsize=(10, 6), fig=fig, colormap=['blue', 'red'],
                                 scatter_kws={'s': 1, 'alpha': 0.1})
        else:
            self.jointplot_w_hue(data=df_pass, x='Length (bp)', y='Phred Score',
                                 hue='flag', figsize=(10, 6), fig=fig, colormap=['blue'],
                                 scatter_kws={'s': 1, 'alpha': 0.1})

        fig.savefig(self.output_folder + "/quality_vs_length_scatter.png")

    def plot_test_old(self, d):
        """
        seaborn jointplot (length vs quality). More manual.
        :param d:
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
        # g.savefig(self.output_folder + "/test.png")
        fig.savefig(self.output_folder + "/test.png")

    def plot_reads_vs_bp_per_sample(self, d):
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
        fig.savefig(self.output_folder + "/reads_vs_bp_per_sample.png")

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
        # fig.savefig(self.output_folder + "/reads_vs_bp_per_sample.png")

        ##########################################
        # df = pd.DataFrame(columns=['Sample', 'Value', 'Info'])
        # for name, size_list in my_sample_dict.items():
        #     df = df.append({'Sample': name, 'Value': sum(size_list), 'Info': 'Total bp'}, ignore_index=True)
        #     df = df.append({'Sample': name, 'Value': len(size_list), 'Info': 'Reads'}, ignore_index=True)
        #
        # g = sns.catplot(x='Sample', y='Value', hue='Info', data=df, kind='bar')
        #
        # plt.tight_layout()
        # g.savefig(self.output_folder + "/reads_vs_bp_per_sample.png")

    def plot_pores_output_vs_time_total(self, d):

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
        fig.savefig(self.output_folder + "/pores_output_vs_time.png")

    def plot_pores_output_vs_time_all(self, d):

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
        # Generate counts for each bin
        hist, edges = np.histogram(time_list_all, bins=x_bins, density=False)
        # Plot the data
        g = sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                        cmap='Greens')  # , marker='o'

        # If not fail, just draw the pass. Else, draw total, fail and pass
        if time_list_fail:  # fail reads might be missing if plotting filtered reads for example.
            # Plot fail
            time_list_fail[:] = [x - time_zero for x in time_list_fail]
            time_list_fail[:] = [x.days * 1440 + x.seconds / 60 for x in time_list_fail]
            time_list_fail[:] = [int(np.round(x)) for x in time_list_fail]
            # nbins = max(time_list_fail) / 15 if max(time_list_fail) % 15 == 0 else int(max(time_list_fail) / 15) + 1
            # x_bins = np.linspace(min(time_list_fail), max(time_list_fail), nbins)
            hist, edges = np.histogram(time_list_fail, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                            cmap='Reds')  # , marker='x'

            # Plot pass - Assume always pass reads present
            time_list_pass[:] = [x - time_zero for x in time_list_pass]
            time_list_pass[:] = [x.days * 1440 + x.seconds / 60 for x in time_list_pass]
            time_list_pass[:] = [int(np.round(x)) for x in time_list_pass]
            # nbins = max(time_list_pass) / 15 if max(time_list_pass) % 15 == 0 else int(max(time_list_pass) / 15) + 1
            # x_bins = np.linspace(min(time_list_pass), max(time_list_pass), nbins)
            hist, edges = np.histogram(time_list_pass, bins=x_bins, density=False)
            sns.scatterplot(data=hist, x_bins=edges, size=3, alpha=0.5, linewidth=0,
                                cmap='Blues')  # , marker='^'

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
        ax.xaxis.set_major_formatter(FuncFormatter(my_formater))
        ax.xaxis.set_major_locator(MultipleLocator(4 * 4))  # 4 block of 15 min per hour. Want every 4 hours

        # Add legend to the plot area
        # https://stackoverflow.com/questions/47391702/matplotlib-making-a-colored-markers-legend-from-scratch

        green_circle = mlines.Line2D([], [], color='green', alpha=0.6, label='All', marker='o',
                                     markersize=5, linestyle='None')
        blue_triangle = mlines.Line2D([], [], color='blue', alpha=0.6, label='Pass', marker='o',
                                      markersize=5, linestyle='None')
        red_x = mlines.Line2D([], [], color='red', alpha=0.6, label='Fail', marker='o',
                              markersize=5, linestyle='None')

        if time_list_fail:
            ax.legend(handles=[green_circle, blue_triangle, red_x], loc='upper right')
        else:
            green_circle = mlines.Line2D([], [], color='green', alpha=0.6, label='Pass', marker='o',
                                         markersize=5, linestyle='None')
            ax.legend(handles=[green_circle], loc='upper right')

        # Add label to axes
        plt.title('Pores output over time')
        plt.ylabel('Reads per 15 minutes')
        plt.xlabel('Sequencing time (hours)')

        plt.tight_layout()  # Get rid of extra margins around the plot
        fig = g.get_figure()  # Get figure from FacetGrid
        fig.savefig(self.output_folder + "/pores_output_vs_time_all.png")

    def plot_channel_output_all(self, d):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :return: png file
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
            layout = self.make_layout(maxval=maxval)
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
        fig.savefig(self.output_folder + "/channel_output_all.png")

    def plot_gc_vs_time(self, d):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :return: png file
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
        fig.savefig(self.output_folder + "/gc_vs_time.png")

    def plot_gc_vs_length_hex(self, d):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :return: png file
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
        g.savefig(self.output_folder + "/gc_vs_length_hex.png")

    def plot_pores_gc_output_vs_time_all(self, d):

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
        fig.savefig(self.output_folder + "/pores_gc_output_vs_time_all.png")

    # @staticmethod
    # def is_prime(n):
    #     return n > 1 and all(n % i for i in islice(count(2), int(sqrt(n) - 1)))

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

    def plot_pores_gc_output_vs_time_per_sample(self, d):
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
        width, height = NanoQC.find_best_matrix(n_sample)
        # print(n_sample, width, height)  # debug

        # Make grid for all samples
        # https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
        fig, ax = plt.subplots(height, width, sharex='col', sharey='row', figsize=(height*5, width*5))
        sample_index = 0
        for i in range(height):
            for j in range(width):
                if sample_index >= len(sample_list):
                    ax[i, j].axis('off')  # don't draw the plot is no more sample for the 'too big' matrix
                    # break
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
        fig.savefig(self.output_folder + "/pores_gc_output_vs_time_all.png")

    # Summary plots

    def parse_summary(self, d):
        """
        Parse "sequencing_summary.txt" file from Albacore
        :param d: Empty summary dictionary
        :return: Dictionary with info about
        """

        print("Parsing summary file...", end='', flush=True)
        start_time = time()
        with open(self.input_summary, 'rb', 1024 * 1024 * 8) as file_handle:
            # fields = list()
            read_counter = 0
            next(file_handle)  # skip first line
            for line in file_handle:
                line = line.rstrip()
                if not line:
                    continue

                fields = line.split(b'\t')

                # Only keep  fields of interest
                length = fields[12]
                if length == b'0':
                    continue  # skip zero-length reads
                seq_id = fields[1]
                channel = fields[3]
                events = fields[6]
                flag = fields[7]
                name = fields[19]
                time_stamp = fields[4]
                average_phred = fields[13]

                seq_summary = SummaryObjects(name, length, channel, events, average_phred, time_stamp, flag)

                d[seq_id] = seq_summary
                read_counter += 1

        # Print read stats
        end_time = time()
        interval = end_time - start_time
        print(" took %s for %d reads" % (self.elapsed_time(interval), read_counter))

    def make_summary_plots(self, d):
        print("\nMaking plots:")

        print('\tPlotting total_reads_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_reads_vs_time_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting total_bp_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_total_bp_vs_time_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_reads_per_sample_vs_time_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting bp_per_sample_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_bp_per_sample_vs_time_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting length_distribution...', end="", flush=True)
        start_time = time()
        self.plot_length_distribution_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting phred_score_distribution...', end="", flush=True)
        start_time = time()
        self.plot_phred_score_distribution_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_length_hex...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_length_hex_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting reads_vs_bp_per_sample...', end="", flush=True)
        start_time = time()
        self.plot_reads_vs_bp_per_sample_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting pores_output_vs_time_all...', end="", flush=True)
        start_time = time()
        self.plot_pores_output_vs_time_all_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting pores_output_vs_time...', end="", flush=True)
        # start_time = time()
        # self.plot_pores_output_vs_time_summary(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        # print('\tPlotting channel_output_total...', end="", flush=True)
        # start_time = time()
        # self.plot_channel_output_total(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))
        #
        # print('\tPlotting channel_output_pass_fail...', end="", flush=True)
        # start_time = time()
        # self.plot_channel_output_pass_fail(d)
        # end_time = time()
        # interval = end_time - start_time
        # print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting channel_output_all...', end="", flush=True)
        start_time = time()
        self.plot_channel_output_all_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

        print('\tPlotting quality_vs_time...', end="", flush=True)
        start_time = time()
        self.plot_quality_vs_time_summary(d)
        end_time = time()
        interval = end_time - start_time
        print(" took %s" % self.elapsed_time(interval))

    def plot_total_reads_vs_time_summary(self, d):
        """
        Plot number of reads against running time. Both Pass and fail reads in the same graph
        :param d: A dictionary to store the relevant information about each sequence
        :return: A png file with the graph
        TODO -> use numpy to handle the plot data, on row per sample?
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()
        t_pass = list()  # time
        t_fail = list()

        for seq_id, seq in d.items():
            t = seq.time_stamp
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/total_reads_vs_time.png")

    def plot_reads_per_sample_vs_time_summary(self, d):
        """
        Plot yield per sample. Just the pass reads
        :param d: Dictionary
        :return: png file
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        # fig, ax = plt.subplots()
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/reads_per_sample_vs_time.png")

    def plot_bp_per_sample_vs_time_summary(self, d):
        """
        Read length per sample vs time
        :param d: Dictionary
        :return: png file
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots(figsize=(10, 6))  # In inches

        # Fetch required information
        my_sample_dict = defaultdict()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/bp_per_sample_vs_time.png")

    def plot_total_bp_vs_time_summary(self, d):
        """
        Sequence length vs time
        :param d: Dictionary
        :return: png file
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()

        # Fetch required information
        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/total_bp_vs_time.png")

        def plot_quality_vs_time_summary(self, d):
            """
            Quality vs time (bins of 1h). Violin plot
            :param d: Dictionary
            :return: png file
            name, length, flag, average_phred, gc, time_string
            """

            fig, ax = plt.subplots(figsize=(10, 4))

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

                df_pass = pd.DataFrame(list(ts_pass3),
                                       columns=['Sequencing time interval (h)', 'Phred score'])  # Convert to dataframe
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
                g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, hue='Flag', split=True,
                                   inner=None)
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
            fig.savefig(self.output_folder + "/quality_vs_time.png")

    def plot_phred_score_distribution_summary(self, d):
        """
        Frequency of phred scores
        :param d: Dictionary
        :return: png file
        # name, length, channel, events, average_phred, time_stamp, flag
        """

        fig, ax = plt.subplots()

        qual_pass = list()
        qual_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
                qual_pass.append(round(float(seq.average_phred), 1))
            else:
                qual_fail.append(round(float(seq.average_phred), 1))

        mean_pass_qual = None
        mean_fail_qual = None
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

    def plot_length_distribution_summary(self, d):
        """
        Frequency of sizes. Bins auto-sized based on length distribution. Log scale x-axis.
        :param d: Dictionary
        :return: png file
        name, length, flag, average_phred, gc, time_string
        """

        size_pass = list()
        size_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
                size_pass.append(int(seq.length))
            else:
                size_fail.append(int(seq.length))

        mean_fail_size = None
        mean_pass_size = None
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
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
        plt.tight_layout()
        fig.savefig(self.output_folder + "/length_distribution.png")

    def plot_quality_vs_length_hex_summary(self, d):
        """
        seaborn jointplot (length vs quality)
        :param d: Dictionary
        :return: png file
        name, length, flag, average_phred, gc, time_string
        """

        sns.set(style="ticks")

        my_dict = dict()
        for seq_id, seq in d.items():
            my_dict[seq_id] = [int(seq.length), float(seq.average_phred), seq.flag.decode('ascii')]

        df = pd.DataFrame.from_dict(my_dict, orient='index', columns=['Length (bp)', 'Phred score', 'flag'])
        df_pass = df.loc[df['flag'] == 'True']
        df_fail = df.loc[df['flag'] == 'False']

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
        g.savefig(self.output_folder + "/quality_vs_length_hex.png")

    def plot_reads_vs_bp_per_sample_summary(self, d):
        # Fetch required information
        my_sample_dict = defaultdict()  # to get the lengths (bp)
        for seq_id, seq in d.items():
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/reads_vs_bp_per_sample.png")

    def plot_pores_output_vs_time_summary(self, d):

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
        fig.savefig(self.output_folder + "/pores_output_vs_time.png")

    def plot_pores_output_vs_time_all_summary(self, d):

        import matplotlib.lines as mlines

        fig, ax = plt.subplots()

        time_list_all = list()
        time_list_pass = list()
        time_list_fail = list()

        for seq_id, seq in d.items():
            time_stamp = seq.time_stamp
            time_list_all.append(time_stamp)
            if seq.flag == b'True':
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
        fig.savefig(self.output_folder + "/pores_output_vs_time_all.png")

    def make_layout(self, maxval):
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

    def plot_channel_output_total(self, d):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :return: png file
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
            # if seq.flag == b'True':
            #     channel_dict[channel_number][0] += 1
            # else:  # seq.flag == b'False':
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
        layout = self.make_layout(maxval=maxval)
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
        fig.savefig(self.output_folder + "/channel_output_total.png")

    def plot_channel_output_pass_fail(self, d):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :return: png file
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
            if seq.flag == b'True':
                channel_dict[channel_number][0] += 1
            else:  # seq.flag == b'False':
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
            layout = self.make_layout(maxval=maxval)
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
        fig.savefig(self.output_folder + "/channel_output_pass_fail.png")

    def plot_channel_output_all_summary(self, d):
        """
        https://github.com/wdecoster/nanoplotter/blob/master/nanoplotter/spatial_heatmap.py#L69
        https://bioinformatics.stackexchange.com/questions/745/minion-channel-ids-from-albacore/749#749

        :param d: Dictionary
        :return: png file
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
            if seq.flag == b'True':
                channel_dict[channel_number][0] += 1
            else:  # seq.flag == b'False':
                channel_dict[channel_number][1] += 1

        # convert to Pandas dataframe
        df = pd.DataFrame.from_dict(channel_dict, orient='index', columns=['Pass', 'Fail'])
        df_all = pd.DataFrame()
        df_all['All'] = df['Pass'] + df['Fail']
        df_pass = df[['Pass']]  # The double square brakets keep the column name
        df_fail = df[['Fail']]

        # Plot
        fig, axs = plt.subplots(nrows=3, figsize=(6, 12))

        for i, my_tuple in enumerate([(df_all, 'All', 'Greens'), (df_pass, 'Pass', 'Blues'), (df_fail, 'Fail', 'Reds')]):
            my_df = my_tuple[0]
            flag = my_tuple[1]
            cmap = my_tuple[2]

            maxval = max(my_df.index)  # maximum channel value
            layout = self.make_layout(maxval=maxval)
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
        fig.savefig(self.output_folder + "/channel_output_all.png")

    def plot_quality_vs_time_summary(self, d):
        """
        Quality vs time (bins of 1h). Violin plot
        :param d: Dictionary
        :return: png file
        name, length, flag, average_phred, gc, time_string
        """

        fig, ax = plt.subplots(figsize=(10, 6))

        ts_pass = list()
        ts_fail = list()
        for seq_id, seq in d.items():
            if seq.flag == b'True':
                ts_pass.append(tuple((seq.time_stamp, round(float(seq.average_phred), 1))))
            else:
                ts_fail.append(tuple((seq.time_stamp, round(float(seq.average_phred), 1))))

        if ts_pass:
            ts_pass[:] = [tuple(((float(x) / 3600), y)) for x, y in ts_pass]  # Convert to hours (float)
            ts_pass.sort(key=lambda x: x[0])  # Sort according to first element in tuple
            ts_pass[:] = [tuple((int(np.round(x)), y)) for x, y in ts_pass]  # Round hours

            df_pass = pd.DataFrame(list(ts_pass), columns=['Sequencing time interval (h)', 'Phred score'])  # Convert to dataframe
            df_pass['Flag'] = pd.Series('pass', index=df_pass.index)  # Add a 'Flag' column to the end with 'pass' value

        if ts_fail:
            ts_fail[:] = [tuple(((float(x) / 3600), y)) for x, y in ts_fail]
            ts_fail.sort(key=lambda x: x[0])
            ts_fail[:] = [tuple((int(np.round(x)), y)) for x, y in ts_fail]

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
            g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, hue='Flag', split=True, inner=None)
            g.figure.suptitle('Sequence quality over time')
        elif ts_pass:
            g = sns.violinplot(x='Sequencing time interval (h)', y='Phred score', data=data, inner=None)
            g.figure.suptitle('Sequence quality over time (pass only)')
        else:  # elif ts_fail:
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
        fig.savefig(self.output_folder + "/quality_vs_time.png")


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

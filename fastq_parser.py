import os
import gzip
from concurrent import futures
from dateutil.parser import parse
import multiprocessing as mp
from itertools import repeat, islice
import functions


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


class FastqParser(object):
    @staticmethod
    def hbytes(num):
        for x in ['bytes', 'KB', 'MB', 'GB']:
            if num < 1024.0:
                return "%3.1f%s" % (num, x)
            num /= 1024.0
        return "%3.1f%s" % (num, 'TB')

    @staticmethod
    def make_chunks(file_handle, size):
        while True:
            chunk = list(islice(file_handle, size))
            if not chunk:
                break
            yield chunk

    @staticmethod
    def read_entry(entry, name, flag):
        my_dict = dict()
        lines = list()
        for line in entry:
            line = line.rstrip()
            lines.append(line)
            if len(lines) == 4:
                my_dict.update(FastqParser.single_fastq_entry_to_dict(lines, name, flag))
                lines = list()
        return my_dict

    @staticmethod
    def single_fastq_entry_to_dict(fastq_entry, name, flag):
        header, seq, extra, qual = fastq_entry  # get each component of list in a variable

        # Split header line to get items
        items = header.split()

        # Store items info in dictionary
        read_dict = dict()
        seq_id = items[0][1:]  # remove the leading "@"
        read_dict['seq_id'] = seq_id
        for i in items[1:]:
            key, value = i.split('=')
            read_dict[key] = value

        # Read length
        length = len(seq)

        # Average phred score
        phred_list = [ord(letter) - 33 for letter in qual]
        average_phred = functions.compute_average_quality(phred_list, length)  # cython

        # %GC
        g_count = float(seq.count('G'))
        c_count = float(seq.count('C'))
        gc = round((g_count + c_count) / float(length) * 100, 2)

        # Time stamp
        try:
            time_string = parse(read_dict['start_time'])
        except KeyError:
            time_string = ''

        # Channel
        try:
            channel = read_dict['ch']
        except KeyError:
            channel = ''

        seq = FastqObjects(name, length, flag, average_phred, gc, time_string, int(channel))
        my_dict = dict()
        my_dict[seq_id] = seq

        return my_dict

    @staticmethod
    def iterate_fastq_parallel(input_fastq, cpu, parallel):
        # Name
        name = os.path.basename(input_fastq).split('.')[0].split('_')[0]
        name = name.replace('_pass', '')
        name = name.replace('_fail', '')

        # Flag
        flag = 'pass'  # Default value
        if 'fail' in input_fastq:
            flag = 'fail'  # Check in path for the word "fail"

        # Chunk fastq files and run chunks in parallel
        fastq_dict = dict()
        with gzip.open(input_fastq, "rt") if input_fastq.endswith('.gz') else open(input_fastq, "r") as f:
            pool = mp.Pool(int(cpu / parallel))
            jobs = [pool.apply_async(FastqParser.read_entry, [chunk, name, flag]) for chunk in FastqParser.make_chunks(
                f, 4000)]
            results = [j.get() for j in jobs]
            pool.close()
            pool.join()
            pool.terminate()  # Needed to do proper garbage collection?

            # Update self.sample_dict with results from every chunk
            for d in results:
                fastq_dict.update(d)  # Do the merge

        return fastq_dict

    @staticmethod
    def parallel_process_fastq(fastq_list, cpu, parallel):
        fastq_dict = dict()
        with futures.ProcessPoolExecutor(max_workers=parallel) as executor:
            for results in executor.map(FastqParser.iterate_fastq_parallel, fastq_list, repeat(cpu), repeat(parallel)):
                fastq_dict.update(results)

        return fastq_dict

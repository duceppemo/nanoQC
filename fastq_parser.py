import os
import gzip
from concurrent import futures
from dateutil.parser import parse
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
    def parse_fastq_to_dict(l, my_dict, name, flag):
        header, seq, extra, qual = l  # get each component of list in a variable

        # Guppy v3.1.5
        # @9ab0b453-9b87-4c7f-925d-e700919d955f runid=2223752dbeaf60b17877fd21486a51f7974218ca sampleid=PSM read=5089 ch=226 start_time=2019-05-30T00:14:52Z
        # @092edb67-90bb-4abe-976e-30ee11798c38 runid=2223752dbeaf60b17877fd21486a51f7974218ca read=46236 ch=388 start_time=2019-05-31T01:50:16Z flow_cell_id=FAK82962 protocol_group_id=PSMWGS-2 sample_id=PSM barcode=unclassified

        # Sequence ID
        items = header.split()
        seq_id = items[0]

        # Read Time stamp
        time_string = next((t for t in items if b'start_time=' in t), None)
        time_string = time_string.split(b'=')[1]
        time_string = parse(time_string)

        # Sequence length
        length = len(seq)

        # Average phred score
        phred_list = [letter - 33 for letter in qual]
        average_phred = functions.compute_average_quality(phred_list, length)  # cython

        # GC percentage
        g_count = float(seq.count(b'G'))
        c_count = float(seq.count(b'C'))
        gc = round((g_count + c_count) / float(length) * 100, 2)

        # Channel
        channel = next((c for c in items if b'ch=' in c), None)
        channel = channel.split(b'=')[1]

        seq = FastqObjects(name, length, flag, average_phred, gc, time_string, channel)

        my_dict[seq_id] = seq

    @staticmethod
    def parse_file(f):
        name = os.path.basename(f).split('.')[0]
        name = name.replace('_pass', '')
        name = name.replace('_fail', '')

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
                    FastqParser.parse_fastq_to_dict(lines, my_dict, name, flag)
                    lines = []
                lines.append(line)
            FastqParser.parse_fastq_to_dict(lines, my_dict, name, flag)

        return my_dict

    @staticmethod
    def parse_fastq_parallel(l, d, cpu):
        with futures.ProcessPoolExecutor(max_workers=cpu) as pool:
            results = pool.map(FastqParser.parse_file, l, chunksize=1)

        # Update dictionary with results from every chunk
        read_counter = 0
        for dictionary in results:
            read_counter += len(dictionary.keys())
            d.update(dictionary)  # Do the merge
        return read_counter

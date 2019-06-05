from time import time


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


class SummaryParser(object):
    @staticmethod
    def parse_summary(input_summary, d):
        """
        Parse "sequencing_summary.txt" file demultiplexed while basecalling with Albacore
        :param d: Empty summary dictionary
        :return: Dictionary with info about
        """
        read_counter = 0
        with open(input_summary, 'rb', 1024 * 1024 * 8) as file_handle:
            # fields = list()
            next(file_handle)  # skip first line
            for line in file_handle:
                line = line.rstrip()
                if not line:
                    continue

                fields = line.split(b'\t')

                # Only keep  fields of interest
                length = fields[11]
                if length == b'0':
                    continue  # skip zero-length reads
                seq_id = fields[1]
                channel = fields[3]
                events = fields[6]
                flag = fields[7]
                name = fields[19]  # barcode
                time_stamp = fields[4]
                average_phred = fields[12]

                seq_summary = SummaryObjects(name, length, channel, events, average_phred, time_stamp, flag)

                d[seq_id] = seq_summary
                read_counter += 1
        return read_counter

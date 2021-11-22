import pandas as pd


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
    def parse_summary(input_summary):
        df = pd.read_csv(input_summary, sep='\t')

        # Keep only column of interest
        df = df.loc[:, ('barcode_arrangement', 'sequence_length_template', 'passes_filtering',
                        'mean_qscore_template', 'start_time', 'channel')]
        # Rename columns
        # ['Name', 'Length', 'Flag', 'Qual', 'Time', 'Channel']
        df.rename(columns={'barcode_arrangement': 'Name',
                           'sequence_length_template': 'Length',
                           'passes_filtering': 'Flag',
                           'mean_qscore_template': 'Qual',
                           'start_time': 'Time',
                           'channel': 'Channel'}, inplace=True)

        # Rename True for pass and False for Fail
        flag_dict = {True: 'pass', False: 'fail'}
        df['Flag'] = df['Flag'].map(flag_dict)

        return df

    @staticmethod
    def parse_summary_og(input_summary, d):
        """
        Parse "sequencing_summary.txt" file demultiplexed while basecalling with Albacore
        :param input_summary: path to the summary file
        :param d: Empty summary dictionary
        :return: Dictionary with info about
        """
        read_counter = 0
        with open(input_summary, 'rb', 1024 * 1024 * 8) as file_handle:
            foi = dict()  # fields of interest
            header = file_handle.readline()  # skip first line
            # find which column is which field
            items = header.split(b'\t')

            # find the index of the fields of interest.
            # Should make it easier if new columns are added or shuffled in new versions of guppy,
            # unless nanopore decides to change names of the column header... who knows!
            foi['length'] = items.index(b'sequence_length_template')
            foi['seq_id'] = items.index(b'read_id')
            foi['channel'] = items.index(b'channel')
            foi['events'] = items.index(b'num_events')
            foi['flag'] = items.index(b'passes_filtering')
            try:
                foi['name'] = items.index(b'barcode_arrangement')
            except ValueError:
                # The basecalling was done without demultiplexing
                foi['name'] = None
            foi['time_stamp'] = items.index(b'start_time')
            foi['average_phred'] = items.index(b'mean_qscore_template')

            for line in file_handle:
                line = line.rstrip()
                if not line:
                    continue

                fields = line.split(b'\t')

                # Only keep  fields of interest
                length = fields[foi['length']]
                if length == b'0':
                    continue  # skip zero-length reads
                seq_id = fields[foi['seq_id']]
                channel = fields[foi['channel']]
                events = fields[foi['events']]
                flag = fields[foi['flag']]
                if foi['name']:
                    name = fields[foi['name']]  # barcode
                else:
                    name = b'unclassified'
                time_stamp = fields[foi['time_stamp']]
                average_phred = fields[foi['average_phred']]

                seq_summary = SummaryObjects(name, length, channel, events, average_phred, time_stamp, flag)

                d[seq_id] = seq_summary
                read_counter += 1
        return read_counter
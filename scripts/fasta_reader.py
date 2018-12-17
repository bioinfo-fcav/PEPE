class FastaRecord(object):
    def __init__(self, title, sequence):
        self.title = title
        self.sequence = sequence

def read_fasta_record(infile):
    # The first line in a FASTA record is the title line.
    # Examples:
    # >third sequence record
    # >gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene
    line = infile.readline()

    if not line:
        # End of file
        return None

    # Double-check that it's a FASTA file by looking for the '>'.
    if not line.startswith(">"):
        raise TypeError("Not a FASTA file: %r" % line)

    # The title starts after the '>' and doesn't include any
    # trailing whitespace.
    title = line[1:].rstrip()

    # Read the sequence lines up to the blank line.
    sequence_lines = []
    while 1:
        # I know I'm at the end of the sequence lines when I reach the
        # blank line or the end of the file.  I can simplify the test
        # by calling rstrip and checking for the empty string.  If it's
        # a sequence line I'll still call rstrip to remove the final
        # newline and any trailing whitespace so I'll go ahead and
        # always rstring the line.
        line = infile.readline().rstrip()
        if line == "":
            # Reached the end of the record or end of the file
            break
        sequence_lines.append(line)

    # Merge the lines together to make a single string containing the sequence
    # (This is faster than doing "sequence = sequence + line" if there are
    # more than a few lines)
    sequence = "".join(sequence_lines)
    
    return FastaRecord(title, sequence)

def read_fasta_records(input_file):
    records = []
    while 1:
        record = read_fasta_record(input_file)
        if record is None:
            break
        records.append(record)
    return records


##### Self-test code

test_data = [
    ("first sequence record",
     "TACGAGAATAATTTCTCATCATCCAGCTTTAACACAAAATTCGCA"),
    ("second sequence record",
     "CAGTTTTCGTTAAGAGAACTTAACATTTTCTTATGACGTAAATGA" +
     "AGTTTATATATAAATTTCCTTTTTATTGGA"),
    ("third sequence record",
     "GAACTTAACATTTTCTTATGACGTAAATGAAGTTTATATATAAATTTCCTTTTTATTGGA" +
     "TAATATGCCTATGCCGCATAATTTTTATATCTTTCTCCTAACAAAACATTCGCTTGTAAA"),
    ]

# Test the ability to read a single record at a time
def test_single_record():
    infile = open("simple.fasta")
    for (title, sequence) in test_data:
        record = read_fasta_record(infile)
        if title != record.title:
            raise AssertionError("title %r != %r" % (title, record.title))
        if sequence != record.sequence:
            raise AssertionError("sequence %r != %r" % (sequence,
                                                        record.sequence))
    record = read_fasta_record(infile)
    if record is not None:
        raise AssertionError(repr(record))

# Test the ability to read all the records into a list
def test_records():
    infile = open("simple.fasta")

    count = 0
    for record in read_fasta_records(infile):
        title, sequence = test_data[count]
        if title != record.title:
            raise AssertionError("title %r != %r" % (title, record.title))
        if sequence != record.sequence:
            raise AssertionError("sequence %r != %r" % (sequence,
                                                        record.sequence))
        count = count + 1

    if count != len(test_data):
        raise AssertionError("Read %d records, expected %d" %
                             (count, len(test_data)))

# The master test function
def test():
    test_single_record()
    test_records()
    
if __name__ == "__main__":
    test()
    print "Done."

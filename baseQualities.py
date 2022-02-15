import pysam

#bam = "/Users/jimin/PEPPER_clr/HG003.rerun_pacbio.clr.phased.haplotagged.chr20.bam"
#samfile = pysam.AlignmentFile(bam, "rb")
#print("count", samfile.count()) #count 295902

class CommandLine():
    def __init__(self, inOptions = None):
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        self.parser.add_argument('--inputBam', type=str, default=None, required=True, help='bam file to adjust quality scores')
        self.parser.add_argument('--outputBam', type=str, default=None, required=True, help='name of output bam file')
        self.parser.add_argument('--score', type=int, default=0, required=True, help='score to adjust all base quality scores to')

        if inOptions is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOptions)


def replaceQualityScores(bam, outputfilename, score):
    """
    Replace bam file base quality scores with specified score using PySam.
    :param bam: Input bam file
    :param outputfilename: Output file with adjusted base qualities

    :param score: Score to replace all base quality scores to
    :return: None
    """
    pysamfile = pysam.AlignmentFile(bam, "rb")
    #print("header", pysamfile.header)
    # read is AlignedSegment
    outfile = pysam.AlignmentFile(outputfilename, "wb", template=pysamfile)
    for read in pysamfile.fetch():
        #header = read.header
        length = read.query_length
        #print(read.query_qualities)
        #print(read.query_qualities)
        read.query_qualities = [score]*length
        outfile.write(read)
    outfile.close()
    pysamfile.close()



def main(myCommandLine=None):
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(myCommandLine)

    inputBam = myCommandLine.args.inputBam
    outputBam = myCommandLine.args.outputBam
    score = myCommandLine.args.score
    replaceQualityScores(inputBam, outputBam, score)


if __name__ == '__main__':
    main()



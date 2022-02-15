import pysam

class CommandLine():
    def __init__(self, inOptions = None):
        import argparse
        
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
   
    # read is AlignedSegment
    outfile = pysam.AlignmentFile(outputfilename, "wb", template=pysamfile)
    for read in pysamfile.fetch():
        
        length = read.query_length
       
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



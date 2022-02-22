from multiprocessing import Process, Pool, cpu_count
import time
from timeit import default_timer as timer
import pysam
import logging


class CommandLine():
    def __init__(self, inOptions = None):
        import argparse
        self.parser = argparse.ArgumentParser(
            description='This program adjusts base qualities of bam files to a score set by user.',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s --inputBam --outputBam --score'
        )
        self.parser.add_argument('--inputBam', type=str, default=None, required=True, help='bam file to adjust quality scores')
        self.parser.add_argument('--outputBam', type=str, default=None, required=True, help='name of output bam file')
        self.parser.add_argument('--score', type=int, default=0, required=True, help='score to adjust all base quality scores to')

        if inOptions is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOptions)


def replaceQualityScores(read):
    """
    Replace bam file base quality scores with specified score using PySam.
    :param bam: Input bam file
    :param outputfilename: Output file with adjusted base qualities

    :param score: Score to replace all base quality scores to
    :return: None
    """
    #print("header", pysamfile.header)
    # read is AlignedSegment
    score = 0
    #header = read.header
    length = read.query_length
    #print(read.query_qualities)
    #print(read.query_qualities)
    read.query_qualities = [score]*length
    return read



def main(myCommandLine=None):
    start = timer()
    print(f'starting computations on {cpu_count()} cores')

    if myCommandLine is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(myCommandLine)

    inputBam = myCommandLine.args.inputBam
    pysamfile = pysam.AlignmentFile(inputBam, "rb")
    alignedSegList = pysamfile.fetch()
    reads = list(x.query_sequence for x in alignedSegList)
    print("type", type(reads))
    print(reads)

    outputBam = myCommandLine.args.outputBam
    outfile = pysam.AlignmentFile(outputBam, "wb", template=pysamfile)
    score = myCommandLine.args.score

    with Pool() as pool:
        res = pool.starmap(replaceQualityScores, reads)
        outfile.write(res)

    outfile.close()
    pysamfile.close()

    end = timer()
    print(f'elapsed time: {end - start}')


if __name__ == '__main__':
    main()



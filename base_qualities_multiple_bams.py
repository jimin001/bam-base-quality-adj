import sys
import concurrent.futures
from datetime import datetime
import pysam

class CommandLine():
    def __init__(self, inOptions = None):
        import argparse
        self.parser = argparse.ArgumentParser(
            description='This program adjusts base qualities of bam files to a score set by user.',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s --inputBam --outputBam --score'
        )
        self.parser.add_argument('--input', type=str, default=None, required=True, help='bam file to adjust quality scores')
        self.parser.add_argument('--output', type=str, default=None, required=True, help='name of output bam file')
        self.parser.add_argument('--fasta', type=str, default=None, required=True, help='fasta file')
        self.parser.add_argument('--score', type=int, default=0, required=True, help='score to adjust all base quality scores to')

        if inOptions is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOptions)



def bam_by_chunk(bam_file_name, output_file_name, contig, score, process_id):

    bam = pysam.AlignmentFile(bam_file_name, "rb")
    outputfilename = output_file_name + contig
    outfile = pysam.AlignmentFile(outputfilename, "wb", template=bam)
    reads = bam.fetch(contig)

    adjusted_reads = []
    for read in reads:
        length = read.query_length

        read.query_qualities = [score]*length
        # store updated read to list

        # write adjusted read straight to outfile
        outfile.write(read)
        #adjusted_reads.append(read)
    outfile.close()
    bam.close()

    return process_id

def get_list_of_contig(fasta_file_name):
    """
    Create dictionary of contig names and contig lengths.
    :param fasta_file_name:
    :param query_contig_name:
    :return:
    """
    # https://pysam.readthedocs.io/en/latest/api.html
    contig_list = []
    fasta = pysam.FastaFile(fasta_file_name)

    contig_names = fasta.references
    contig_lengths = fasta.lengths

    for contig_name in contig_names:
        contig_list.append(contig_name)
        #contig_list.append(contig_name)
    return contig_list



def run_parallel(total_processes, bam_file_name, output_file_name, score, contig_list):
    total_reads = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=total_processes) as executor:
        futures = [executor.submit(bam_by_chunk, bam_file_name, output_file_name, contig, score, process_id) for process_id, contig in enumerate(contig_list)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                process_id = fut.result()
                # reads_with_fixed_qualities = fut.result()
                #total_reads += total_reads_in_chunk
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESS " + str(process_id) + " FINISHED SUCCESSFULLY.\n")
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None

    # for fixed_read in reads_with_fixed_qualities:
    #     output_bam.write(fixed_read)
    return total_reads


if __name__ == '__main__':
    myCommandLine = CommandLine()
    bam_file = myCommandLine.args.input
    #bam_file = '/Users/jimin/PEPPER_clr/HG003.pacbio.clr.haplotagged.chr20.bam'
    output_name = myCommandLine.args.output
    #fasta_file = "/Users/jimin/PEPPER_clr/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
    fasta_file = myCommandLine.args.fasta
    score = myCommandLine.args.score

    contig_list = get_list_of_contig(fasta_file)

    updated_contig_dict = run_parallel(4, bam_file, output_name, score, contig_list)


# 141993
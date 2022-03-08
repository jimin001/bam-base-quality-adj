import sys
import concurrent.futures
from datetime import datetime
import pysam


def function(bam_file_name, region):
    bam = pysam.AlignmentFile(bam_file_name, "rb")
    reads = bam.fetch(region)
    total_reads = 0
    for read in reads:
        total_reads += 1
    return total_reads


def function_by_chunk(bam_file_name, region, region_start, region_end, process_id):
    bam = pysam.AlignmentFile(bam_file_name, "rb")
    reads = bam.fetch(region, region_start, region_end)
    total_reads = []
    for read in reads:
        read_length = read.query_length
        total_reads.append(read)
        #total_reads += 1
    return total_reads, process_id


def get_len_of_contig(fasta_file_name, query_contig_name):
    # https://pysam.readthedocs.io/en/latest/api.html
    fasta = pysam.FastaFile(fasta_file_name)

    contig_names = fasta.references
    contig_lengths = fasta.lengths

    for contig_name, contig_length in zip(contig_names, contig_lengths):
        if contig_name == query_contig_name:
            return contig_length
    return 0


def run_parallel(total_processes, bam_file_name, query_region, chunks):
    total_reads = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=total_processes) as executor:
        futures = [executor.submit(function_by_chunk, bam_file_name, query_region, start_pos, end_pos, process_id) for process_id, (start_pos, end_pos) in enumerate(chunks)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                total_reads_in_chunk, process_id = fut.result()
                # reads_with_fixed_qualities = fut.result()
                total_reads += total_reads_in_chunk
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESS " + str(process_id) + " FINISHED SUCCESSFULLY.\n")
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None

    # for fixed_read in reads_with_fixed_qualities:
    #     output_bam.write(fixed_read)
    return total_reads


if __name__ == '__main__':
    bam_file = '/Users/jimin/PEPPER_clr/HG003.pacbio.clr.haplotagged.chr20.bam'
    fasta_file = "/Users/jimin/PEPPER_clr/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
    #bam_file = "/Users/kishwar/Kishwar/data/users/common/bam/HG002_guppy_507_2_GRCh38_pass.chr20.30x.bam"
    region = "chr20"
    #fasta_file = "/Users/kishwar/Kishwar/data/users/common/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    length_of_contig = get_len_of_contig(fasta_file, "chr20")
    chunk_size = 1000000
    chunks_of_chr20 = [(start_pos, min(start_pos + chunk_size-1, length_of_contig) ) for start_pos in range(0, length_of_contig, 1000000)]

    total_reads = run_parallel(4, bam_file, region, chunks_of_chr20)
    print(total_reads)
    # total_reads = 0
    # for chunk in chunks_of_chr20:
    #     start_pos, end_pos = chunk
    #     total_reads += function_by_chunk(bam_file, region, start_pos, end_pos)
    # print(total_reads)

# 141993
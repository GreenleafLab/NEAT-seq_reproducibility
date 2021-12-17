# Runs basically instantaneously if the bam files are already indexed
import pysam
import os

snp_chr = "chr17"
snp_base = 75518526 - 1

exclude_bam_flag = 1804
include_bam_flag = 2

# Note: The SRA data is not yet public at this time, but once it is this section
# can be updated to include the bam downloads
input_data_path = "/oak/stanford/groups/wjg/amyfchen/GEO_submission_June2021"
output_path = "fig3_correlation_analysis/outputs/"
os.makedirs(output_path, exist_ok=True)

bam_lane1 = pysam.AlignmentFile(f"{input_data_path}/CD4_lane1/lane1_atac_possorted_bam.bam", "rb")
bam_lane2 = pysam.AlignmentFile(f"{input_data_path}/CD4_lane2/lane2_atac_possorted_bam.bam", "rb")

out_lane1 = open(f"{output_path}/lane1_fragments.tsv", "w")
out_lane2 = open(f"{output_path}/lane2_fragments.tsv", "w")

def write_snp_fragments(bam_in, tsv_out, snp_chr, snp_base):
    for pileupcolumn in bam_in.pileup(snp_chr, snp_base, snp_base+1):
        if pileupcolumn.pos != snp_base:
            continue
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            # query position is None if is_del or is_refskip is set.
            read = pileupread.alignment

            # Mapping quality filters
            if read.flag & exclude_bam_flag != 0 or read.flag & include_bam_flag != include_bam_flag:
                continue
            if read.mapq < 30:
                continue
            
            start = read.reference_start
            end = read.next_reference_start
            start, end = min(start, end), max(start, end)
            
            seq = read.query_sequence[pileupread.query_position]
            qual = read.query_qualities[pileupread.query_position]
            cell = read.get_tag("CB") if read.has_tag("CB") else "UNKONWN_CELL"
            print(snp_chr, start+4, end-5, cell, seq, qual, sep="\t", file=tsv_out)


write_snp_fragments(bam_lane1, out_lane1, snp_chr, snp_base)
write_snp_fragments(bam_lane2, out_lane2, snp_chr, snp_base)
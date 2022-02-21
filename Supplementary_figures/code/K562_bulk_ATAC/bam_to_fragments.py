import sys
import pysam

# Convert bam file to a fragment file format, while adding +4/-5 coordinate adjustment
# Read bam input from stdin, output bed output to stdout.
# Designed to put a bulk bam file into the same format as a single cell sample.
input = pysam.AlignmentFile(sys.stdin, "rb")
sample_name = sys.argv[1]

prev_read = next(input)
while True:
    read = next(input, None)
    if read is None or prev_read is None:
        break
    if read.query_name != prev_read.query_name:
        print("Read with missing mate:", prev_read.query_name, file=sys.stderr)
        prev_read = read
        continue

    if prev_read.is_reverse:
        prev_read, read = read, prev_read

    chromosome = read.reference_name
    start = prev_read.reference_start + 4
    end = start + prev_read.template_length - 5
    cell_barcode = sample_name
    print(chromosome, start, end, cell_barcode, sep="\t")
    prev_read = next(input, None)


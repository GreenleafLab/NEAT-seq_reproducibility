import gzip
import os
import http.cookiejar
import urllib.request
import shutil
import sys

# Download all the data required for reproducing the figures of the NEAT-seq paper

# Note: The SRA data is not yet public at this time, which means that the
# bam data required for the SNP analysis is not publicly accessible yet

gse = "GSE178707"
token = input("Enter reviewer token: ")

geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&token={token}"

cj = http.cookiejar.CookieJar() 
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
opener.open(geo_url)

data_files = [
    "GSM5396328_barnyard_atac_fragments.tsv.gz",
    "GSM5396329_barnyard_raw_feature_bc_matrix.h5",
    "GSM5396330_ADT_counts_lane1.csv.gz",
    "GSM5396330_corr_gene_adt.csv.gz",
    "GSM5396330_corr_peak_adt.csv.gz",
    "GSM5396330_putative_targets.csv.gz",
    "GSM5396332_lane1_atac_fragments.tsv.gz",
    "GSM5396333_CD4_RNA_counts.rds.gz",
    "GSM5396333_lane1_features.tsv.gz",
    "GSM5396334_ADT_counts_lane2.csv.gz",
    "GSM5396336_CD4_Peak_matrix.rds.gz",
    "GSM5396336_lane2_atac_fragments.tsv.gz",
]

gzipped_rds = [
    "GSM5396333_CD4_RNA_counts.rds.gz",
    "GSM5396336_CD4_Peak_matrix.rds.gz",
] 

# Download files, and for .rds.gz files additionally decompress to just .rds
output_dir = "geo_download"
os.makedirs("geo_download", exist_ok=True)
for i, f in enumerate(data_files):
    output_path = f if f not in gzipped_rds else f.replace(".gz", "")
    
    print("Downloading", f, f"({i+1}/{len(data_files)})", end="")
    if os.path.isfile(os.path.join(output_dir, output_path)):
        print(" Skipping (file exists)")
        continue
    else:
        print("")
    
    accession = f[:f.index("_")]

    try:
        download = opener.open(f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file&file={f}")
    except urllib.error.HTTPError:
        print("Error downloading file (check reviewer token)")
        sys.exit()

    if f in gzipped_rds:
        download = gzip.open(download)
    
    shutil.copyfileobj(
        download,
        open(os.path.join(output_dir, output_path), "wb")
    )


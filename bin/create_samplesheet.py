import pandas as pd
from glob import glob

# Get list of JSON files and extract sample names
samples = sorted(list({i.split("/")[-1].split(".")[0] for i in glob("../taxprofiler/taxprofiler_results/fastp/*.gz")}))

samples = sorted(list({i.split("/")[-1].split(".")[0].split(".merged")[0].replace("_1", "").replace("_2", "") for
                       i in glob("../taxprofiler/taxprofiler_results/fastp/*.fastq.gz")}))
# Extract only the prefix part (before "_") of the sample names
# samples = sorted(list({"_".join(sample.split("_")[:-1]) if sample.startswith("Arc") else sample.split("_")[0] for sample in samples}))

# Get list of fastq_1 and fastq_2 files
fastq_1 = sorted([i for i in glob("../taxprofiler/taxprofiler_results/fastp/*1.fastp.fastq.gz")])
fastq_2 = sorted([i for i in glob("../taxprofiler/taxprofiler_results/fastp/*2.fastp.fastq.gz")])
# runmerg = sorted([i for i in glob("../taxprofiler/taxprofiler_results/run_merging/*.fastq.gz")])

# Create a DataFrame
try:
    samplesheet = pd.DataFrame({"sample": samples, "fastq_1": fastq_1, "fastq_2": fastq_2})
except:
    breakpoint()

# Save DataFrame to a CSV file
samplesheet.to_csv("maple_samplesheet.csv", index=False)


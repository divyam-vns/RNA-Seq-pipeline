# Install FastQC if not already installed
sudo apt-get install fastqc

# Create a directory for QC reports
mkdir -p /path/to/qc_reports

# Run FastQC on all FASTQ files
fastqc /path/to/raw_reads/*.fastq.gz -o /path/to/qc_reports
# Review the generated HTML reports to identify any issues that may need to be addressed before alignment.

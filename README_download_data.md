
# Download Data Script - `download_data.sh`

## Description

This script is designed to download sequencing data from the Sequence Read Archive (SRA) using the SRA Toolkit's `prefetch` and `fasterq-dump` commands. It automates the process of downloading specific SRA accession numbers, converting them into FASTQ files, and storing them in a specified directory.

### Features:
- Downloads SRA data for a predefined list of SRA accession numbers.
- Converts the downloaded SRA data to FASTQ format using `fasterq-dump`.
- Stores the FASTQ files in the `/home/projects/data` directory.

## Requirements

This script requires the following:
- **SRA Toolkit**: Ensure that `prefetch` and `fasterq-dump` are available on your system.
- **Linux environment**: The script is designed for use on a Linux system.
  
To check if the SRA toolkit is installed, run:
```bash
which prefetch
which fasterq-dump
```

## Usage

1. **Make the script executable**:
   After downloading or creating the script file, make it executable by running the following command:
   ```bash
   chmod +x /home/projects/download_data.sh
   ```

2. **Run the script**:
   To execute the script and download the data, run:
   ```bash
   /home/projects/download_data.sh
   ```

3. **Check the output**:
   The script will download the data and save it as FASTQ files in the `/home/projects/data` directory.

## Customization

- You can modify the list of **SRA accession numbers** in the script to download different datasets.
  - Example: To add new accessions, update the `ACCESSIONS` array in the script:
    ```bash
    ACCESSIONS=("SRR3734796" "SRR3734797" "SRR3734798" "SRR1234567")
    ```

## Notes

- **Storage**: Ensure that you have enough disk space in the `/home/projects/data` directory to store the downloaded data.
- **Network**: This script requires a stable internet connection to download the data from SRA.

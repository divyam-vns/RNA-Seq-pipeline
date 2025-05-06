# RNA-Seq Environment Setup Script

This repository contains `setup_env.sh`, a shell script that sets up the environment for RNA-Seq analysis by configuring paths to bioinformatics tools such as:

- STAR
- SAMtools
- sratoolkit (prefetch, fasterq-dump)
- Minimap2
- Cufflinks
- IRFinder
- psiclass
- rmats-turbo

Location:
The tools are expected to be located inside `/home/software/`.

What `setup_env.sh` Does:
- Exports PATH variables to include installed tools
- Makes the environment ready to use for RNA-Seq data analysis

How to Use:

1. Make the script executable:
   ```bash
   chmod +x /root/setup_env.sh
   ```

2. Run the script:
   ```bash
   source /root/setup_env.sh
   ```

This will update the PATH for the session, allowing you to use the tools.

Tip:
To persist these settings across terminal sessions, add the contents of `setup_env.sh` to your `.bashrc` or `.profile`.

Author:
Divya Mishra  
GitHub Profile: https://github.com/divyam-vns

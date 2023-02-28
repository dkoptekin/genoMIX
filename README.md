![logo](test/genoMIX_fig00.png)

A simple python script to generate artificial admixed individual(s) from given source individuals or populations by using the Eigenstrat dataset. For more detail about the Eigenstrat format, you can check [here](https://reich.hms.harvard.edu/software/InputFileFormats).

The tool splits snp file into chunks with a fixed number of SNPs, then creates an artificial admixed genome by randomly selecting chunks from given source populations based on their admixture proportions.


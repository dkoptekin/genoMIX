![logo](test/genoMIX_fig00.png)

A simple python script to generate artificial admixed individual(s) from given source individuals or populations by using the Eigenstrat dataset. For more detail about the Eigenstrat format, you can check [here](https://reich.hms.harvard.edu/software/InputFileFormats).

The tool splits the ```input.snp``` file into chunks with a fixed number of SNPs then creates an artificial admixed genome by randomly selecting chunks from one individual of given source populations based on their admixture proportions.

The script can generate admixed individuals as a 2, 3 and/or 4-way mixture using parameters specified in a YAML file like below.

```example.yaml```

```yaml
data_prefix: /path/to/input/eigenstrat/input # Eigenstrat path without extension
output_prefix: /path/to/output/eigenstrat/output
data_type: packed # packed, unpacked
chunk_size: 1000
2-way:
  - sources:
    - Anatolia_EN
    - IronGates_HG
  - rates:
    - [0.70, 0.30]
    - [0.60, 0.40]
3-way:
  - sources:
    - Anatolia_EN
    - IronGates_HG
    - CHG
  - rates:
    - [0.34, 0.33, 0.33]
4-way:
  - sources:
    - Anatolia_EN
    - IronGates_HG
    - CHG
    - Levant_N
  - rates:
    - [0.34, 0.22, 0.22, 0.22]

```

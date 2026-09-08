# Software environment

The analyses require R 4.4 or later. Package names are listed in `DESCRIPTION`; each stage also checks its direct dependencies before reading controlled inputs.

Every locked analysis stage writes `sessionInfo.txt` beside its results. The session record, input/output MD5 manifests and stage locks form the software and data lineage for the executed analysis.

The genomic workflow additionally uses LDSC, FUMA, MAGMA, PLINK and SMR. Record the executable version, command line, reference panel and input checksum in the corresponding local results directory. Provider software and licensed reference files are not redistributed here.

For a new computing environment:

1. install R 4.4 or later;
2. install the CRAN and Bioconductor packages listed in `DESCRIPTION`;
3. install GenomicSEM, LAVA, cfdr.pleio and ldscr from their documented distributions;
4. configure external executables in `config/config.R`;
5. run the required stage from the repository root;
6. retain the generated `sessionInfo.txt`, lock and manifest with the analysis results.

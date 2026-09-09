[DATA_ACCESS.md](https://github.com/user-attachments/files/31987535/DATA_ACCESS.md)
# Data access

This repository distributes analysis code and aggregate reporting templates. Participant-level cohort data, genotypes, controlled summary statistics, licensed molecular resources and reference panels remain with their original providers.

## Ageing cohorts

| Resource | Access route | Repository use |
|---|---|---|
| CHARLS | https://charls.charlsdata.com/ | Longitudinal, mortality and five-year prediction analyses |
| ELSA | https://www.elsa-project.ac.uk/accessing-elsa-data | Longitudinal, mortality, prediction sensitivity and dementia-related analyses |
| HRS | https://hrs.isr.umich.edu/data-products | Longitudinal, mortality and prediction development/internal validation |
| MHAS | https://www.mhasweb.org/ | Longitudinal, mortality and external prediction evaluation |
| SHARE | https://share-eric.eu/data/data-access | Longitudinal, mortality, external prediction and dementia-related analyses |

## Genomic and molecular resources

| Resource | Access route | Repository use |
|---|---|---|
| IEU OpenGWAS | https://gwas.mrcieu.ac.uk/ | Input GWAS discovery and metadata |
| FUMA | https://fuma.ctglab.nl/ | Locus definition and functional annotation |
| FinnGen | https://www.finngen.fi/en/access_results | External endpoint mapping |
| Taiwan Precision Medicine Initiative | https://tpmi.ibms.sinica.edu.tw/ | External endpoint mapping under applicable governance |
| GTEx | https://gtexportal.org/home/datasets | Tissue eQTL follow-up |
| eQTL Catalogue | https://www.ebi.ac.uk/eqtl/Data_access/ | Regional eQTL follow-up |
| MetaBrain | https://www.metabrain.nl/ | Brain eQTL colocalisation |
| 1000 Genomes | https://www.internationalgenome.org/data | Linkage-disequilibrium reference |

## Local data contract

Required file schemas are documented in `data/README.md`. Copy `config/config.example.R` to `config/config.R` and point each configuration variable to an approved local resource. `config/config.R`, `results/`, controlled inputs and row-level intermediate files are excluded from Git.

## Redistribution

Do not commit participant identifiers, dates, row-level longitudinal records, genotype files, controlled GWAS or QTL files, credentials, access tokens or licensed reference panels. Aggregate tables, figures, source data and audit manifests may be shared only when permitted by the relevant data-use agreements.

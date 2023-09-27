# This workflow is deprecated and has been replaced by [PacificBiosciences/HiFI-human-WGS-WDL](https://github.com/PacificBiosciences/HiFI-human-WGS-WDL).

# Disclaimer

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.

# PacBio Human WGS Workflow - snakemake implementation

This repo contains a data analysis pipeline to comprehensively detect and prioritize variants in human genomes using PacBio HiFi reads. It consists of three [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows designed to run sequentially, with PacBio HiFi BAMs or FASTQs as the primary input.

Documentation for the PacBio Human WGS Workflow (snakemake implementation) is available [here](Tutorial.md).

## Contributors

- Juniper Lake ([@juniper-lake](https://github.com/juniper-lake))
- William Rowell ([@williamrowell](https://github.com/williamrowell))
- Aaron Wenger ([@amwenger](https://github.com/amwenger))

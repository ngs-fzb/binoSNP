<img src="https://github.com/ngs-fzb/binoSNP/blob/master/images/binoSNP_logo.png" height="153" width="400">


*****

SNP detection based on binomial test procedure. Perl scripts with R integration and usage of bam-readcount algorithm. 

<p align="center">
 <img src="https://github.com/ngs-fzb/binoSNP/blob/master/images/workflow.web.png">
</p>

*****
 
The script accepts any BAM-file as input, but ideally duplicates (PCR artefacts) have been removed and base quality scores have been recalibrated. Additionally, the script requires an interval list where the positions to be examined are named as well as a RefAlt table defining reference and the alternative allele for those positions. Together with binoSNP we provide one example list specifically for the analysis of Tuberculosis bacteria. The list (Resisnps) includes positions known to be associated with drug resistance. As a first step the `bam-readcount`algorithm from Larson [(bam-readcount)](https://github.com/genome/bam-readcount) is executed to extract information about the number and quality of reference and alternative alleles at the positions named in the interval list and stores this information in a text file. In a second step the resulting txt-file is read into R and for each position a p-value is calculated by using the binomial test procedure. In the next step a table is produced containing all information including the calculated p-value for each position named in the interval list. The last step implies the user-defined filtering, e.g. report variants with a p-value below 5 % (standard value for statistical significance).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Requirements

You need to have the following programs installed and available in your PATH to run binoSNP, for installation please refer to the respective manuals:

* R | [R project page](https://www.r-project.org/)

* Perl | [Perl homepage](https://www.perl.org/)
  * Perl module Statistics::R | install with `cpan Statistics::R`

* bam-readcount | [bam-readcount homepage](https://github.com/genome/bam-readcount)

### Installing binoSNP
The script does not need to be specifically installed or compiled. Once the repository is downloaded or cloned (`git clone https://github.com/ngs-fzb/binoSNP`) you can run binoSNP from the console.
You may need to change the permissions of the binoSNP executable with `chmod 755 binoSNP`. Optionally add binoSNP to your path or run it directly from the binoSNP directory.


## Usage

```
binoSNP --help
```
Will print the following help message:
```
binoSNP 1.0.0 - Copyright (C) 2019  Viola Dreyer, Christian Utpatel
   
   [USAGE]: binoSNP [--OPTION PARAMETER] <.bam file>
   
   Available OPTIONS and default PARAMETERS:
   -i [--interval]      List of intervals to be analyzed
                        Default [Resisnps_Master.v28.interval_list.tsv]

   -m [--mut]           Mutation table (RefAlt table)
                        Default [Resisnps_Master.v28_RefuAlt.tsv]

   -o [--outdir]        Output directory
                        Default [./Low_Freq]

   -r [--ref]           Reference sequence used for aligment in fasta format
                        Default [M._tuberculosis_H37Rv_2015-11-13.fasta]

   -p [--pvalue]        p-value used to filter the results
                        Default [0.05]

   -h [--help]          This help message

   -v [--version]       Version of binoSNP
   ```
You can either run a single file by using:
```
binoSNP /path/to/your/filename.bam
```
Or you run all bamfiles of your directory as a batch:
```
binoSNP path/to/your/bamfiles/*.bam
```

Specify the output directory with `-o path\to\output` otherwise the results will be written to `./Low_Freq`at the location from which the script was invoked.


In case you need other positions analyzed, you can create and use your own intervals. The list of intervals should follow the structure below:
```
M.tuberculosis_H37Rv \t start1 \t stop1 \n
M.tuberculosis_H37Rv \t start2 \t stop2 \n
M.tuberculosis_H37Rv \t start3 \t stop3 \n
M.tuberculosis_H37Rv \t start4 \t stop4 \n
```
For annotation of the positions, the user should also provide a mutation list in the following format:
```
Pos \t REF \t ALT \t Annotation \n
761140 \t A \t G \t RMP resistance \n
```
Usage of an own list:
```
binoSNP -i path/to/your/interval.list -m path/to/your/mutations.list
```



## Authors

* **Viola Dreyer** - *Initial work* - [binoSNP](https://github.com/ngs-fzb/binoSNP)
* **Christian Utpatel** - *Code contribution* 
* Stefan Niemann - *Head*


## License
Copyright (C) 2019 Viola Dreyer, Christian Utpatel, Stefan Niemann. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see [LICENCES](https://www.gnu.org/licenses/gpl-3.0).

## Publication  
In review.

## Acknowledgments

Parts of this work have been supported by the German Center for Infection Research (DZIF).


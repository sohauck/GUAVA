## Introduction
GUAVA is a simple application with two independent parts, one for processing and one for visualisation, to rapidly analyse large MLST datasets.
The name is short for Gene-by-Gene Unified Analysis and Visualising Application.  
It can calculate GC content, locus lengths (minimum, maximum and average and alignment lengths), as well as a handful of measures of diversity and selection. 
It is designed to have as few dependencies as possible, to run quickly even on desktop computers, and to export the results in a portable format. 

## What is necessary for these scripts to work correctly
+ Perl: tested on 5.18.2, including List::Util and LWP::Simple; to download, access https://www.perl.org/get.html
+ If alignments are requested, MAFFT: tested on v7.221; download from http://mafft.cbrc.jp/alignment/software/ 
+ If visualisation is run locally: 
  + R: tested on 3.3.2; to download, access https://cran.r-project.org/
  + R package ggplot2: tested on 2.2.1; to download, run ```install.packages("ggplot2")``` from within R
  + R package shiny: tested on 1.0.4; to download, run ```install.packages("shiny")``` from within R

## To use:
##### Processing phase (Perl)
+ Start the Perl script in the command line, with "perl" followed by the file address. You can choose to use the command line options listed below, or simply answers the prompts given by the script. 
##### Visualising phase (R shiny)
+ To run the visualisation app locally, run "R -e shiny::runApp('/*directory*/')" on the command line, with *directory* replaced by the address of where the *app.R* file is contained. Alternatively, load the script on an R editor such as R Studio and launch directly from the editor. 
+ To access the online visualisation app, visit https://shauck-zoo.shinyapps.io/guava/. 

## Options which may be added to command line:
+ **-help**
Entering this command will exit the script and call up the Usage instructions.
+ **-table**
The table containing allele designations, as a tab-separated table of loci vs. isolates, with allele IDs in the cells. Multiple alleles can be represented with semi-colons, for example *"1;5"*.  
+ **-locuscat**
A tab-separated file with locus names in the first column and categories on the second, used in making graphs.
+ **-transpose**
'Yes' if the table has loci as columns and isolates as rows; 'No' otherwise.  
+ **-dboption**
One of three options for how to access sequence *A* if a directory of FASTA sequences already exists, *B* if grabbing a complete directory from PubMLST, *C* if grabbing selectively from PubMLST
+ **-FASTA**
If FASTA reference sequences already exist (Option *A*), the directory with complete FASTA files. Files must be in the format *locusname.FAS*, for example MYCO001234.FAS  
+ **-dbname**
If grabbing FASTA sequences from PubMLST (Options *B* or *C*), the name of the database, usually in the format "pubmlst\_*genusname*\_seqdef".
+ **-out**
The directory where all results will be saved.
+ **-dup**
The minimum frequency required for an allele to be included in the analysis. Default is '1', where all alleles are included.
+ **-mafft**
This can be used to add more parameters to MAFFT, for example **"--auto"**.


## Parameters in the Results table:
+ **Locus**: The name of locus as it appears in the allele designation table. 
+ **Missing**: The number of isolates for which no allele designations was present in the table. 
+ **MultipleCopies**: The number of isolates in which the locus had at least two tagged alleles.
+ **CountNuc**: The count of unique nucleotide sequences for specified locus.  
+ **CountAA**: The count of unique amino acid sequences for specified locus.
+ **MinLength**: The length in nucleotide letters of shortest nucleotide sequence.
+ **MaxLength**: The length in nucleotide letters of longest nucleotide sequence.
+ **AvgLength**: The average length in nucleotide letters of unique nucleotide sequences.
+ **AllelicDiv**: The number of unique nucleotide sequences found divided by average length (greater if more diverse). 
+ **RatioCount**: The ratio of unique nucleotide to unique amino acid sequences.
+ **AliLenNuc**: The length of the alignment of all sequences in this locus in nucleotide format. 
+ **AliLenAA**: The length of the alignment of all sequences in this locus in amino acid format. 
+ **VSitesNuc**: The proportion of sites in the alignment of nucleotide sequences which had any variation. 
+ **VSitesAA**: The proportion of sites in the alignment of amino acid sequences which had any variation. 
+ **RatioVS**: The ratio of proportions of unique sites in nucleotide to amino acid sequences. 


## Features of the Visualisation Shiny app:
+ **Dynamic plot**: All plots, labels and axes change immediately to reflect changes in the Control Panel options.
+ **Choose plotting options**: Drop-down menus allow for multiple option for the horizontal and vertical axes as well as the point colours.
+ **Flexible designation cut-off**: A slider on the left-hand control panels that removes loci that have designations in fewer than a certain percentage of isolates. 
+ **Labelling loci**: Selecting loci by dragging and double-clicking allows them to be "labelled", causing their name to appear on the graph, and allowing 
+ **Removing loci**: Labelled loci can be removed from the analysis with a single click, allowing the analysis to be customised directly
+ **Interactive table**: An interactive table allows for filtering and sorting by any parameters, including label status. 


## Known Issues
+ The MAFFT alignment tool will not accept file addresses with spaces 
+ R may encounter errors if unusual characters are used in locus or catergory names (for example, a single apostrophe, ')  

# MutantHunter README

MutantHunter is a pipeline to identify NLR-type resistance genes using RenSeq and EMS mutagenesis screens.

In a nutshell, you have a plant with a single resistance gene. You do an EMS mutant screen and challenge individuals from independent M2 families with the pathogen. In rare cases you find loss of resistance. You take DNA from the susceptible mutants and DNA from the wildtype. 
RenSeq is enrichment sequencing targeted on NLR genes. This reduces the complexity of sequencing.
RenSeq is performed individually on each mutant and the wildtype. A denovo assembly is performed on the wildtype data. This assembly is then used as a reference for mapping the data from mutants as well as the wildtype as positive control. The assembly is annotated with NLR-Parser to filter for contigs associated with NLRs. Contigs of the assembly are also aligned to the bait sequences of the RenSeq library. Data from mapping and the annotation is fed into the MutantHunter and candidate contigs are identified based on the number of mutants that have a mutation in the contig.


## Prerequisites
### MEME suite version 4.9.1
The MEME suite is available at [http://meme.nbcr.net/meme/](http://meme.nbcr.net/meme/)

Please note that the most actual version of meme is not compatible with NLR Parser. Use meme 4.9.1.

Don't worry about setting up the Apache webserver. You just need MAST, so the quick install is sufficient. 


### JRE 1.6
Make sure you have the Java Runtime Environments 1.6 or higher. Download from [http://java.com](http://java.com)

### BWA
Download from [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)

### Samtools
Download from [http://samtools.sourceforge.net/](http://samtools.sourceforge.net/)

### Blast+
Download Blast+ from [NCBI](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### A *denovo* assembly software.
Use your favourite assembler to assemble the rawdata of your wildtype. Note that the wrong toold might result in contigs representing a consensus of two or more NLRs. For Illumina data we had nice results with [CLC assembly cell](http://www.clcbio.com/products/clc-assembly-cell/) and default parameters. Please note that this one is a comercial software. Feel free to experiment with free software.

## Preprocessing

### Clean raw data
Use [sickle](http://bioinformatics.ucdavis.edu/software/) or whatever tool you like to clean your rawdata.

Example

```
sickle pe -t sanger -f read1.fq -r read2.fq -o read1.clean.fq -p read2.clean.fq -s useless.fq
``` 

### De novo assembly of wildtype

Use your favourite assembler to assemble the rawdata of your wildtype. Note that the wrong toold might result in contigs representing a consensus of two or more NLRs. For Illumina data we had nice results with [CLC assembly cell](http://www.clcbio.com/products/clc-assembly-cell/) and default parameters.

### Mappings

Map (cleaned) rawdata of your wild type and the mutants against the denovo assembly. We reccomend using [bwa](http://bio-bwa.sourceforge.net/) and [samtools](http://samtools.sourceforge.net/).

```
bwa index assembly.fasta
bwa aln assembly.fasta read1.fastq > read1.aln
bwa aln assembly.fasta read2.fastq > read2.aln
bwa sampe assembly.fasta read1.aln read2.aln read1.fastq read2.fastq > raw.sam
samtools view -f2 -Shub -o raw.bam raw.sam
samtools sort raw.bam sorted
samtools rmdup sorted.bam rmdup.bam
samtools index rmdup.bam
samtools faidx assembly.fasta
samtools mpileup -f assembly.fasta -BQ0 rmdup.bam > pileup.txt

```

### Annotation of denovo assembly
Use the NLR-Parser version2 to generate an xml file with the motif annotations. Download the NLR-Parser.jar from MutantHunter release.
The meme.xml is [here](https://github.com/steuernb/MutantHunter/)

```
java -jar NLR-Parser.jar -x meme.xml -y /path/to/meme/bin/mast -i assembly.fasta -c assembly.mast.xml

```

### Define on-target regions in your assembly

Use a set of NLR sequences, e.g. the bait sequences to annotate those regions on contigs that are associated with NLRs. The bait library can be found here: [Triticea_RenSeq_Baits_V1.fasta](https://github.com/steuernb/MutantHunter/)


```
makeblastdb -in baits.fasta -dbtype nucl
blastn -query assembly.fasta -db baits.fasta -outfmt 5 -out assembly_vs_baits.blastn.xml

```


## Hunter pipeline



### Parse pileups
Parse individual pileups and convert to an xml format. This has been separated from the actual MutantHunter because the pielups can be parsed in parallel. Download Pileup2XML.jar from MutantHunter release.

```
java -jar Pileup2XML.jar -m assembly.mast.xml -i pileup.txt -o pileup.xml

```

**Please not that you have to add the `-w` option for running Pileup2XML on wild type**



### MutantHunter

The MutantHunter will provide the candidate contigs. You will have to allocate a bit of memory to the java vm. To add e.g. 16 Gb or RAM you would write `java -Xmx16000M -jar ...`

```
java -jar MutantHunter.jar -w wiltype.pileup.xml -m mutant1.pileup.xml mutant2.plieup.xml [...] -b assembly_vs_baits.blastn.xml -o output.txt 

```


#### Input parameters
 
parameter | argument        | description
---       |   ---           | ---
**-w**    | *STR*           | The xml generated from the wildtype pileup using Pileup2XML.jar
**-m**    | *STR* [*STR*]+  | The xml files of each mutant generated using Pileup2XML.jar
**-b**    | *STR*           | The alignment of contigs to bait sequences.
**-o**    | *STR*           | Outputfile
**-a**    | *float*         | Maximum reference allele frequency of a SNP to be reported. Default is 0.1
**-c**    | *int*           | Minimum coverage for position to be regarded. Default is 10
**-n**    | *int*           | Minimum number of mutants to report a contig. Default is 2
**-z**    | *int*           | Number of coherent positions with zero coverage to call a deletion mutant. Default is 50


## Tips

### Manual validation
Load your bam file in a genome browser, e.g. [Savant](http://www.genomesavant.com/p/home/index) and check your candidate contigs


## Contact
If there are any issues with the tool or if you would like to collaborate with us, please don't hesitate to contact [us](mailto:burkhard.steuernagel@jic.ac.uk).
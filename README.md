#AnnoSNP

##Description
AnnoSNP is a light weight tool for predicting the effect of SNPs on genes, which was designed for newly asssembled genomes.

###Advatages of AnnoSNP:

    1. You do not need to build a database for fasta file.
    2. Besides Chromosome name, loci, ReferenceGenotype and SNPGenotype, no additional information is needed.


##Prerequisite
Perl version (>5.12.0) is required


##Install
$ git clone https://github.com/lhui2010/AnnoSNP

##Usage
perl AnnoSNP.pl gff fasta snp

Using the test dataset
perl AnnoSNP.pl data/test.gff data/test.fasta  data/test.snp

##Input
GFF3 formated annotation (must be consistent with the fasta file)
Fasta formated sequence file
slim VCF formated SNP


slim VCF format:

ChromosomeName | Coordinate | Reference |  AlternativeSNP
-------------  |---------   |---------- | -------------  
chromosome01   |  412917 | G    |   T
chromosome01   |  422620 | T    |   G
chromosome01   |  424391 | A    |   G
chromosome01   |  429519 | T    |   G
chromosome01   |  433888 | A    |   G


##Output
___*.snp_spec___ contained the annotation result for each SNP

| GeneName |  Chromosome | Strand | Loci | Reference | SNP  | LociCDS | Phase |  ReferenceCodon | SNPCodon  |  MutationType |
| --------------   | ---------   | ------ | ---- | --------  | --- | ------- | ----- | --------------- | --------  |  ------------- |
| Os01t0104000-01 | chromosome01  |    -   |  205867  |  A  |   G  |   1251  |    2  |   TAT |  TAC |  Synonymous |
| Os01t0108000-01  | chromosome01   |   +  |   422620  |  T  |   G  |   621 |  2   |  CAT  | CAG |  Radical |
| Os01t0108000-01 |  chromosome01  |    +  |   424391  |  A   |  G   |  1032   |   2  |   AAA |  AAG |  Synonymous |
| Os01t0108500-00 |  chromosome01   |   +   |  460096  |  T  |   C   |  690 |  2  |   CTT |  CTC  | Synonymous |
| Os01t0110400-01 |  chromosome01  |    +  |   549270  |  A   |  G  |   619  | 0  |   ACA  | GCA  | Radical |



___*.snp_stat___ is the statistics of mutation on each genes

GeneName  |  Nonsynonymous |  Synonymous | Radical
---------  | ---------   | ---------- | -----------  
Os01t0104000-01 | 0   | 1 |   0
Os01t0108000-01 | 1  |  1  |  1
Os01t0108500-00 | 0 |   1  |  0


##Benchmark

Tested on Quad-Core AMD Opteron (tm) Processor 8374 HE 800MHz using rice resequencing data. The average speed was `19ms/SNP`.


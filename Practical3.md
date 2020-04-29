# Gentetics&Genomics 2020

# Introduction to genomics data editing and manipulation

## MAIN GOALS 
The objective of these two sessions is to learn some basic but useful practical tricks for genome-scale data mining. To do that, you will practice with a few of the many bioinformatics utilities commonly used to genome data handling and manipulation. The lesson is divided in three main parts. The first part (session 1, the weed of May 4th), will be dedicated to install all required packages and software, to obtain the genomics data files to work with, and to see some simple examples of how to extract features and retrieve valuable information form these files. The second part (homework), will consist of solving a practical exercise (see below) based on the same tools and data retrieved in the first session; students must deliver (June 5th) a small report explaining briefly the workflow and the scripts used to find the solution. Finally, in the last part (session 2, June 10th) we will provide the appropriate feedback for the correct solution of the exercise.

## SESSION 1
## Installing packages and software utilities
The first step consists in to download and install the packages and bioinformatics tools necessary to work with genomics data. It is worth noting that here you will use only a very small representation of the enormous variety of packages and software utilities (both in `Python` and other programming languages) currently available to work with genome sequences, annotations and variants.


Installing `gffutils` package (you must have installed a `Python` distribution). In the command-line terminal, type::

```bash
pip install gffutils
```
		
> `gffutils` is a `Python` package for working with GFF and GTF files in a hierarchical manner [https://pythonhosted.org/gffutils/index.html]. 

Installing `BioPython` (required for some `gffutils` utilities):

```bash
pip install biopython
```

> `BioPython` is a set of freely available tools for biological computation written in `Python` (https://biopython.org/).

Installing `VCFtools` from GitHub platform (required to work with genomic variants):

Linux:

```bash
sudo apt install git
```

Mac (you must to have installed Homebrew; https://brew.sh/):

```bash
brew install git
```

Linux and Mac (it may be necessary to install [autoconf](autoconf.md)):

```bash
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make
sudo make install
cd ..
```

> `VCFtools` (https://github.com/vcftools/vcftools) is a program package designed for working with VCF files (https://vcftools.github.io/specs.html) [You many need sudo permissions to run make install].

## Genomics data files
Before working with genomic files, we need to install `wget`, a software package for retrieving remote files, and `samtools`, a suit of programs for interacting with high-throughput sequencing data:

Linux:

```bash
sudo apt install wget samtools
```

Mac: 

```bash
brew install wget samtools
```

You will download genome sequences (in FASTA format), and feature annotations (in GFF3 format; https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) of the human chromosome 19 (=GRCh37 assembly) from Ensembl FTP server (http://www.ensembl.org/info/data/ftp/index.html) [assembly version must be exactly the same for both files!]. Create a new working directory and use the wget and bgzip (from `samtools` suite) commands to download and decompress the genomic files:

```bash
mkdir myworkdir
cd myworkdir
wget -qO- ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz | bgzip -d > Homo_sapiens.GRCh37.dna.chromosome.19.fa
wget -qO- ftp://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chromosome.19.gff3.gz | bgzip -d > Homo_sapiens.GRCh37.87.chromosome.19.gff3
```

You can take a first look the first part (=the head; in this case the first 1500 lines) of these files using the UNIX command head:

```bash
head -1500 Homo_sapiens.GRCh37.87.chromosome.19.gff3
head -1500 Homo_sapiens.GRCh37.dna.chromosome.19.fa
```

## Examples

**EXAMPLE 1**

Imagine that you are interested in obtaining all the proteins (the products of all mRNAs) encoded by a gene of interest (e.g. the human gene APOE). [For simplicity, we will assume that we know that the gene is on the human chromosome 19]. You can write a Python script for this task:

Create a new text file (within the /myworkdir folder) and rename it as “APOE_proteins.py”. Then add this expression to the first line of this file:

```python
#!/usr/bin/env python
```

An easy way to quickly and efficiently to store and access the structural annotations of the chromosome 19 is to create a database of the features and relationships in the GFF3 file (Homo_sapiens.GRCh37.87.chromosome.19.gff3). We can use the package gffutils to build this database. Add these lines to the script to import this Python package:

```python
import gffutils
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import BiopythonWarning
```

> In addition to `gffutils`, we also import some practical `Biopython` methods to manage sequences. 

Before creating the database, and for practical purposes, it is recommended to store data file names into variables:
	
```python
myGFF="Homo_sapiens.GRCh37.87.chromosome.19.gff3"
myFASTA="Homo_sapiens.GRCh37.dna.chromosome.19.fa"
myGENE="APOE"
```

Now we are ready to build the database with the genomic features of the human chromosome 19. To do that, we will use the function `gffutils.create_db()`

```python
db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)
```

> The database (db) is stored in the memory of the computer (`':memory:'`). The `merge_strategy` command ensures that each line of the GFF file maintain separate accessible keys (recommendable here to have all genomics information in the database). 

To learn some useful commands that allow us to access the database for the features of the gene of interest we can open a `Python` interactive session (>>>).  Open a command-line terminal in your computer and then type in `python`. Then type the above commands to create de data base:

```
>>> import gffutils
>>> myGFF="Homo_sapiens.GRCh37.87.chromosome.19.gff3"
>>> myFASTA="Homo_sapiens.GRCh37.dna.chromosome.19.fa"
>>> myGENE="APOE"
>>> db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)
```
	
Each annotation feature (gene, mRNA, exon, CDS, …), in specified in a line of the GFF file and has in its 9th column a list of attributes defining the feature, in the format tag=value. To access the different attributes of a feature, new need to know the first of these tags, the tag ID. For a gene feature, for instance, the tag ID is in the form ID=gene:ENSG00000116032 (the attribute identifier of the gene is the Ensembl gene ID). For transcripts, ID=transcript:ENST00000263097, and for CDS, ID=CDS:ENSP00000263097. Note that the attribute identifiers of these features are the Ensembl transcript and protein IDs, respectively. 

To access the gene ENSG00000116032, type:

```
>>> mygene=db['gene:ENSG00000116032']
```

To see the complete feature, print the variable:

```
>>> print(mygene)
19	ensembl_havana	gene	1000418	1009731	.	+	.	ID=gene:ENSG00000116032;Name=GRIN3B;biotype=protein_coding;description=glutamate receptor%2C ionotropic%2C N-methyl-D-aspartate 3B [Source:HGNC Symbol%3BAcc:16768];gene_id=ENSG00000116032;logic_name=ensembl_havana_gene;version=5
```

	
To directly access to some the fields of this feature, you can type:

```
>>> mygene.source
'ensembl_havana'
>>> mygene.start
1000418
>>> mygene.end
1009731
>>> mygene.strand
'+'
```

You can also access to the attributes stored in the 9th column, using the function `feature.attributes()`:

```
>>> mygene.attributes['Name']
['GRIN3B']
>>> mygene.attributes['description']
'glutamate receptor, ionotropic, N-methyl-D-aspartate 3B [Source:HGNC Symbol;Acc:16768]']
```

The same for a transcript feature:

```
>>> mytr=db['transcript:ENST00000607316']
>>> print(mytr)
19	havana	mRNA	1011851	1013339	.	-	.	ID=transcript:ENST00000607316;Name=TMEM259-018;biotype=protein_coding;version=1;Parent=gene:ENSG00000182087;havana_transcript=OTTHUMT00000471311;havana_version=1;transcript_id=ENST00000607316
>>> mytr.start
1000418
>>> mytr.end
1009731
```

As you know, some features, such as transcripts, CDS or exons, have parent features. Transcripts, for example, have parent genes, while exons and CDS have parent transcripts. These ontologies are also specified in attributes and can be accessed with the same function:

```
>>> mytr.attributes[‘Parent’]
['gene:ENSG00000182087']
```

Nevertheless, sometimes you need to access all, or part of the children features of the same parent at the same time. The package `gffutils` implements specific functions to do that: 

```
>>> mygene=db['gene:ENSG00000064666']
>>> print(mygene)
19	ensembl_havana	gene	1026298	1039068	.	+	.	ID=gene:ENSG00000064666;Name=CNN2;biotype=protein_coding;description=calponin 2 [Source:HGNC Symbol%3BAcc:2156];gene_id=ENSG00000064666;logic_name=ensembl_havana_gene;version=10
>>> for t in db.children(mygene, featuretype='mRNA', order_by='start'):
...		print t
...
19	ensembl_havana	mRNA	1026298	1039068	.	+	.	ID=transcript:ENST00000263097;Name=CNN2-001;biotype=protein_coding;version=4;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420293;havana_version=3;transcript_id=ENST00000263097;tag=basic;ccdsid=CCDS12053.1
19	ensembl_havana	mRNA	1026610	1039063	.	+	.	ID=transcript:ENST00000348419;Name=CNN2-003;biotype=protein_coding;version=3;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420296;havana_version=2;transcript_id=ENST00000348419;tag=basic;ccdsid=CCDS12054.1
19	havana	mRNA	1026618	1038026	.	+	.	ID=transcript:ENST00000565096;Name=CNN2-004;biotype=protein_coding;version=2;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420297;havana_version=2;transcript_id=ENST00000565096;tag=basic
19	havana	mRNA	1026628	1038377	.	+	.	ID=transcript:ENST00000562958;Name=CNN2-005;biotype=protein_coding;version=2;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420298;havana_version=1;transcript_id=ENST00000562958;tag=basic
19	havana	mRNA	1026690	1036676	.	+	.	ID=transcript:ENST00000562075;Name=CNN2-006;biotype=protein_coding;version=2;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420299;havana_version=2;transcript_id=ENST00000562075
19	havana	mRNA	1026690	1036687	.	+	.	ID=transcript:ENST00000607102;Name=CNN2-015;biotype=protein_coding;version=1;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000471312;havana_version=1;transcript_id=ENST00000607102
19	havana	mRNA	1031136	1038063	.	+	.	ID=transcript:ENST00000568865;Name=CNN2-009;biotype=protein_coding;version=1;Parent=gene:ENSG00000064666;havana_transcript=OTTHUMT00000420302;havana_version=1;transcript_id=ENST00000568865
```

> The function `db.children()` extracts all the children features of a parent. In the case above, we have iterated over all mRNAs of the gene ENSG00000064666. As you can see, this gene has seven different mRNAs. We can also print, for instance, all CDS included in the first of these mRNA (ENST00000263097) by typing:

```
>>> first_tr=db['transcript:ENST00000263097']
>>> for cds in db.children(first_tr, featuretype='CDS', order_by='start'):
...		print cds
...
19	ensembl_havana	CDS	1026661	1026723	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1031070	1031191	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1032391	1032457	.	+	1	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1032558	1032695	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1036129	1036245	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1036415	1036561	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
19	ensembl_havana	CDS	1037624	1037899	.	+	0	ID=CDS:ENSP00000263097;Parent=transcript:ENST00000263097;protein_id=ENSP00000263097
```

> Note that, in contrast with gene and transcript features, in CDS the field “phase” shows a value and attributes contain information about the protein ID.

At this point we are ready to return with our `Python` script. Remember that our objective is to obtain all the proteins encoded by the gene APOE, and in fact, we know now how to use `gffutils` to iterate over all mRNAs of a gene of interest. However, we need first the Ensembl ID of APOE to be able to access the features of this gene in the database. Of course, you can easely obtain this ID from genomic databases with a simple search in the web. Nevertheless, we can take advantage of our database at the same time that we learn another interesting `gffutils` function:

```
>>> for g in db.features_of_type('gene'):
...	name=g.attributes[‘Name’]
...	if name == ['APOE']:
...		id=g.attributes[‘ID’]
...		print id
['gene:ENSG00000130203']
```

> With the function `db.features_of_type()` you can iterate over all features of the same type in the database. 

-----------------------------------------------------------------------------

**Let’s continue here with our script**

It would be interesting to create a file named “APOE_proteins.fas” to save the proteins encoded by the APOE gene (=the output file for our script). Add this new line to the script:

```python
file = open(myGENE+'_proteins.fasta','w')
```

We can now use gffutils to iterate over all mRNAs of the APOE gene:

```python
apoe='gene:ENSG00000130203'
for mrna in db.children(apoe, featuretype='mRNA', order_by='start'):
	print >>file, '>'+apoe+';'+ mrna.attributes['ID']
```

> Note that we have added an instruction so that, for each mRNA, the script writes a line in the output file that begins with a symbol ">" and that contains the Ensembl IDs of the gene and the transcript. As you know, this line will allow the output format to be FASTA.

To obtain the proteins encoded by these mRNAs, we first need to extract their coding sequences (CDS). We must therefore introduce a nested loop to iterate over the CDS features of each transcript:

```python
	string_cds=''
	for cds in db.children(mrna, featuretype='CDS', order_by='start'):
```

At this point, the instructions we add will apply to each CDS of each mRNA of the APOE gene. First, and in order, to all the CDS of the first RNA, then to the second, and so on.

And which actions we need to obtain the complete CDS of each transcript?
First we need to extract the DNA sequence of each children CDS feature of the same transcript from the file “Homo_sapiens.GRCh37.dna.chromosome.19.fa” using the `gffutils` function `feature.sequence()`, and then connect them in the same sequence. 

```python
		seq = cds.sequence(myFASTA)
		string_cds += seq
		complete_cds = Seq(string_cds, generic_dna)
```

> The variable `string_cds` contains just a `Python` string with A, T, C, and G symbols. We used the Biopython object `Seq()` to combine the string with a biological alphabet, in this case the generic alphabet for DNA sequences.

The last line introduced in our script closes the for loop over the CDS of then same transcript. Next step would be to translate this sequence. We can use one of the many methods incorporated in `Seq()`,the method `translate()`:

```python
		protein = complete_cds.translate()
```

Finally, we print the protein in the output file, just after the FASTA formatted line with the IDs, and close the file:

```python
		for aa in range(0, len(protein), 60):
			print >>file, protein[aa:aa+60]
		file.close()
```

> Using this `for` loop, we are printing the protein sequences in lines of equal length (60 amino acids per line). 

The script is now ready. Run the script in the same directory in which you have downloaded the genomic data files `/myworkdir` in the terminal and check results. You can validate them in: 
http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000130203;r=19:44905791-44909393.

**EXAMPLE 2**

In this second example, you will see how to extract information from a VCF file. We will use as an example, the VCF of the human chromosome 19 from the 1000 genomes project phase 3; https://www.internationalgenome.org/category/phase-3/) therefore, information about variants. The specific objective is to estimate SNP diversity in the APOE gene and obtain the frequency of its two common allelic variants (rs429358 and rs7412) in African and European populations.

Of course, first we should obtain the VCF file, or at least the part of this file that corresponds to the APOE gene. The Htslib project, which is distributed with samtools, includes `tabix`, an interesting command-line utility to download ONLY the region of interest from a VCF file. To do that, `tabix` needs the specific coordinates of the region in the scaffold or chromosome. Again, we can use `gffutils` (e.g. directly in the `Python` console; still in the /myworkdir directory) to obtain these coordinates:

```python
>>> import gffutils
>>> myGFF="Homo_sapiens.GRCh37.87.chromosome.19.gff3"
>>> db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)
>>> mygene=db['gene:ENSG00000130203']
>>> mygene.start
45409011
>>> mygene.end
45412650
```

Close the `Python` console, open a new command-line terminal in your computer and type:

```bash
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 19:45409011-45412650 > APOE.vcf 
```

> The variants annotated in the region specified in the coordinates 19:45409011-45412650 are now in the file “APOE.vcf” in our computer. 

Next step is to estimate nucleotide diversity in the APOE gene region. We can obtain, for instance, the distribution of this diversity across the gene using VCFtools, another interesting ngs program that includes a flag to calculate SNP density:

```bash
vcftools --vcf APOE.vcf --SNPdensity 100 --out APOE
```

> In this particular case, we are estimating SNP density across the APOE gene in non-overlapping windows of 100 pb. The number of SNP both per 100 bp window and per kb, are saved in the file “APOE.snpden” [The program prints some warnings in the stdout about INFO entries, don’t worry about that].

A very illustrative way to visualize SNP density across the APOE gene is to plot `VCFtools` results.  One possibility is to use a plotting function in `R`(>). Open a command-line or terminal in your computer and then type `R. Use these `R` commands to generate a very simple plot:

```R
> pdf("SNPdensity.pdf")
> data<-read.table(file="APOE.snpden", header=T)
> plot(data[,3],type = "l", col="blue", ylab="SNPs", xlab="APOE gene (bp)",xaxt="n")
> xtick <- seq(0, 37, by=5)
> axis(side=1, at=xtick, labels = seq(0,3700, by=500))
> dev.off()
```

Most VCF files also contain information on the frequency and geographic distribution of genetic variants. If we want to know, for example, the frequency of the two SNPs in the APOE gene associated with Alzheimer's disease. The simplest way to do that is to create a new VCF file that includes only the information of these SNPs using `VCFtools` and extract allele frequency information using the vcf-query (a utility that is distributed with `VCFtools`):

```bash
vcftools --vcf APOE.vcf --snps SNP.txt --recode --recode-INFO-all --out SNPs_only
```

> The file “SNPs.txt” is a file with the list of SNPs (dbSNP identifiers). Note that when we create a new VCF file, we must apply the `–recode` and `--recode-INFO-all` option to write out the variants that pass-through filters and to include all data from the INFO fields in the output.
  
Using the new generated VCF file, we can extract the frequency of the SNPs in the populations of iterest:

```bash
vcf-query -f 'The frequency of %ID in African and European populations is %INFO/AFR_AF and %INFO/EUR_AF, respectively\n' SNPs_only.recode.vcf
The frequency of rs429358 in African and European populations is 0.2678 and 0.1551, respectively
The frequency of rs7412 in African and European populations is 0.1029 and 0.0626, respectively
```

In `vcf-query`, we can access to the `INFO` field (%INFO) of the VCF file and print a specific entry (in this case `AFR_AF` and `EUR_AF`; to know these specific tags you can print the header of the VCF file using the utility `bcftools`:

```bash
bcftools view -h APOE.vcf
```



# CRISPR-wheat
Make it easy to design CRISPR/Cas9 guide RNAs in wheat and other polyploid species.

## OVERVIEW

Design group specific guide RNAs (gRNAs) in wheat. You can set the targeting sequence names, so the gRNAs are only target these sequences.

There are good websites for doing the same job, but right now they do not have wheat genome. If they provided wheat genome someday, then you can ignore this tool.

Two good websites:

**CRISPR-P 2.0**: http://crispr.hzau.edu.cn/CRISPR2/

**CRISPR-GE**: http://skl.scau.edu.cn/

------

## Software dependences

- **python 2**

- **BatMis**: https://code.google.com/archive/p/batmis/ (Chandana Tennakoon, Rikky W. Purbojati, Wing-Kin Sung; BatMis: a fast algorithm for k-mismatch mapping, Bioinformatics, Volume 28, Issue 16, 15 August 2012, Pages 2122–2128, https://doi.org/10.1093/bioinformatics/bts339)

- **bedtools**: bedtools.readthedocs.io/

------

## Input

A fasta file with all the homeolog sequences. Both cDNA and genomic DNA sequences are okay, but I would recommend cDNAs.


## Parameters
``` sh
get_guide_RNA.py
	-i <sequence.fa>
	-t <seqID1,seqID2>
	-o <output file name>
	-b <blast all produced gRNAs against the genomes: 1 for YES and 0 for NO, default is NO>
	-l <sgRNA length, default is 20>
	-a <whether to target all homeologs: 1 for YES and 0 for NO. default is 1>
	-p <PAM sequence. default is NGG>
	-q <PAM position: left or right>
	-h <print the help>
	-r <reference location: BatMis processed files>
	-c <Cas9 etc cut position, for testing restriction enzymes>
	-f <gff3 file for gene locations: a bedtools sorted and filtered gff3 file which only kept "gene", "exon", "five_prime_UTR", and "three_prime_UTR".>
```
------

## Example
``` sh
cd GW2mRNAs
../bin/get_guide_RNA.py -i GW2_reference.fas -t ALL -p NGG -q right -l 20 -c 17 -o selected_gRNAs.txt -b 1
```
------

## Output

**ID**: gRNA sequence ID

**Start**: start position on the template

**End**: end position on the template

**Strand**: forward or reverse strand of the template

**Length**: length of the gRNAs

**Sequence (5' -> 3')**: gRNA 5' -> 3'sequence without PAM

**GC_content_All**: GC content of the whole gRNA sequence

**GC_content_first_10nt**: GC content of the first 10nts from the PAM of the whole gRNA sequence

**Reverse Complement (5' -> 3')**: Reverse complement sequence of the gRNAs

**PAM**: PAM sequence on the template

**Template used**: the sequence that was used to search for potential gRNAs

**Restriction Enzyme**: restriction enzyme sites that may be destroyed by the deletion.

and more for off target checking (top 10 potential targets in the genome and their gene location: exon or intron etc)

------

## Off target consideration

Based on the review paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4877526/), for Cas9:

1. In most cases, Cas9/gRNA cannot recognize a DNA site bearing more than three mismatches;

2. Cas9/gRNA cannot recognize and edit a DNA site with any number of mismatches near a PAM (within 10–12 bp);

Cpf1 has low off-target editing rates. Based on Kim et al. (http://www.ncbi.nlm.nih.gov/pubmed/27272384) and Kleinstiver et al. (http://dx.doi.org/10.1101/057802):

1. Double mismatches ablated Cpf1 activity, except when they were present in the 3’ end of the target sequence (bases 19-23)

2. Kleinstiver et al. reported that Cpf1 can tolerate mismatches at gRNA positions 1, 8, 9, and 19-23.

3. Kim et al. found that Cpf1 could tolerate single or double mismatches in the 3' PAM-distal region, but not in the 5' PAM-proximal region.

So overall, it seems that mismatches close to the PAM are not tolerable, but distal of PAM is for both Cas9 and Cpf1.

------
## Commands used to prepare BatMis input and gff3 file
I installed **BatMis** and **Bedtools** in my PATH. Please also make sure they are in your PATH, otherwise you need to update the commands in the script with the full path of the two tools.
```sh
## BatMis
# wheat genome is too big for BatMis to index, so I have to split them into chromosomes
splitfasta.pl 161010_Chinese_Spring_v1.0_pseudomolecules.fasta
# then index each chromosome: some chromosomes did not complete the whole process, but we already get the files needed later. So no worry.
cat chrlist.txt | xargs -n 1 -P 5  -I % sh -c 'build_index %' > tempout.txt

## Prepare gff3 used for the tool
# sort first
bedtools sort -i iwgsc_refseqv1.0_HighConf_2017Mar13.gff3 > sorted_HC.gff3
# filter out CDS and mRNA terms
awk '!/CDS/ && !/mRNA/' sorted_HC.gff3 > filtered_sorted_HC.gff3
```
------
## Acknowledgements

I borrowed a lot of ideas from CRISPR-P 2.0. Free software used in the package includes:

- **Muscle**: http://www.drive5.com/muscle/

- **BatMis**: https://code.google.com/archive/p/batmis/ (Chandana Tennakoon, Rikky W. Purbojati, Wing-Kin Sung; BatMis: a fast algorithm for k-mismatch mapping, Bioinformatics, Volume 28, Issue 16, 15 August 2012, Pages 2122–2128, https://doi.org/10.1093/bioinformatics/bts339)

- **bedtools**: bedtools.readthedocs.io/

- **blast+**: https://blast.ncbi.nlm.nih.gov/Blast.cgi

- **Cas9 sgRNA on target and off target scoring scripts**: John G. Doench*, Nicolo Fusi*, Meagan Sullender*, Mudra Hegde*, Emma W. Vaimberg*, Katherine F. Donovan, Ian Smith, Zuzana Tothova, Craig Wilen , Robert Orchard, Herbert W. Virgin, Jennifer Listgarten*, David E. Root. Optimized sgRNA design to maximize activity and minimize off-target effects for genetic screens with CRISPR-Cas9. Nature Biotechnology, 2016.












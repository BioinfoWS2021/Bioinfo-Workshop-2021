## Tutorial/Practical

##### *Pierre Ramond, Swan LS Sow*
##### *Thursday, 27th May 2021*

---

### 0. Background
Besides the ease, simplicity, speed and relatively lower cost of long-reads generated from third generation sequencing tech, one of the key advantages of longer amplicons lies in its potential of increased taxonomic resolution that cannot be achieved by targeted 16S/18S rRNA gene sub-region sequencing used in short-read sequencing platforms<sup>1</sup>. Longer reads also allow for significant improvements of genome assemblies<sup>2</sup>.

The primary aim of this tutorial is to illustrate the benefits of longer amplicons by comparing the quality of taxonomic assignments of long versus short amplicons of the same sequences.

<p>&nbsp;</p>

### 1. Dataset
We will work with long-read 16S and 18S rRNA gene amplicon datasets generated from samples taken from the Wadden Sea time series and Rothera Oceanographic and Biological Time Series (RaTS - Southern Ocean) in this practical session.

![WS_SO_Together](https://user-images.githubusercontent.com/25445935/119269708-5a3b0400-bbf9-11eb-8827-19c438bab66c.png)
Fig. 1: Map of sampling sites

16S amplicons were generated with the primer set A519F-1492R-pB-3771<sup>3</sup>, while the 18S amplicons were generated with the primer set Euk528F-U1391R<sup>4</sup>.

First let's get some data into the project folder. Login to ADA, and *cd* to the following directory where all files required for this tutorial are located:

```
cd /export/lv4/projects/workshop_2021/S13_LongRead/
```

Make a symbolic link of the sequence data fastq files to your working folders. You are encouraged to figure out how to do this based on unix commands covered in the previous workshop sessions, but for simplicity we have also provided the code below.

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;View code</a>    
</summary>
<pre><code>cp -rsv /export/lv4/projects/workshop_2021/S13_LongRead/ /export/lv3/scratch/bioinformatics_workshop/Users/your_username</code></pre>
<pre><code>cd /export/lv3/scratch/bioinformatics_workshop/Users/your_username/S13_LongRead/</code></pre>
</details>


<p>&nbsp;</p>

**Programs/packages required for this tutorial:**

1. cutadapt v1.16
2. seqkit v0.14.0
3. MOTHUR v.1.40.4

The programs should already installed on ADA, but if *seqkit* is not available to you, you can install it with:

```
conda install -c bioconda seqkit
```
<p>&nbsp;</p>

### 2. Extracting specific sub-regions of the 16S & 18S rRNA gene
The original reads generated from the MinION sequencing are ~1100 bp for the 16S amplicons and ~1200 bp for the 18S amplicons. We will use *cutadapt* to trim the sequences to the desired fragment lengths and extract specific 16S and 18S rRNA gene sub-regions.

First let's create a directory for the output of *cutadapt*:

```
cd /export/lv3/scratch/bioinformatics_workshop/Users/your_username/S13_LongRead/reads
mkdir cutadapt_reads
```

This is how cutadapt works:

cutadapt -j 0 -e 0.3 -O 12 \ # parameters on the errors tolerated to recognize the primers
  --discard-untrimmed \ # we want only the reads that were cut
  -a CCAGCASCYGCGGTAATTCC...TYRATCAAGAACGAAAGT \ # sequence 3-5' of forward primer ... reverse complement of the sequence 3-5' of reverse primer
  -a ACTTTCGTTCTTGATYRA...GGAATTACCGCRGSTGCTGG \ # sequence 3-5' of reverse primer ... reverse complement of the sequence 3-5' of forward primer
  -M 600 \ # to make sure we get the subregions targeted by the primers we add a read length selection
  -o cutadapt_reads/18S_sub_V4_STOECK.fasta \ # output directory
18S.fastq # file to work on

To extract the 18S V4 region, we use the primer sequences that were developed by Stoeck et.al. (2010) as the adapter sequence parameter in *cutadapt* as follows:

```
cutadapt -j 0 -e 0.3 -O 12 --discard-untrimmed -a CCAGCASCYGCGGTAATTCC...TYRATCAAGAACGAAAGT -a ACTTTCGTTCTTGATYRA...GGAATTACCGCRGSTGCTGG -M 600 -o cutadapt_reads/18S_sub_V4_STOECK.fasta 18S.fastq
```

The 16S V4 region can be extracted by using the 515F-806R primer sequences developed by Caporaso et.al. (2011):

```
cutadapt -j 0 -e 0.3 -O 12 --discard-untrimmed -a GTGCCAGCMGCCGCGGTAA...ATTAGAWADDDBDGTAGTCC -a GGACTACHVHHHTWTCTAAT...TTACCGCGGCKGCTGGCAC -M 600 -o cutadapt_reads/16S_sub_V4_806R.fasta 16S.fastq
```

Or alternatively, with the 515F-926R primer pair<sup>7</sup>, which targets the V4-V5 region:

```
cutadapt -j 0 -e 0.3 -O 12 --discard-untrimmed -a GTGYCAGCMGCCGCGGTAA...AAACTYAAAKRAATTGRCGG -a CCGYCAATTYMTTTRAGTTT...TTACCGCGGCKGCTGRCAC -M 600 -o cutadapt_reads/16S_sub_V4_926R.fasta 16S.fastq
```

Now, let's make a list of the reads that matched the adapter(primer) sequences from the *cutadapt* step above.

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;View code</a>    
</summary>
<pre><code>grep ">" cutadapt_reads/18S_sub_V4_STOECK.fasta | sed 's/>//' | sed 's/\s.*$//' > cutadapt_reads/18S_reads_ID.txt
grep ">" cutadapt_reads/16S_sub_V4_806R.fasta | sed 's/>//' | sed 's/\s.*$//' > cutadapt_reads/16S_806R_reads_ID.txt
grep ">" cutadapt_reads/16S_sub_V4_926R.fasta | sed 's/>//' | sed 's/\s.*$//' > cutadapt_reads/16S_926R_reads_ID.txt</code></pre>
</details>

<p>&nbsp;</p>

We'll then extract the long reads that came through the *cutadapt* pipeline based on the list of reads with *seqkit*'s *grep* function:

```
seqkit grep -f cutadapt_reads/18S_reads_ID.txt 18S.fastq -o cutadapt_reads/18S_og_reads.fastq
seqkit grep -f cutadapt_reads/16S_806R_reads_ID.txt 16S.fastq -o cutadapt_reads/16S_og_reads_806R.fastq
seqkit grep -f cutadapt_reads/16S_926R_reads_ID.txt 16S.fastq -o cutadapt_reads/16S_og_reads_926R.fastq
```

Finally, we'll remove the adapters, primers and Unique Molecular Identifiers (UMIs)<sup>9</sup> from the long reads by trimming the first and last 80 bp of each sequence with *seqkit*'s 'subseq' function. Assigning taxonomy to reads with adapters or UMIs would probably lead to spurious taxonomic annotation.

```
cd cutadapt_reads
seqkit subseq -r 80:-80  18S_og_reads.fastq > 18S_og_reads_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_806R.fastq > 16S_og_reads_806R_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_926R.fastq > 16S_og_reads_926R_trimm.fastq
```
<p>&nbsp;</p>

### 3. Generating taxonomy of the V4 sub-region fragments

First, make a directory for the files to annotate, i.e. the read files for the sub-regions (_sub_V4_) and the original trimmed reads (_og_*_trimm.fastq)

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;View code</a>    
</summary>
<pre><code>mkdir sub_regions

cp 16S_og_reads_806R_trimm.fastq sub_regions
cp 16S_og_reads_926R_trimm.fastq sub_regions
cp 16S_sub_V4_806R.fasta sub_regions
cp 16S_sub_V4_926R.fasta sub_regions
cp 18S_og_reads_trimm.fastq sub_regions
cp 18S_sub_V4_STOECK.fasta sub_regions</code></pre>
</details>

<p>&nbsp;</p>

We'll use *MOTHUR* to classify our sequences with the silva.nr_v138 database for the 16S sequences and pr2_version_4.13.0_18S for the 18S sequences. *MOTHUR* works better with sequences in the fasta format, so we'll first convert all fastq sequences to fasta format with *seqtk*:

```
cd sub_regions

seqtk seq -A 16S_og_reads_926R_trimm.fastq > `ls 16S_og_reads_926R_trimm.fastq | sed 's/\.fastq/\.fasta/'`
seqtk seq -A 16S_og_reads_806R_trimm.fastq > `ls 16S_og_reads_806R_trimm.fastq | sed 's/\.fastq/\.fasta/'`
seqtk seq -A 18S_og_reads_trimm.fastq > `ls 18S_og_reads_trimm.fastq | sed 's/\.fastq/\.fasta/'`
```

Classify the sequences with *MOTHUR* using the default *wang* method<sup>8</sup>. This should take ~ 5-10 minutes. The cutoff value indicates that only classified sequences with bootstrap values ≥80% will be returned.

The mothur command is composed of the directory (set.dir) of the database of annotation used: here we used mothur-formated versions of silva.nr_v138 (contain species annotation up to the genus level only) and pr2_version_4.13.0_18S (up to species level).

The file to classify: our fasta files

The cutoff for bootstrap values ≥80%: taxonomic levels under this thresholds are considered unclassified for all reads (you wonder what is a bootstrap value really? https://www.researchgate.net/post/What-is-a-bootstrap-value#58ef7ba6b0366d8cd460abfe).

Dont' forget to change *username* with your username!

```
# 16S files
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/16S_og_reads_806R_trimm.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/16S_og_reads_926R_trimm.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/16S_sub_V4_806R.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/16S_sub_V4_926R.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"
```
```
# 18S files
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/18S_og_reads_trimm.fasta, reference=pr2_version_4.13.0_18S_mothur.fasta,taxonomy=pr2_version_4.13.0_18S_mothur.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv3/scratch/workshop_2021/Users/*username*/S13_LongRead/reads/cutadapt_reads/sub_regions/18S_sub_V4_STOECK.fasta, reference=pr2_version_4.13.0_18S_mothur.fasta,taxonomy=pr2_version_4.13.0_18S_mothur.tax, cutoff=80)"
```
<p>&nbsp;</p>

### 4. Generating length gradient fragments from the trimmed long reads

To further test the theory, we'll generate sequence fragments from the original long reads with 100bp length variations (i.e. 100bp, 200bp, 300bp... 1000bp, original). Make a new folder for the length gradient fragments, and then subfolders within that folder for the 16S and 18S fragments. Make a copy of the fastq file containing the trimmed long reads.

```
cd /export/lv3/scratch/bioinformatics_workshop/Users/*username*/S13_LongRead/reads/cutadapt_reads/
mkdir Length_gradients
mkdir Length_gradients/18S
mkdir Length_gradients/16S
```

We copy and rename the trimmed orginal reads (we use only the subset of reads that matched the 515F-806R primers for 16S)

```
cp 18S_og_reads_trimm.fastq Length_gradients/18S/18S_trim_original.fastq
cp 16S_og_reads_806R_trimm.fastq Length_gradients/16S/16S_trim_original.fastq
```

<p>&nbsp;</p>

Then use *cutadapt* to trim the fragment to different lengths:

```
# 16S
cd Length_gradients/16S/
cutadapt -l 100 -o 16S_trim_100bp.fastq 16S_trim_original.fastq
cutadapt -l 200 -o 16S_trim_200bp.fastq 16S_trim_original.fastq
cutadapt -l 300 -o 16S_trim_300bp.fastq 16S_trim_original.fastq
cutadapt -l 400 -o 16S_trim_400bp.fastq 16S_trim_original.fastq
cutadapt -l 500 -o 16S_trim_500bp.fastq 16S_trim_original.fastq
cutadapt -l 600 -o 16S_trim_600bp.fastq 16S_trim_original.fastq
cutadapt -l 700 -o 16S_trim_700bp.fastq 16S_trim_original.fastq
cutadapt -l 800 -o 16S_trim_800bp.fastq 16S_trim_original.fastq
cutadapt -l 900 -o 16S_trim_900bp.fastq 16S_trim_original.fastq
cutadapt -l 1000 -o 16S_trim_1000bp.fastq 16S_trim_original.fastq
```

```
# 18S
cd Length_gradients/18S/
cutadapt -l 100 -o 18S_trim_100bp.fastq 18S_trim_original.fastq
cutadapt -l 200 -o 18S_trim_200bp.fastq 18S_trim_original.fastq
cutadapt -l 300 -o 18S_trim_300bp.fastq 18S_trim_original.fastq
cutadapt -l 400 -o 18S_trim_400bp.fastq 18S_trim_original.fastq
cutadapt -l 500 -o 18S_trim_500bp.fastq 18S_trim_original.fastq
cutadapt -l 600 -o 18S_trim_600bp.fastq 18S_trim_original.fastq
cutadapt -l 700 -o 18S_trim_700bp.fastq 18S_trim_original.fastq
cutadapt -l 800 -o 18S_trim_800bp.fastq 18S_trim_original.fastq
cutadapt -l 900 -o 18S_trim_900bp.fastq 18S_trim_original.fastq
cutadapt -l 1000 -o 18S_trim_1000bp.fastq 18S_trim_original.fastq
```

We'll add the fragment size to the sequence header so we can easily identify the different size fragments of the same sequence:

```
# 16S:
sed -i '1~4s/\s\+/_16S_trim_100bp.fastq /' 16S_trim_100bp.fastq
sed -i '1~4s/\s\+/_16S_trim_200bp.fastq /' 16S_trim_200bp.fastq
sed -i '1~4s/\s\+/_16S_trim_300bp.fastq /' 16S_trim_300bp.fastq
sed -i '1~4s/\s\+/_16S_trim_400bp.fastq /' 16S_trim_400bp.fastq
sed -i '1~4s/\s\+/_16S_trim_500bp.fastq /' 16S_trim_500bp.fastq
sed -i '1~4s/\s\+/_16S_trim_600bp.fastq /' 16S_trim_600bp.fastq
sed -i '1~4s/\s\+/_16S_trim_700bp.fastq /' 16S_trim_700bp.fastq
sed -i '1~4s/\s\+/_16S_trim_800bp.fastq /' 16S_trim_800bp.fastq
sed -i '1~4s/\s\+/_16S_trim_900bp.fastq /' 16S_trim_900bp.fastq
sed -i '1~4s/\s\+/_16S_trim_1000bp.fastq /' 16S_trim_1000bp.fastq
sed -i '1~4s/\s\+/_original /' 16S_trim_original.fastq
```

```
# 18S:
sed -i '1~4s/\s\+/_18S_trim_100bp.fastq /' 18S_trim_100bp.fastq
sed -i '1~4s/\s\+/_18S_trim_200bp.fastq /' 18S_trim_200bp.fastq
sed -i '1~4s/\s\+/_18S_trim_300bp.fastq /' 18S_trim_300bp.fastq
sed -i '1~4s/\s\+/_18S_trim_400bp.fastq /' 18S_trim_400bp.fastq
sed -i '1~4s/\s\+/_18S_trim_500bp.fastq /' 18S_trim_500bp.fastq
sed -i '1~4s/\s\+/_18S_trim_600bp.fastq /' 18S_trim_600bp.fastq
sed -i '1~4s/\s\+/_18S_trim_700bp.fastq /' 18S_trim_700bp.fastq
sed -i '1~4s/\s\+/_18S_trim_800bp.fastq /' 18S_trim_800bp.fastq
sed -i '1~4s/\s\+/_18S_trim_900bp.fastq /' 18S_trim_900bp.fastq
sed -i '1~4s/\s\+/_18S_trim_1000bp.fastq /' 18S_trim_1000bp.fastq
sed -i '1~4s/\s\+/_original /' 18S_trim_original.fastq
```

Again, we'll need to convert all the the fastq files to fasta before the taxonomic annotation step. We'll use a loop to repeat the conversion for all the fastq files in the folder. (remember, change *username*).

```
for i in /export/lv3/scratch/bioinformatics_workshop/Users/*username*/S13_LongRead/reads/cutadapt_reads/Length_gradients/*/*fastq;
do
seqtk seq -A $i > `ls $i | sed 's/\.fastq/\.fasta/'`;
done
```

Finally, let's assign taxonomy for all sequences of the different length gradients:

```
# Assign taxonomy for all fasta files 16S folder (I swear to god, change *username* or I will track you, find you and kill you)

for i in /export/lv3/scratch/bioinformatics_workshop/Users/*username*/S13_LongRead/reads/cutadapt_reads/Length_gradients/16S/*.fasta;
do
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=$i, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"
done
```

```
# Assign taxonomy for all fasta files 18S folder (you know what to do by now... or you are dead)

for i in /export/lv3/scratch/bioinformatics_workshop/Users/*username*/S13_LongRead/reads/cutadapt_reads/Length_gradients/18S/*.fasta;
do
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=$i, reference=pr2_version_4.13.0_18S_mothur.fasta, taxonomy=pr2_version_4.13.0_18S_mothur.tax, cutoff=80)"
done
```

<p>&nbsp;</p>

### 5. Compare taxonomic assignment quality

We'll now graphically compare the effects of long versus short amplicons of the same sequences on the quality of taxonomic assignments. We'll do this in R studio, with an R script that can be found in:

```
/export/lv3/scratch/bioinformatics_workshop/Users/*username*/S13_LongRead/scripts
```

To launch R studio in ADA: [http://ada.nioz.nl:8787/](http://ada.nioz.nl:8787/)
<p>&nbsp;</p>

![16S_plot](https://user-images.githubusercontent.com/25445935/119827287-40e7cf80-bef9-11eb-81e1-eb27bedf7fc5.png)
Fig. 2

<p>&nbsp;</p>

![18S_plot](https://user-images.githubusercontent.com/25445935/119827412-62e15200-bef9-11eb-851a-4c45f9ee7a48.png)
Fig. 3

---
###### Refs:

<sup>1. Johnson, J.S., Spakowicz, D.J., Hong, BY. et al. Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis. Nat Commun 10, 5029 (2019). https://doi.org/10.1038/s41467-019-13036-1</sup>

<sup>2. Koren, S., and Phillippy, A.M. (2015). One chromosome, one contig: complete microbial genomes from long-read sequencing and assembly. Current Opinion in Microbiology 23, 110-120.</sup>

<sup>3. Martijn, J., Lind, A.E., Schon, M.E., et al. (2019). Confident phylogenetic identification of uncultured prokaryotes through long read amplicon sequencing of the 16S-ITS-23S rRNA operon. Environ Microbiol 21, 2485-2498.

<sup>4. Edgcomb, V., Orsi, W., Bunge, J., et al. (2011) Protistan microbial observatory in the Cariaco Basin, Caribbean. I. Pyrosequencing vs Sanger insights into species richness. ISME J 5, 1344–1356. https://doi.org/10.1038/ismej.2011.6</sup>

<sup>5. Stoeck, T., Bass, D., Nebel, M., et al. (2010), Multiple marker parallel tag environmental DNA sequencing reveals a highly complex eukaryotic community in marine anoxic water. Molecular Ecology, 19: 21-31. https://doi.org/10.1111/j.1365-294X.2009.04480.x </sup>

<sup>6. Caporaso, J.G., Lauber, C.L., Walters, W.A., et al. (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. Proceedings of the National Academy of Sciences 108, 4516.</sup>

<sup>7. Parada, A.E., Needham, D.M. and Fuhrman, J.A. (2016), Primers for marine microbiome studies. Environ Microbiol, 18: 1403-1414. https://doi.org/10.1111/1462-2920.13023</sup>

<sup>8. Wang, Q., Garrity, G.M., Tiedje, J.M., et al. (2007). Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology 73, 5261-5267.</sup>

<sup>9. Hiatt, J.B., Patwardhan, R.P., Turner, E.H., et.al. (2010). Parallel, tag-directed assembly of locally derived short sequence reads. Nat Methods 7, 119-122.</sup>

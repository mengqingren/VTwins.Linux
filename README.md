## VTwins.Linux

**This standalone R script can perform phenotype-enriched features including species or pathways for metagenomic sequencing or 16S sequencing.**

### Installation
<pre>
git clone https://github.com/mengqingren/VTwins.Linux.git
#or download directly
unzip VTwins.Linux-main.gz
</pre>

### Description
- PairFind.Linux.R : R script to perform phynotype-enriched features
- Shuffle.WilcoxonScore.R : R script to process data and plot
- test.data : Directory including mock datasets and real datas to test the R scripts
- Unique.Pair.Permutation: Processed data and figures
- Species.EightDatasets : Directory of CRC integrated datasets for species
- Pathway.EightDatasets : Directory of CRC integrated datasets for pathways


### Dependence

**R package** 
- R version 4.0.5
- tidyverse_1.3.1
- optparse_1.6.6

### Basic Usage
<pre>
Rscript PairFind.Linux.R -d test.data/test.data.txt -g test.data/test.phenodata.txt -a test.data/ -m euclidean -s 10000 -u 0.8 -n 0.2 -w ShuffleWstat -b BoundarySample -p BoundaryPair -o Results.txt -c 0.01
</pre>

### Parameter Description
- `-d`: must be a data frame with columns representing features and rows representing samples.
- `-g`: must be a data frame with two columns, and colnames are `id` (column 1) and `grp` (column 2). column id represent the sample id, column grp consist of `grp1` and `grp2`, representing the ctrl and disease, repsectively.
- `-m`: distance calculating method. it must bu consist with the method in `dist` function.
- `-a`: filename of output directory. Default: ./
- `-w`: filename of shuffle W stats
- `-b`: filename of output boundary samples with distance
- `-p`: filename of output final pairs 
- `-s`: shuffle times
- `-n`: lower percentage of shuffle pairs
- `-u`: higher percentage of shuffle pairs
- `-c`: p value cutoff of `Incre.aveRank.P` and `Decre.aveRank.P`

**Note:** if the sample pairs are less than 10, it will return nothing.

### Output
- **Decre.aveRank.P** : P value based on averange rank. Generally, to confirm the significant features, we usually use variable `Decre.aveRank.P` to evaluate the disease-enriched features and `Incre.aveRank.P` to evaluate the control enriched features. 
- **Incre.aveRank.P** : P value based on averange rank. Generally, to confirm the significant features, we usually use variable `Decre.aveRank.P` to evaluate the disease-enriched features and `Incre.aveRank.P` to evaluate the control enriched features. 
- Decre.minRank. : P value based on min rank.
- Incre.minRank.P : P value based on min rank.
- Decre.maxRank.P : P value based on max rank.
- Incre.maxRank.P : P value based on max rank.
- **Species** : Feature
- **Enriched** : phenotype-enriched groups with provided p value cutoff
- Decre.maxRank.P.FDR : P.adjust of Decre.maxRank.P
- Decre.minRank.P.FDR : P.adjust of Decre.minRank.P
- **Decre.aveRank.P.FDR** : P.adjust of Decre.aveRank.P
- **Incre.aveRank.P.FDR** : P.adjust of Incre.aveRank.P
- Decreasing.Rank.Max : shuffle W stats rank based on max method by decreasing
- Increasing.Rank.Max : shuffle W stats rank based on max method by increasing
- Decreasing.Rank.Min : shuffle W stats rank based on min method by decreasing
- Increasing.Rank.Min : shuffle W stats rank based on min method by increasing
- **Decreasing.Rank.Average** : shuffle W stats rank based on average method by decreasing
- **Increasing.Rank.Average** : shuffle W stats rank based on average method by increasing
- **Ctlmean**: mean value of relativa abundance for paired control samples
- **Dismean**: mean value of relativa abundance for paired disease samples

**Important** : Must keep the `Same Order` of samples in  relative abundance dataframe and phenotype dataframe. And we also keep the samples `clustered and placed` according the `grp1 and grp2` of variable `grp` in phenotype data. 

**You can refer to the examples data format showing in [mengqingren/VTwins](https://github.com/mengqingren/VTwins)**

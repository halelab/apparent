# apparent

### Description
The 'apparent' R package performs parentage analysis based on a test of genetic identity between expected progeny (EPij), built using homozygous loci from all pairs of possible parents (i and j), and all potential offspring (POk). Using the Gower Dissimilarity metric (GD) (Gower 1971), genetic identity between EPij and POk is taken as evidence that individuals i and j are the true parents of offspring k.  Evaluation of triad (two parents + offspring) significance is based on the distribution of all GD (EPij|k) values. Specifically, a Dixon test (Dixon 1950; 1951) is used to identify a gap-based threshold that separates true triads and from spurious associations.  
For any offspring not successfully assigned to a pair of parents, perhaps due to the absence of one parent from the test population, a non-mandatory Dyad analysis can be employed to identify a likely single parent for a given offspring. In this analysis, a two-stage test is applied to discriminate an offspring's true parent from its other close relatives (e.g. siblings) that may also be present in the population. In the first stage, 'apparent' calculates the mean GD (GDM) between a POk and all expected progeny arising from the j possible triads involving potential parent i. In the second stage, it calculates a coefficient of variation (GDCV) among the pairwise GD's between POk and each expected progeny arising from the j triads involving potential parent i. An individual that is simultaneously a low outlier in the first test and a high outlier in the second is identified as a likely parent of POk. In an effort to facilitate interpretation, results of both the triad and optional dyad analyses are presented in tabular and graphical form.

**User note**: As the size of the exploratory triad space inflates (i.e. as the population size grows), the likelihood of low-probability instances of low GD values, occurring purely by chance, increases.  Such outlier values can "bleed" into the gap between true and spurious triads, in effect making the gap disappear and leading to a null result, even if true triads exist in the population (Type II error).  Such a dilution of a true gap is generally due to the presence of one or more problematic genotypes in the population, problematic due to the fact that the available number of loci for their analysis is much for some reason lower than for other individuals, leading to artificially low GDij|k values.  To address such problems the following two "best practices" are recommended to users:
1.  Always include at least one known (true) triad in the analysis, as a reference/control.
2.  After running the triad analysis, inspect the summary file "apparent-Triad-Summary.txt" to see if any individuals exhibit unusually low mean GDij|k values or mean number of loci usable for analysis.  If such individuals exist, remove them and re-run the analysis.

### Usage
- Download the [apparent_TestData.txt](https://github.com/halelab/apparent/blob/master/apparent_TestData.txt) for your R working directory.  
- Download the [apparent R script](https://github.com/halelab/apparent/blob/master/apparent.R) and load it into R.
```
In R
# Install/Load the required "outliers" R package
install.packages("outliers")
library(outliers)
# Set the path for your working directory
setwd("/path/to/work/dir")
# Load the apparent_TestData
InputFile <- read.table(file="apparent_TestData.txt",sep="\t",h=F)
# Run apparent
apparent (InputFile, MaxIdent = 0.1, alpha = 0.01, nloci = 100, FullList = TRUE, self = TRUE, plot = TRUE, Dyad = FALSE)
```

### Arguments
- **InputFile**: A tab-delimited text file with n rows (n = number of individuals in the population) and m+2 columns (m = number of SNP loci). Column one contains the ID's of the individuals in the population. Column two contains a classification key which assigns each individual to one of five possible classes from analysis (see below). The third and all subsequent columns contain genotype calls, with one SNP per column and the alleles separated by "/". A missing genotype is represented by "-/-".  An example input file with 5 individuals and 5 SNP loci is shown below:  

|Genotype|key|loci01|loci02|loci03|loci04|loci05|
|---|---|---|---|---|---|---|
|Individual01|All|A/A|A/T|C/A|C/C|T/C|
|Individual02|Pa|T/A|A/A|-/-|C/C|T/T|
|Individual03|Mo|A/A|A/A|C/A|C/G|T/T|
|Individual04|Fa|T/T|A/T|A/A|G/G|-/-|
|Individual05|Off|A/A|T/T|C/C|C/C|T/C|

The keys values allowed in the second column are:  
*All*: "All" - The individual will be tested as a potential Mother, Father, and Offspring;  
*Pa*: "Parent" - The individual will be tested as a potential Mother and Father;  
*Mo*: "Mother" - The individual will be tested only as a potential Mother;  
*Fa*: "Father" - The individual will be tested only as a potential Father;  
*Off*: "Offspring" - The individual will be tested only as a potential Offspring.  

- **MaxIdent**: Sets the maximum triad GDij|k to be considered for outlier significance testing. This parameter directly impacts computation time. By default, MaxIdent is set to 0.1.  
- **alpha**: The alpha level for all significance testing (triad and dyad analyses).  
- **nloci**: The minimum acceptable number loci to be used when computing the pairwise GDij|k. The default value of 100 is suggested, based on previous investigations. All triads for which the number of usable SNPs falls below nloci will be excluded from the analysis.  
- **FullList**: Logical value for creating an output file containing all pairwise comparisons from the triad analysis. The default value is TRUE.  
- **self**: Logical value for instructing 'apparent' whether or not to consider self-crossing (parent i = parent j). The default value is TRUE.  
- **plot**: Logical value for plotting the results of both the triad and dyad analyses. The default value is TRUE.   
- **Dyad**: Logical value for instructing 'apparent' to perform a dyad analysis, following the triad analysis. The default value is FALSE.  

### Outputs
The pedigree information for the 12 progenies (four families of three full-sibs each) from controlled crosses used to assess the performance of apparent can be found [here](https://github.com/halelab/apparent/blob/master/Families.xlsx).  

- [**apparent-Triad-All.txt**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Triad-All.txt) contains the full table of results from the triad analysis, including the cross type (self or out-cross), the number of usable SNPs, and the GDij|k, of all tested triads.  
- [**apparent-Triad-Summary.txt**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Triad-Summary.txt) reports useful summary statistics from the triad analysis.  For the population as a whole: Overall mean GDij|k and its standard deviation, as well as overall mean number of usable SNPs and its standard deviation. For each individual genotype in the analysis: Mean GDij|k, GDij|k range, and mean usable number of SNPs for all comparisons involving that genotype.  Genotypes exhibiting a mean GDij|k or number of usable SNPs less than 2 SD's below the population means are flagged with an asterisk.   
- [**apparent-Triad-Sig.txt**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Triad-Sig.txt) reports only those triads found to be significant.  
- [**apparent-Triad-Plot.png**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Triad-Plot.png) is a plot of the distribution of GDij|k values, annotated with the gap-based threshold that separates true triads from spurious associations.  The plot is useful in interpretating the results of the triad analysis, the full details of which are found in **apparent-Triad-All.txt** and **apparent-Triad-Sig.txt**.  
- [**apparent-Dyad-Sig.txt**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Dyad-Sig.txt) reports only those parent-offspring dyads found to be significant.  
- [**apparent-Dyad-Plot.pdf**](https://github.com/halelab/apparent/blob/master/outputs/apparent-Dyad-Plot.pdf) is a two-paneled figure showing the distributions of GDM and GDCV values upon which the dyad analysis is based. These plots are useful in interpretating the results of the dyad analysis, the full details of which are found in **apparent-Dyad-Sig.txt**.

### Reference
Gower JC. A general coefficient of similarity and some of its properties. Biometrics. 1971;27:857-871.  
Dixon WJ. Analysis of extreme values. Ann. Math. Stat. 1950;21(4):488-506.  
Dixon WJ. Ratios involving extreme values. Ann. Math. Stat. 1951;22(1):68-78.


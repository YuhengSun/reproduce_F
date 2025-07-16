# reproduce_F
This is my workflow to reproduce Francesco's analysis.

## PCA
Load plink module for analysis:

```
module load PLINK/1.9b_6.13-x86_64
```

#### Linkage pruning
I started with ```VCF=wholegenome_sparrows_variants_norm.vcf.gz```. According to Francesco, it's a variant only .vcf file of all the concatenated chromosomes (incl. sex chromosome Z) and scaffolds. It was normalized by the command ```bcftools norm```.

```
plink --vcf $VCF \
  --double-id \
  --allow-extra-chr \
  --chr-set 28 \
  --set-missing-var-ids @:# \
  --indep-pairwise 10 10 0.1 \
  --out sparrows_ys
```

The parameters ```10 10 0.1``` were found in Francesco's thesis: “All loci within a 10kb window that showed linkage above 0.1 were removed”.

This produced ```sparrows_ys.prune.in``` and ```sparrows_ys.prune.out```.

#### Perform PCA

```
plink --vcf $VCF \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --extract sparrows_ys.prune.in \
  --make-bed \
  --pca \
  --out sparrows_ys \
  --mind 0.99 \
  --chr-set 28
```

```--mind 0.99``` is from Francesco's note.

This produced ```sparrows_ys.bed```, ```sparrows_ys.bim```, ```sparrows_ys.fam```, ```sparrows_ys.eigenval``` and ```sparrows_ys.eigenvec```.
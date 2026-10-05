parallel -j 32 'bedtools shuffle -i non_genic_circles.bed -g hg38.genome -excl ENCFF356LFX.bed -chrom > random/permutation_{}.bed' ::: {1..1000}

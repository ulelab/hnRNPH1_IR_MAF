zcat PSI-Metadata.tab.gz | awk -F'\t' 'NR > 1 {
  split($3, a, "[:-]");
  strand = substr($5, length($5), 1);  # last character of FullCO
  print a[1] "\t" a[2]-1 "\t" a[3] "\t" $1 "\t0\t" strand;
}' > introns.bed

cut -f1,2 /scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38.p14-GencodeRelease44/GRCh38.primary_assembly.genome.fa.fai > /scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38.p14-GencodeRelease44/hg38.genome

#bedtools flank -i introns.bed -g /scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38.p14-GencodeRelease44/hg38.genome -l 50 -r 50 -s > introns_flank50.bed 

awk '{
  start = ($2 - 50 >= 0) ? $2 - 50 : 0;  # avoid negative start
  end = $3 + 50;
  print $1 "\t" start "\t" end "\t" $4 "\t" $5 "\t" $6;
}' introns.bed > introns_flank50.bed

bedtools getfasta \
  -fi /scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38.p14-GencodeRelease44/GRCh38.primary_assembly.genome.fa \
  -bed introns_flank50.bed \
  -s -name \
  -fo sequences.fa

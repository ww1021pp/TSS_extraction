#####################local computer##########
###
cd ~/Document/mm10/
awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$7,$9}}' gencode.vM25.primary_assembly.annotation.gff3 | grep "gene_name" | sed 's/ID=//'|sed 's/;/\t/'| sed 's/gene_id.*gene_name=//g'|sed 's/;.*//g' > GenCode.vM25.mm10.Gene_bed/Gencode.vM25.gene.bed
cd ./GenCode.vM25.mm10.Gene_bed
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$NF"\t"$4}' Gencode.vM25.gene.bed | sort -k1,1 -k2,2n > genCode.vM25.gene.sorted.bed
scp -r GenCode.vM25.mm10.Gene_bed pliu@10.14.18.147:/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/
cd /workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/

####awk forward and revers #######
awk '($6 == "+") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 2500), ($2 + 1000), $4, $5, $6  }' > genCode.vM25.gene.tssUp2.5kbDn1kb.for.padded.bed
awk '($6 == "-") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($3 - 1000), ($3 + 2500), $4, $5, $6  }' > genCode.vM25.gene.tssUp2.5kbDn1kb.tss.rev.padded.bed

bedops --everything genCode.vM25.gene.tssUp2.5kbDn1kb.for.padded.bed genCode.vM25.gene.tssUp2.5kbDn1kb.tss.rev.padded.bed > GenCodeVM25.gene.tssup2.5kbDn1kb.padded.bed
bedtools intersect -a GenCodeVM25.gene.tssup2.5kbDn1kb.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.tssup2.5kbDn1kb.padded.filtered.bed


######################genebody +- 100 kb#########
 cat genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 100000), ($3 + 100000), $4, $5, $6  }' |sed 's/-[0-9]\+/0/' >genCode.vM25.gene.geneBodyUPDn100kb.padded.bed
 bedtools intersect -a genCode.vM25.gene.geneBodyUPDn100kb.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.geneBodyUPDn100kb.padded.filtered.bed


 ##########TSS 5kb, 10kb, 20kb###

 awk '($6 == "+") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 5000), ($2 + 5000), $4, $5, $6  }' > genCode.vM25.gene.tss5K.for.padded.bed
 awk '($6 == "-") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($3 - 5000), ($3 + 5000), $4, $5, $6  }' > genCode.vM25.gene.tss5K.rev.padded.bed
 bedops --everything genCode.vM25.gene.tss5K.for.padded.bed genCode.vM25.gene.tss5K.rev.padded.bed > GenCodeVM25.gene.tss5K.padded.bed
 bedtools intersect -a GenCodeVM25.gene.tss5K.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.tss5K.filtered.bed

 awk '($6 == "+") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 2000), ($2 + 2000), $4, $5, $6  }' > genCode.vM25.gene.tss2K.for.padded.bed
 awk '($6 == "-") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($3 - 2000), ($3 + 2000), $4, $5, $6  }' > genCode.vM25.gene.tss2K.rev.padded.bed
 bedops --everything genCode.vM25.gene.tss2K.for.padded.bed genCode.vM25.gene.tss5K.rev.padded.bed > GenCodeVM25.gene.tss2K.padded.bed
 bedtools intersect -a GenCodeVM25.gene.tss2K.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.tss2K.filtered.bed


 awk '($6 == "+") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 10000), ($2 + 10000), $4, $5, $6  }' > genCode.vM25.gene.tss10K.for.padded.bed
 awk '($6 == "-") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($3 - 10000), ($3 + 10000), $4, $5, $6  }' > genCode.vM25.gene.tss10K.rev.padded.bed
 bedops --everything genCode.vM25.gene.tss10K.for.padded.bed genCode.vM25.gene.tss10K.rev.padded.bed > GenCodeVM25.gene.tss10K.padded.bed
 bedtools intersect -a GenCodeVM25.gene.tss10K.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.tss10K.filtered.bed


 awk '($6 == "+") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($2 - 20000), ($2 + 20000), $4, $5, $6  }' > genCode.vM25.gene.tss20K.for.padded.bed
 awk '($6 == "-") { print $0 }' genCode.vM25.gene.sorted.bed | awk 'BEGIN{ OFS="\t" }{ print $1, ($3 - 20000), ($3 + 20000), $4, $5, $6  }' > genCode.vM25.gene.tss20K.rev.padded.bed
 bedops --everything genCode.vM25.gene.tss20K.for.padded.bed genCode.vM25.gene.tss20K.rev.padded.bed > GenCodeVM25.gene.tss20K.padded.bed
 bedtools intersect -a GenCodeVM25.gene.tss20K.padded.bed -b ../mm10.bound.bed > GenCodeVM25.gene.tss20K.filtered.bed

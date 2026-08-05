Code availability
All analysis steps were performed on a Linux system, using R version 4.3.1 and Python version 3.12.9. The codes and software parameters used in the pipeline are as follows:
(1)Post-Sequencing Data Quality Control
fastp Version: 0.23.4
Usage:
fastp \
        -i ${fq1} -I ${fq2} \
        -o ${sample_id}_clean_1.fq.gz -O ${sample_id}_clean_2.fq.gz \
        --cut_tail --cut_tail_window_size 1 \
        --cut_tail_mean_quality 20 --average_qual 20 \
        --length_required 50 --thread 4 \
(2)Mapping
Bwa Version: 0.7.18-r1243-dirty
Samblaster Version: 0.1.26
Samtools Version:1.17
Usage:
        bwa mem -k 18 -T 25 -t 8 -M ref_genome.fa \
        ${sample_id}_clean_1.fq.gz ${sample_id}_clean_2.fq.gz \
        | samblaster -M \
        | samtools sort \
        -O bam -T  ${sample_id}_temp -@ 4 > ${sample_id}.bam && \
        samtools index ${sample_id}.bam
(3)Variant Detection
Gatk Version:gatk-4.5.0.0
Usage:
gatk --java-options "-Xmx15G" HaplotypeCaller -R $fa -ERC GVCF -I ${sample_id}.sorted.rmdup.bam -O ${sample_id}.sorted.rmdup.bam.gvcf.gz --tmp-dir tmp --native-pair-hmm-threads 6
gatk --java-options "-Xmx15G" GenomicsDBImport \
        -R $fa --variant ${sample1}.g.vcf.gz --variant ${sample2}.g.vcf.gz \
        --genomicsdb-workspace-path genomicsdb --tmp-dir tmp --merge-input-intervals true
gatk --java-options "-Xmx15G" GenotypeGVCFs \
        -R $fa -V gendb://genomicsdb -O cohort.raw.vcf.gz --tmp-dir tmp
gatk SelectVariants -R $fa -V cohort.raw.vcf.gz --select-type-to-include SNP -O cohort.raw.snp.vcf.gz
gatk SelectVariants -R $fa -V cohort.raw.vcf.gz --select-type-to-include INDEL -O cohort.raw.indel.vcf.gz
gatk VariantFiltration -R $fa -V cohort.raw.snp.vcf.gz \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        -O hardfilted.snp.vcf.gz
gatk VariantFiltration -R $fa -V cohort.raw.indel.vcf.gz \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "FS > 200.0" --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --filter-expression "SOR > 10.0" --filter-name "SOR10" \
        -O hardfilted.indel.vcf.gz
bcftools view -f PASS hardfilted.snp.vcf.gz | bgzip > hardfilted.snp.pass.vcf.gz
bcftools view -f PASS hardfilted.indel.vcf.gz | bgzip > hardfilted.indel.pass.vcf.gz
(4)Variant Filtering
Perl version: v5.26.2
Usage:
Perl mis_filter_vcf_V7.pl -vcf hardfilted.snp.pass.vcf.gz -out ./ -gq_all 20 -dp_min1 3 -dp_max1 999999999999999 -DP_min 50 -DP_max 999999999999999 -gq_single 0 -miss 0.1 -maf 0.05 -num 240 -dp_pos 2 -gq_pos 3 -name snp.filter-dp3-miss0.1-maf0.05.nodouble
Perl mis_filter_vcf_V7.pl -vcf hardfilted.indel.pass.vcf.gz -out ./ -gq_all 20 -dp_min1 3 -dp_max1 999999999999999 -DP_min 50 -DP_max 999999999999999 -gq_single 0 -miss 0.3 -maf 0.05 -num 240 -dp_pos 2 -gq_pos 3 -name indel.filter-dp3-miss0.3-maf0.05.nodouble
(5)Variant Annotation
ANNOVAR Version: 2020-06-07
Usage:
gzip -dc snp.filter-dp3-miss0.1-maf0.05.nodouble.vcf.gz > out.filted.snp.vcf
convert2annovar.pl -format vcf4 -snpqual 0 out.filted.snp.vcf > out.snp.avinput
annotate_variation.pl -buildver out -geneanno out.snp.avinput ${annoLIBRARY}
statSNPAnno.pl out.snp.avinput.variant_function out.snp.avinput.exonic_variant_function out.filted.snp.vcf out.SNP_Annotation_Statistics.xls

gzip -dc indel.filter-dp3-miss0.3-maf0.05.nodouble.vcf.gz > out.filted.indel.vcf
convert2annovar.pl -format vcf4 -snpqual 0 out.filted.indel.vcf > out.indel.avinput
annotate_variation.pl -buildver out -geneanno out.indel.avinput ${annoLIBRARY}
statindelAnno.pl -outDir ${proj_dir}/anno/INDEL_anno -projDir ${proj_dir}/anno/INDEL_anno -gsize ${genomesize}

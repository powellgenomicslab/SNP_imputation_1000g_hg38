
#!/usr/bin/env python
shell.executable('bash')



# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pgen",
        psam = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.psam",
        pvar = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pvar"
    output:
        bed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmap_input.bed",
        outbed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmap_output.bed",
        excluded_ids = output_dict["output_dir"] + "/crossmapped/{ancestry}_excluded_ids.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"]
    threads:
        imputation_dict["crossmap_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        in_plink = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset",
        out = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink",
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        singularity exec --bind {params.bind} {params.sif} CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        singularity exec --bind {params.bind} {params.sif} awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.in_plink} --exclude {output.excluded_ids} --make-bed --output-chr MT --out {params.out}
        singularity exec --bind {params.bind} {params.sif} awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        pgen = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bed",
        psam = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bim",
        pvar = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"]
    threads:
        imputation_dict["sort_bed_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink",
        out = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.infile} --make-bed --max-alleles 2 --output-chr MT --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi"
    output:
        bed = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bed",
        bim = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bim",
        fam = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_memory"],
        java_mem = lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_java_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_memory"]
    threads:
        imputation_dict["harmonize_hg38_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}",
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -jar {params.jar}\
            --input {params.infile}\
            --inputType PLINK_BED\
            --ref {input.vcf}\
            --refType VCF\
            --update-id\
            --output {params.out}
        """


rule plink_to_vcf:
    input:
        bed = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bed",
        bim = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bim",
        fam = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.fam"
    output:
        data_vcf_gz = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"]
    threads:
        imputation_dict["plink_to_vcf_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}",
        out = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.infile} --recode vcf id-paste=iid --chr 1-22 --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.data_vcf_gz}
        """


rule vcf_fixref_hg38:
    input:
        fasta = fasta,
        vcf = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi",
        data_vcf = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf = output_dict["output_dir"] + "/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = output_dict["output_dir"] + "/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"]
    threads:
        imputation_dict["vcf_fixref_hg38_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} | \
        singularity exec --bind {params.bind} {params.sif} bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}

        #Index
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = output_dict["output_dir"] + "/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"]
    threads:
        imputation_dict["filter_preimpute_vcf_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        maf = lambda wildcards: float(maf_df["MAF"][maf_df.Ancestry == wildcards.ancestry].values),
        missing = imputation_dict["snp_missing_pct"],
        hwe = imputation_dict["snp_hwe"]
    shell:
        """
        #Add tags
        singularity exec --bind {params.bind} {params.sif} bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        singularity exec --bind {params.bind} {params.sif} bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} |\
        singularity exec --bind {params.bind} {params.sif} bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        singularity exec --bind {params.bind} {params.sif} bcftools filter -e "ALT='.'" |\
        singularity exec --bind {params.bind} {params.sif} bcftools norm -d all |\
        singularity exec --bind {params.bind} {params.sif} bcftools norm -m+any |\
        singularity exec --bind {params.bind} {params.sif} bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        #Index the output file
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.filtered_vcf}
        """

rule het:
    input:
        vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf = temp(output_dict["output_dir"] + "/het/{ancestry}_filtered_temp.vcf"),
        inds = output_dict["output_dir"] + "/het/{ancestry}_het_failed.inds",
        het = output_dict["output_dir"] + "/het/{ancestry}_het.het",
        passed = output_dict["output_dir"] + "/het/{ancestry}_het_passed.inds",
        passed_list = output_dict["output_dir"] + "/het/{ancestry}_het_passed.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"]
    threads: imputation_dict["het_threads"]
    params:
        het_base = output_dict["output_dir"] + "/het/{ancestry}_het",
        script = "/opt/SNP_imputation_1000g_hg38/Imputation/scripts/filter_het.R",
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het/{ancestry}_het",
        sif = input_dict["singularity_image"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}
        singularity exec --bind {params.bind} {params.sif} vcftools --vcf {output.tmp_vcf} --het --out {params.het_base}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list = output_dict["output_dir"] + "/het/{ancestry}_het_passed.txt",
        vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz",
        index = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"]
    threads: imputation_dict["het_filter_threads"]
    params:
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}

        #Index the output file
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule calculate_missingness:
    input:
        filtered_vcf = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals = output_dict["output_dir"] + "/genotype_donor_annotation/{ancestry}_individuals.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"]
    threads:
        imputation_dict["calculate_missingness_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_genotypes"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.filtered_vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}

        singularity exec --bind {params.bind} {params.sif} vcftools --gzvcf {output.tmp_vcf} --missing-indv --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bcftools query -l {input.filtered_vcf} >> {output.individuals}
        """


rule split_by_chr:
    input:
        filtered_vcf = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = output_dict["output_dir"] + "/split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        index = output_dict["output_dir"] + "/split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"]
    threads:
        imputation_dict["split_by_chr_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = output_dict["output_dir"] + "/split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = genetic_map,
        phasing_file = phasing_dir + "/chr{chr}.bcf"
    output:
        vcf = output_dict["output_dir"] + "/eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"]
    threads: imputation_dict["eagle_prephasing_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/eagle_prephasing/{ancestry}_chr{chr}_phased"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        """


rule minimac_imputation:
    input:
        vcf = output_dict["output_dir"] + "/eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file = impute_dir + "/chr{chr}.m3vcf.gz"
    output:
        output_dict["output_dir"] + "/minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"]
    threads:
        imputation_dict["minimac_imputation_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/minimac_imputed/{ancestry}_chr{chr}",
        minimac4 = "/opt/bin/minimac4",
        chunk_length = imputation_dict["chunk_length"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.minimac4} --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix {params.out} \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads} \
            --ChunkLengthMb {params.chunk_length}
        """


rule combine_vcfs_ancestry:
    input:
        vcfs = lambda wildcards: expand(output_dict["output_dir"] + "/minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", chr = chromosomes, ancestry = ancestry_subsets)
    output:
        combined = output_dict["output_dir"] + "/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        ind = output_dict["output_dir"] + "/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_ancestry_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_ancestry_memory"]
    threads: imputation_dict["combine_vcfs_ancestry_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        files_begin = output_dict["output_dir"] + "/minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools concat -Oz {params.files_begin} > {output.combined}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.combined}
        """


rule combine_vcfs_all:
    input:
        vcfs = lambda wildcards: expand(output_dict["output_dir"] + "/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets)
    output:
        combined = output_dict["output_dir"] + "/vcf_all_merged/imputed_hg38.vcf.gz",
        ind = output_dict["output_dir"] + "/vcf_all_merged/imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_memory"]
    threads: imputation_dict["combine_vcfs_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
    shell:
        """
        if [[ $(ls -l {input.vcfs} | wc -l) > 1 ]]
        then
            singularity exec --bind {params.bind} {params.sif} bcftools merge -Oz {input.vcfs} > {output.combined}
        else
            cp {input.vcfs} {output.combined}
        fi
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.combined}
        """

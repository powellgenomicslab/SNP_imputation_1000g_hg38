#!/usr/local/envs/py36/bin python3

rule bed2pgen:
    input:
        output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand-updated-chr23.vcf"
    output:
        output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pgen"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2pgen_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2pgen_memory"]
    threads: final_vcf_dict["bed2pgen_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        bfile = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand-updated-chr{chr}",
        out = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.bfile} --make-pgen --out {params.out}
        """

rule fix_pvar:
    input:
        pgen = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pgen",
        pvar = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pvar"
    output:
        old = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.old",
        temp = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.temp",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["fix_pvar_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["fix_pvar_memory"]
    threads: final_vcf_dict["fix_pvar_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        pfile = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cp {params.pfile}.pvar {output.old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.old} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$3,$4,$5)}}' > {output.temp}
        singularity exec --bind {params.bind} {params.sif} grep "#" {input.pvar} | sed 's/, /,/g' > {params.pfile}.pvar
        singularity exec --bind {params.bind} {params.sif} awk -F"\t" 'BEGIN{{OFS=FS = "\t"}} NR==FNR{{a[$1 FS $2 FS $3] = $0; next}} {{ind = $1 FS $2 FS $3}} ind in a {{print a[ind], $6, $7}}' {output.temp} {input.pvar} >> {params.pfile}.pvar
        """


rule pgen2vcf:
    input:
        output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.temp"
    output:
        # vcf = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf",
        gvcf = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["pgen2vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["pgen2vcf_memory"]
    threads: final_vcf_dict["pgen2vcf_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        pfile = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}",
        out = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.pfile} --recode vcf --out {params.out}
        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} tabix -p vcf {output.gvcf}
        """

if os.path.exists(output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv"):
    ancestry_file = pd.read_csv(output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv", sep = "\t", header = None, names = ["Ancestry"])

rule combine_ancestries:
    input:
        vcfs = lambda wildcards: expand(output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf.gz", chr = wildcards.chr, ancestry = ancestry_file["Ancestry"])
    output:
        output_dict["output_dir"] + "/vcf/combined/QC_filtered_chr{chr}.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_ancestries_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_ancestries_memory"]
    threads: final_vcf_dict["combine_ancestries_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        bfile = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} vcf-merge {input.vcfs} > {output}
        """


rule vcf_sort:
    input:
        output_dict["output_dir"] + "/vcf/combined/QC_filtered_chr{chr}.vcf"
    output:
        output_dict["output_dir"] + "/vcf/combined_sorted/QC_filtered_sorted_chr{chr}.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["vcf_sort_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["vcf_sort_memory"]
    threads: final_vcf_dict["vcf_sort_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools sort {input} -Oz -o {output}
        """

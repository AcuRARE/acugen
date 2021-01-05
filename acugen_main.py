import sys, os, shutil
import docker
import json
import gzip
from pathlib import Path
import subprocess

def main(fa_read1, fa_read2, out_dir='/workspace/acugen/out/'):

    #DEFAULT DIRECTORIES

    working_dir = '/workspace/acugen/'
    ref_file_dir = '/workspace/GRCh38/GRCh38.p13.genome.fa'
    dbsnp_file_dir = '/workspace/dbSNP/00-All.vcf.gz'
    docker_version_list = '/workspace/pipe/config/version_list.json'

    output_dir = os.path.abspath(out_dir)
    if not os.path.isdir(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    os.chdir(working_dir)
    fa_read1_dir = fa_read1
    fa_read2_dir = fa_read2

    read1_bname = os.path.basename(fa_read1_dir)
    read2_bname = os.path.basename(fa_read2_dir)

    client = docker.from_env()
    with open(docker_version_list) as handle:
        ver = json.loads(handle.read())

    bwa_thread_count = 20
    picard_xmx_val = 64 # in Gig's
    gatk_xmx_val = 64

    def get_tool_key(tool_name):
        # Will return if N > 1 tool with same suffix
        return([key + ':' + val for key, val in ver.items() if key.endswith(tool_name)])

    #READ-GROUP DATA

    LB = 'NULL' #Read-Group library
    PL = 'ILLUMINA' #Read-Group platorm
    with gzip. open(fa_read1_dir, 'r') as f:
        RGID = '_'.join(f.readline().decode('ascii').split(':')[0:4])[1:] # Read Group ID / Platform Unit
    PU = RGID + '.' + LB
    fa_read1_bname = os.path.basename(fa_read1_dir)
    SM = fa_read1_bname[0:fa_read1_bname.index('R1')-1] #sample ID
    RG_str = f'@RG\\tID:{RGID}\\tSM:{SM}\\tPL:{PL}\\tLB:{LB}\\tPU:{PU}'

    # Temp file dirs
    base_temp_dir = working_dir + f'data/temp/{SM}/'
    Path(base_temp_dir).mkdir(parents=True, exist_ok=True)

    bwa_sam_dir = base_temp_dir + 'aln.sam'
    samtools_sorted_bam_dir = base_temp_dir + 'aln.sorted.bam'
    samtools_sorted_bai_dir = base_temp_dir + 'aln.sorted.bam.csi'
    picard_dup_bam_dir = base_temp_dir + 'aln.dup.bam'
    picard_metrics_dir = base_temp_dir + 'aln_metrics.txt'
    gatk_recal_table_dir = base_temp_dir + 'aln_recal.table'
    gatk_recal_bam_dir = base_temp_dir + 'aln_recal.bam'
    gatk_recal_bai_dir = base_temp_dir + 'aln_recal.bam.csi'
    samtools_mpileup_dir = base_temp_dir + 'pre-call_mp.vcf'

    Path(os.path.join(output_dir, f'individual_vcfs/{SM}/')).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(output_dir, 'merged_vcfs/')).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(output_dir, 'bam/')).mkdir(parents=True, exist_ok=True)

    #Output file dirs
    samtools_bam_dir = os.path.join(output_dir, f'bam/{SM}.bam')
    gatk_hc_vcf_dir = os.path.join(output_dir, f'individual_vcfs/{SM}/{SM}_hc.vcf')
    platypus_vcf_dir = os.path.join(output_dir, f'individual_vcfs/{SM}/{SM}_pl.vcf')
    bcf_mp_vcf_dir = os.path.join(output_dir, f'individual_vcfs/{SM}/{SM}_mp.vcf')
    freebayes_vcf_dir = os.path.join(output_dir, f'individual_vcfs/{SM}/{SM}_fb.vcf')
    survivor_merged_vcf_dir = os.path.join(output_dir, f'merged_vcfs/{SM}_srv_merged.vcf')
    bcftools_merged_vcf_dir = os.path.join(output_dir, f'merged_vcfs/{SM}_btls_merged.vcf')

    indiv_vcfs_tree = os.path.join(output_dir, f'individual_vcfs/{SM}')


    for tool in ver:
        print(tool + ':' + ver[tool], end=' ---')
        try:
            container = client.containers.run(tool + ':' + ver[tool])
        except:
            print(' ERR')
            print(container.logs())
        else:
            print(' DONE')

    print(f'Sample ID: {SM}.')

    #BWA

    bwa_key = get_tool_key('bwa')
    bwa_sam_bname = os.path.basename(bwa_sam_dir)
    bwa_sam_tree = os.path.split(bwa_sam_dir)[0:-1][0]
    ref_file_bname = os.path.basename(ref_file_dir)
    ref_file_tree = os.path.split(ref_file_dir)[0:-1][0]

    Path(bwa_sam_dir).touch()
    print('1. Running "bwa mem"')
    client.containers.run(*bwa_key, f"bash -c 'bwa mem -M -t {bwa_thread_count} -R \"{RG_str}\" ./ref/{ref_file_bname} ./{read1_bname} ./{read2_bname} > ./out/{bwa_sam_bname}'", working_dir='/data/',  volumes={fa_read1_dir:{'bind':f'/data/{read1_bname}', 'mode':'rw'}, fa_read2_dir:{'bind':f'/data/{read2_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref/', 'mode':'rw'}, bwa_sam_dir:{'bind':f'/data/out/{bwa_sam_bname}', 'mode':'rw'}})

    #SAMTOOLS

    sam_bcftools_key = get_tool_key('samtools_bcftools')
    samtools_bam_bname = os.path.basename(samtools_bam_dir)
    samtools_bam_tree = os.path.split(samtools_bam_dir)[0:-1][0]
    samtools_sorted_bam_bname = os.path.basename(samtools_sorted_bam_dir)
    samtools_sorted_bam_tree = os.path.split(samtools_sorted_bam_dir)[0:-1][0]
    samtools_sorted_bai_bname = os.path.basename(samtools_sorted_bai_dir)
    samtools_sorted_bai_tree = os.path.split(samtools_sorted_bai_dir)[0:-1][0]

    Path(samtools_bam_dir).touch()
    print('2. Running "samtools view"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools view -bSh {bwa_sam_bname} > {samtools_bam_bname}"', working_dir='/data/',  volumes={bwa_sam_dir:{'bind':f'/data/{bwa_sam_bname}', 'mode':'rw'}, samtools_bam_dir:{'bind':f'/data/{samtools_bam_bname}', 'mode':'rw'}})

    Path(samtools_sorted_bam_dir).touch()
    print('3. Running "samtools sort"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools sort {samtools_bam_bname} > {samtools_sorted_bam_bname}"', working_dir='/data/',\
        volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, samtools_bam_dir:{'bind':f'/data/{samtools_bam_bname}', 'mode':'rw'}})

    Path(samtools_sorted_bai_dir).touch()
    print('4. Running "samtools index"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools index -bc {samtools_sorted_bam_bname}"', working_dir='/data/',\
        volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, samtools_sorted_bai_dir:{'bind':f'/data/{samtools_sorted_bai_bname}', 'mode':'rw'}})

    # PICARD

    picard_key = get_tool_key('picard')
    picard_dup_bam_bname = os.path.basename(picard_dup_bam_dir)
    picard_dup_bam_tree = os.path.split(picard_dup_bam_dir)[0:-1][0]
    picard_metrics_bname = os.path.basename(picard_metrics_dir)
    picard_metrics_tree = os.path.split(picard_metrics_dir)[0:-1][0]

    Path(picard_dup_bam_dir).touch()
    Path(picard_metrics_dir).touch()
    print('5. Running "picard MarkDuplicates"')
    client.containers.run(*picard_key, f'picard MarkDuplicates I={samtools_sorted_bam_bname} O={picard_dup_bam_bname} M={picard_metrics_bname} \
        REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE', working_dir='/data/', volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, \
        picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, picard_metrics_dir:{'bind':f'/data/{picard_metrics_bname}', 'mode':'rw'}})

    # GATK

    gatk_key = get_tool_key('gatk')
    gatk_recal_table_bname = os.path.basename(gatk_recal_table_dir)
    gatk_recal_table_tree = os.path.split(gatk_recal_table_dir)[0:-1][0]
    gatk_recal_bam_bname = os.path.basename(gatk_recal_bam_dir)
    gatk_recal_bam_tree = os.path.split(gatk_recal_bam_dir)[0:-1][0]
    dbsnp_file_bname = os.path.basename(dbsnp_file_dir)
    dbsnp_file_tree = os.path.split(dbsnp_file_dir)[0:-1][0]
    gatk_hc_vcf_bname = os.path.basename(gatk_hc_vcf_dir)
    gatk_recal_bai_bname = os.path.basename(gatk_recal_bai_dir)

    Path(gatk_recal_table_dir).touch()
    print('6. Running "gatk BaseRecalibrator"')
    client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" BaseRecalibrator -I {picard_dup_bam_bname} --known-sites ./dbsnp/{dbsnp_file_bname} -R ./ref/{ref_file_bname} -O {gatk_recal_table_bname}', working_dir='/data/', volumes={dbsnp_file_tree:{'bind':f'/data/dbsnp/', 'mode':'rw'}, \
        picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, gatk_recal_table_dir:{'bind':f'/data/{gatk_recal_table_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}})

    Path(gatk_recal_bam_dir).touch()
    print('7. Running "gatk ApplyBQSR"')
    client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" ApplyBQSR -I {picard_dup_bam_bname} -R ./ref/{ref_file_bname} --bqsr-recal-file {gatk_recal_table_bname} -O {gatk_recal_bam_bname}', working_dir='/data/', volumes={\
        picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, gatk_recal_table_dir:{'bind':f'/data/{gatk_recal_table_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}})

    Path(gatk_recal_bai_dir).touch()
    print('8. Running "samtools index" (2)')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools index -bc {gatk_recal_bam_bname}"', working_dir='/data/',\
        volumes={gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}})

    Path(gatk_hc_vcf_dir).touch()
    print('9. Running "gatk HapplotypeCaller"')
    client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" HaplotypeCaller -I {gatk_recal_bam_bname} -R ./ref/{ref_file_bname} --dbsnp ./dbsnp/{dbsnp_file_bname} -O ./out/{gatk_hc_vcf_bname}', working_dir='/data/', volumes={dbsnp_file_tree:{'bind':f'/data/dbsnp/', 'mode':'rw'}, \
        ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}, gatk_hc_vcf_dir:{'bind':f'/data/out/{gatk_hc_vcf_bname}', 'mode':'rw'}})

    #PLATYPUS

    platypus_key = get_tool_key('platypus-variant')
    platypus_vcf_bname = os.path.basename(platypus_vcf_dir)

    Path(platypus_vcf_dir).touch()
    print('10. Running "platypus callVariants"')
    client.containers.run(*platypus_key, f'platypus callVariants --bamFiles={gatk_recal_bam_bname} --refFile=./ref/{ref_file_bname} --output=./out/{platypus_vcf_bname}', working_dir='/data/', volumes={\
        ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}, platypus_vcf_dir:{'bind':f'/data/out/{platypus_vcf_bname}', 'mode':'rw'}})

    #MPILEUP

    bcftools_key = get_tool_key('bcftools')
    bcf_mp_vcf_bname = os.path.basename(bcf_mp_vcf_dir)
    #samtools_mpileup_bname = os.path.basename(samtools_mpileup_dir)
    #samtools_mpileup_tree = os.path.split(samtools_mpileup_dir)[0:-1][0]

    Path(bcf_mp_vcf_dir).touch()
    print('11. Running "samtools mpileup | bcftools call"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools mpileup -uf ./ref/{ref_file_bname} {gatk_recal_bam_bname} | bcftools call -mv --threads 20 > ./out/{bcf_mp_vcf_bname}"', working_dir='/data/',  volumes={ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, bcf_mp_vcf_dir:{'bind':f'/data/out/{bcf_mp_vcf_bname}', 'mode':'rw'}})

    #FREEBAYES

    freebayes_key = get_tool_key('freebayes')
    freebayes_vcf_bname = os.path.basename(freebayes_vcf_dir)

    Path(freebayes_vcf_dir).touch()
    print('12. Running "freebayes"')
    client.containers.run(*freebayes_key, f'bash -c "freebayes -f ./ref/{ref_file_bname} {gatk_recal_bam_bname} > ./out/{freebayes_vcf_bname}"', working_dir='/data/', volumes={\
        ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, freebayes_vcf_dir:{'bind':f'/data/out/{freebayes_vcf_bname}', 'mode':'rw'}})

    #VCF MERGE - SURVIVOR

    survivor_key = get_tool_key('survivor')
    survivor_merged_vcf_bname = os.path.basename(survivor_merged_vcf_dir)

    Path(survivor_merged_vcf_dir).touch()
    print('13. Running "survivor"')
    client.containers.run(*survivor_key, f'bash -c "ls *vcf > sample_files; SURVIVOR merge sample_files 100 1 1 1 0 1 /out/{survivor_merged_vcf_bname}"', working_dir='/data/', volumes={\
        indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'}, survivor_merged_vcf_dir:{'bind':f'/out/{survivor_merged_vcf_bname}', 'mode':'rw'}})

    #VCF MERGE - VCFTOOLS
    #vcftools_key = get_tool_key('vcftools')

    print('14. Running "bgzip & tabix"')

    bcftools_merged_vcf_bname = os.path.basename(bcftools_merged_vcf_dir)


    indiv_vcfs_dir = os.listdir(os.path.join(output_dir, f'/individual_vcfs/{SM}/'))
    indiv_filtered_vcfs_dir = [os.path.splitext(vcf)[0] + '.filter' + os.path.splitext(vcf)[1] for vcf in indiv_vcfs_dir]
    #indiv_bgzipped_tree = os.path.join(indiv_vcfs_tree, 'bgzipped/')
    #Path(indiv_bgzipped_tree).mkdir(parents=True, exist_ok=True)
    for vcf in indiv_vcfs_dir:
        if vcf == SM+'_pl.vcf':
            continue
        if os.path.isfile(os.path.join(indiv_vcfs_tree, vcf)) and os.path.splitext(vcf)[1] == '.vcf':
            vcf_filtered_bname = os.path.splitext(vcf)[0] + '.filter' + os.path.splitext(vcf)[1]
            #vcf_filtered_dir = os.path.join(indiv_bgzipped_tree, vcf_filtered_bname)
            vcf_filtered_dir = os.path.join(indiv_vcfs_tree, vcf_filtered_bname)
            #print(vcf_filtered_dir)
            Path(vcf_filtered_dir).touch()
            client.containers.run(*bcftools_key, f'bcftools view -i \"%FILTER=\'PASS\' | %FILTER=\'.\'\" -o /out/{vcf_filtered_bname} {vcf}', working_dir='/data/', volumes={\
                indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'},\
                vcf_filtered_dir:{'bind':f'/out/{vcf_filtered_bname}', 'mode':'rw'}})
            os.remove(os.path.join(f'./out/individual_vcfs/{SM}/', vcf))
            subprocess.run(['bgzip', vcf_filtered_bname], cwd=indiv_vcfs_tree)

    for gz in os.listdir(indiv_vcfs_tree):
        if os.path.splitext(gz)[1] == '.gz':
            subprocess.run(['tabix', '-p', 'vcf', gz], cwd=indiv_vcfs_tree)

    Path(bcftools_merged_vcf_dir).touch()
    print('15. Running "bcftools merge"')

    freebayes_filtered_vcf_gz_bname = os.path.splitext(freebayes_vcf_bname)[0] + '.filter' + os.path.splitext(freebayes_vcf_bname)[1] + '.gz'
    bcf_mp_filtered_vcf_gz_bname = os.path.splitext(bcf_mp_vcf_bname)[0] + '.filter' + os.path.splitext(bcf_mp_vcf_bname)[1] + '.gz'
    gatk_hc_filtered_vcf_gz_bname = os.path.splitext(bcf_mp_vcf_bname)[0] + '.filter' + os.path.splitext(bcf_mp_vcf_bname)[1] + '.gz'

    client.containers.run(*bcftools_key, f'bash -c "bcftools merge --force-samples {freebayes_filtered_vcf_gz_bname} {bcf_mp_filtered_vcf_gz_bname} {gatk_hc_filtered_vcf_gz_bname} > /out/{bcftools_merged_vcf_bname}"', working_dir='/data/', volumes={\
        #freebayes_vcf_dir:{'bind':f'/data/{freebayes_vcf_bname}', 'mode':'rw'},\
        #bcf_mp_vcf_dir:{'bind':f'/data/{bcf_mp_vcf_bname}', 'mode':'rw'},\
        #platypus_vcf_dir:{'bind':f'/data/{platypus_vcf_bname}', 'mode':'rw'},\
        #gatk_hc_vcf_dir:{'bind':f'/data/{gatk_hc_vcf_bname}', 'mode':'rw'},\
        indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'},\
        bcftools_merged_vcf_dir:{'bind':f'/out/{bcftools_merged_vcf_bname}', 'mode':'rw'}})

    print('14. Removing temporary files in ' + base_temp_dir)
    shutil.rmtree(base_temp_dir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

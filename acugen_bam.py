import sys, os, shutil, glob
import docker
import json
import gzip
from pathlib import Path
import subprocess
import pysam
from datetime import datetime as dt

def ts():
    print('[', end='')
    print(dt.now(), end='')
    print('] - ', end='')

def main(bam_file_dir, out_dir):

    #DEFAULT DIRECTORIES

    working_dir = '/workspace/acugen/'
    docker_version_list = '/workspace/acugen/config/version_list.json'
    references_dir = '/workspace/References/'
    dbsnp_file_dir = '/workspace/References/dbSNP/00-All.vcf.gz'

    output_dir = os.path.abspath(out_dir)
    if not os.path.isdir(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    os.chdir(working_dir)
    bam_file = pysam.AlignmentFile(bam_file_dir, 'rb')
    bam_file_bname = os.path.basename(bam_file_dir)

    client = docker.from_env()
    with open(docker_version_list) as handle:
        ver = json.loads(handle.read())

    bwa_thread_count = 20
    samtools_thread_count = 20
    picard_xmx_val = 64 # in Gig's
    gatk_xmx_val = 64

    def get_tool_key(tool_name):
        # Will return if N > 1 tool with same suffix
        return([key + ':' + val for key, val in ver.items() if key.endswith(tool_name)])

    #READ-GROUP DATA

    LB = bam_file.header['RG'][0]['LB'] #Read-Group library
    PL = bam_file.header['RG'][0]['ID'] #Read-Group platorm
    RGID = bam_file.header['RG'][0]['ID']
    PU = RGID + '.' + LB
    SM = bam_file.header['RG'][0]['SM'] #sample ID
    RG_str = f'@RG\\tID:{RGID}\\tSM:{SM}\\tPL:{PL}\\tLB:{LB}\\tPU:{PU}'

    # Temp file dirs
    base_temp_dir = working_dir + f'data/temp/{SM}/'
    Path(base_temp_dir).mkdir(parents=True, exist_ok=True)

    #bwa_sam_dir = os.path.join(base_temp_dir, f'{SM}.sam')
    #samtools_bam_dir = os.path.join(base_temp_dir, f'{SM}.bam')
    samtools_sorted_bam_dir = os.path.join(base_temp_dir, f'{SM}_sorted.bam')
    samtools_sorted_bai_dir = os.path.join(base_temp_dir, f'{SM}_sorted.bam.bai')
    picard_dup_bam_dir = os.path.join(base_temp_dir, f'{SM}_dup.bam')
    picard_metrics_dir = os.path.join(base_temp_dir, f'{SM}_metrics.txt')
    gatk_recal_table_dir = os.path.join(base_temp_dir, f'{SM}_recal.table')
    gatk_recal_bam_dir = os.path.join(base_temp_dir, f'{SM}_recal.bam')
    gatk_recal_bai_dir = os.path.join(base_temp_dir, f'{SM}_recal.bam.bai')

    Path(os.path.join(output_dir, f'individual_vcfs/{SM}/')).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(output_dir, 'merged_vcfs/')).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(output_dir, f'fastqc/{SM}/')).mkdir(parents=True, exist_ok=True)
    #Path(os.path.join(output_dir, 'bam/')).mkdir(parents=True, exist_ok=True)

    #Output file dirs
    samtools_bam_dir = bam_file_dir
    fastqc_dir = os.path.join(output_dir, f'fastqc/{SM}')
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

    #FASTQC

    Path(fastqc_dir).touch()
    ts(); print('1. Running "fastqc"')
    fastqc_key = get_tool_key('fastqc')
    client.containers.run(*fastqc_key, f"bash -c 'fastqc -t 20 -f bam /input/{bam_file_bname} -o ./'", working_dir='/data/',  volumes={bam_file_dir:{'bind':f'/input/{bam_file_bname}', 'mode':'rw'}, fastqc_dir:{'bind':'/data/', 'mode':'rw'}, })

    #BWA

    # bwa_key = get_tool_key('bwa')
    # bwa_sam_bname = os.path.basename(bwa_sam_dir)
    # bwa_sam_tree = os.path.split(bwa_sam_dir)[0:-1][0]

    # Path(bwa_sam_dir).touch()
    #ts(); print('1. Skipping "bwa mem"')
    # client.containers.run(*bwa_key, f"bash -c 'bwa mem -M -t {bwa_thread_count} -R \"{RG_str}\" ./ref/{ref_file_bname} ./{read1_bname} ./{read2_bname} > ./out/{bwa_sam_bname}'", working_dir='/data/',  volumes={fa_read1_dir:{'bind':f'/data/{read1_bname}', 'mode':'rw'}, fa_read2_dir:{'bind':f'/data/{read2_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref/', 'mode':'rw'}, bwa_sam_dir:{'bind':f'/data/out/{bwa_sam_bname}', 'mode':'rw'}})

    #SAMTOOLS

    sam_bcftools_key = get_tool_key('samtools_bcftools')
    samtools_bam_bname = os.path.basename(samtools_bam_dir)
    samtools_bam_tree = os.path.split(samtools_bam_dir)[0:-1][0]
    samtools_sorted_bam_bname = os.path.basename(samtools_sorted_bam_dir)
    samtools_sorted_bam_tree = os.path.split(samtools_sorted_bam_dir)[0:-1][0]
    samtools_sorted_bai_bname = os.path.basename(samtools_sorted_bai_dir)
    samtools_sorted_bai_tree = os.path.split(samtools_sorted_bai_dir)[0:-1][0]

    #TODO: check if index (.idx) exists
    #TODO: check if sorted
    #TODO: check if index up to date

    # Path(samtools_bam_dir).touch()
    #ts(); print('2. Skipping "samtools view"')
    # client.containers.run(*sam_bcftools_key, f'bash -c "samtools view -bSh {bwa_sam_bname} > {samtools_bam_bname}"', working_dir='/data/',  volumes={bwa_sam_dir:{'bind':f'/data/{bwa_sam_bname}', 'mode':'rw'}, samtools_bam_dir:{'bind':f'/data/{samtools_bam_bname}', 'mode':'rw'}})

    Path(samtools_sorted_bam_dir).touch()
    ts(); print('2. Running "samtools sort"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools sort -@ {samtools_thread_count} {samtools_bam_bname} > {samtools_sorted_bam_bname}"', working_dir='/data/',\
        volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, samtools_bam_dir:{'bind':f'/data/{samtools_bam_bname}', 'mode':'rw'}})

    Path(samtools_sorted_bai_dir).touch()
    ts(); print('3. Running "samtools index"')
    client.containers.run(*sam_bcftools_key, f'bash -c "samtools index -@ {samtools_thread_count} -b {samtools_sorted_bam_bname}"', working_dir='/data/',\
        volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, samtools_sorted_bai_dir:{'bind':f'/data/{samtools_sorted_bai_bname}', 'mode':'rw'}})

    # PICARD

    picard_key = get_tool_key('picard')
    picard_dup_bam_bname = os.path.basename(picard_dup_bam_dir)
    picard_dup_bam_tree = os.path.split(picard_dup_bam_dir)[0:-1][0]
    picard_metrics_bname = os.path.basename(picard_metrics_dir)
    picard_metrics_tree = os.path.split(picard_metrics_dir)[0:-1][0]

    Path(picard_dup_bam_dir).touch()
    Path(picard_metrics_dir).touch()
    ts(); print('4. Running "picard MarkDuplicates"')
    client.containers.run(*picard_key, f'picard MarkDuplicates I={samtools_sorted_bam_bname} O={picard_dup_bam_bname} M={picard_metrics_bname} \
        REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE', working_dir='/data/', volumes={samtools_sorted_bam_dir:{'bind':f'/data/{samtools_sorted_bam_bname}', 'mode':'rw'}, \
        picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, picard_metrics_dir:{'bind':f'/data/{picard_metrics_bname}', 'mode':'rw'}})

    #GATK variable names
    gatk_recal_table_bname = os.path.basename(gatk_recal_table_dir)
    #gatk_recal_table_tree = os.path.split(gatk_recal_table_dir)[0:-1][0]
    gatk_recal_bam_bname = os.path.basename(gatk_recal_bam_dir)
    gatk_recal_bam_tree = os.path.split(gatk_recal_bam_dir)[0:-1][0]
    dbsnp_file_bname = os.path.basename(dbsnp_file_dir)
    dbsnp_file_tree = os.path.split(dbsnp_file_dir)[0:-1][0]
    gatk_hc_vcf_bname = os.path.basename(gatk_hc_vcf_dir)
    gatk_recal_bai_bname = os.path.basename(gatk_recal_bai_dir)
    #platypus variable names
    platypus_vcf_bname = os.path.basename(platypus_vcf_dir)
    #mpileup variable names
    bcf_mp_vcf_bname = os.path.basename(bcf_mp_vcf_dir)
    #samtools_mpileup_bname = os.path.basename(samtools_mpileup_dir)
    #samtools_mpileup_tree = os.path.split(samtools_mpileup_dir)[0:-1][0]
    #freebayes variable names
    freebayes_vcf_bname = os.path.basename(freebayes_vcf_dir)

    def annotate(ref_file_bname, ref_file_tree):

        print('Using reference file: ' + ref_file_bname)

        # GATK

        gatk_key = get_tool_key('gatk')

        Path(gatk_recal_table_dir).touch()
        ts(); print('5. Running "gatk BaseRecalibrator"')
        client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" BaseRecalibrator -I {picard_dup_bam_bname} --known-sites ./dbsnp/{dbsnp_file_bname} -R ./ref/{ref_file_bname} -O {gatk_recal_table_bname}', working_dir='/data/', volumes={dbsnp_file_tree:{'bind':f'/data/dbsnp/', 'mode':'rw'}, \
            picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, gatk_recal_table_dir:{'bind':f'/data/{gatk_recal_table_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}})

        Path(gatk_recal_bam_dir).touch()
        ts(); print('6. Running "gatk ApplyBQSR"')
        client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" ApplyBQSR -I {picard_dup_bam_bname} -R ./ref/{ref_file_bname} --bqsr-recal-file {gatk_recal_table_bname} -O {gatk_recal_bam_bname}', working_dir='/data/', volumes={\
            picard_dup_bam_dir:{'bind':f'/data/{picard_dup_bam_bname}', 'mode':'rw'}, gatk_recal_table_dir:{'bind':f'/data/{gatk_recal_table_bname}', 'mode':'rw'}, ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}})

        Path(gatk_recal_bai_dir).touch()
        ts(); print('7. Running "samtools index" (2)')
        client.containers.run(*sam_bcftools_key, f'bash -c "samtools index -@ {samtools_thread_count} -b {gatk_recal_bam_bname}"', working_dir='/data/',\
            volumes={gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}})

        Path(gatk_hc_vcf_dir).touch()
        ts(); print('8. Running "gatk HaplotypeCaller"')
        client.containers.run(*gatk_key, f'gatk --java-options "-Xmx{gatk_xmx_val}G" HaplotypeCaller -I {gatk_recal_bam_bname} -R /ref/{ref_file_bname} --dbsnp /dbsnp/{dbsnp_file_bname} -O /out/{gatk_hc_vcf_bname}', working_dir='/data/', volumes={dbsnp_file_tree:{'bind':f'/dbsnp/', 'mode':'rw'}, \
            ref_file_tree:{'bind':f'/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}, gatk_hc_vcf_dir:{'bind':f'/out/{gatk_hc_vcf_bname}', 'mode':'rw'}})


        #PLATYPUS

        #platypus_key = get_tool_key('platypus-variant')

        #Path(platypus_vcf_dir).touch()
        ts(); print('9. Skipping "platypus callVariants"')
        #client.containers.run(*platypus_key, f'platypus callVariants --bamFiles={gatk_recal_bam_bname} --refFile=./ref/{ref_file_bname} --output=./out/{platypus_vcf_bname}', working_dir='/data/', volumes={\
        #    ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, gatk_recal_bai_dir:{'bind':f'/data/{gatk_recal_bai_bname}', 'mode':'rw'}, platypus_vcf_dir:{'bind':f'/data/out/{platypus_vcf_bname}', 'mode':'rw'}})

        #MPILEUP

        bcftools_key = get_tool_key('bcftools')

        Path(bcf_mp_vcf_dir).touch()
        ts(); print('10. Running "bcftools mpileup | bcftools call"')
        client.containers.run(*bcftools_key, f'bash -c "bcftools mpileup -a FORMAT/DP -Ov -f ./ref/{ref_file_bname} {gatk_recal_bam_bname} | bcftools call -cv -f GQ --threads {samtools_thread_count} > ./out/{bcf_mp_vcf_bname}"', working_dir='/data/',  volumes={ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, bcf_mp_vcf_dir:{'bind':f'/data/out/{bcf_mp_vcf_bname}', 'mode':'rw'}})

        #FREEBAYES

        freebayes_key = get_tool_key('freebayes')

        Path(freebayes_vcf_dir).touch()
        ts(); print('11. Running "freebayes"')
        client.containers.run(*freebayes_key, f'bash -c "freebayes --genotype-qualities -f ./ref/{ref_file_bname} {gatk_recal_bam_bname} > ./out/{freebayes_vcf_bname}"', working_dir='/data/', volumes={\
            ref_file_tree:{'bind':f'/data/ref', 'mode':'rw'}, gatk_recal_bam_dir:{'bind':f'/data/{gatk_recal_bam_bname}', 'mode':'rw'}, freebayes_vcf_dir:{'bind':f'/data/out/{freebayes_vcf_bname}', 'mode':'rw'}})

    ref_fa_dir_list = [os.path.abspath(f) for f in glob.iglob(os.path.join(references_dir, '**/*'), recursive=True) if f.endswith('.fa')]
    print(ref_fa_dir_list)
    for ref_fa_dir in ref_fa_dir_list:
        ref_fa_bname = os.path.basename(ref_fa_dir)
        ref_fa_tree = os.path.split(ref_fa_dir)[0:-1][0]

        try:
            annotate(ref_fa_bname, ref_fa_tree)
            break
        except docker.errors.ContainerError as exc:
            errmsg = 'Input files reference and reads have incompatible contigs'
            if errmsg in exc.__str__():
                print('Reference mismatch.')
                continue
            else:
                raise exc


    #VCF MERGE - VCFTOOLS
    #vcftools_key = get_tool_key('vcftools')

    ts(); print('12. Running "bgzip & tabix"')

    bcftools_merged_vcf_bname = os.path.basename(bcftools_merged_vcf_dir)

    indiv_vcfs_dir = os.listdir(os.path.join(output_dir, f'individual_vcfs/{SM}/'))
    indiv_filtered_vcfs_dir = [os.path.splitext(vcf)[0] + '.filter' + os.path.splitext(vcf)[1] for vcf in indiv_vcfs_dir]
    #indiv_bgzipped_tree = os.path.join(indiv_vcfs_tree, 'bgzipped/')
    #Path(indiv_bgzipped_tree).mkdir(parents=True, exist_ok=True)
    for vcf in indiv_vcfs_dir:
        if vcf == SM+'_pl.vcf' or vcf.endswith('.filter.vcf'):
            continue
        if os.path.isfile(os.path.join(indiv_vcfs_tree, vcf)) and os.path.splitext(vcf)[1] == '.vcf':
            vcf_filtered_bname = os.path.splitext(vcf)[0] + '.filter' + os.path.splitext(vcf)[1]
            #vcf_filtered_dir = os.path.join(indiv_bgzipped_tree, vcf_filtered_bname)
            vcf_filtered_dir = os.path.join(indiv_vcfs_tree, vcf_filtered_bname)
            #print(vcf_filtered_dir)
            Path(vcf_filtered_dir).touch()
            bcftools_key = get_tool_key('bcftools')
            client.containers.run(*bcftools_key, f'bcftools view -i \"(%FILTER=\'PASS\' | %FILTER=\'.\') & (MIN(FMT/DP)>10 & MIN(FMT/GQ)>15)\" -o /out/{vcf_filtered_bname} {vcf}', working_dir='/data/', volumes={\
                indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'},\
                vcf_filtered_dir:{'bind':f'/out/{vcf_filtered_bname}', 'mode':'rw'}})
            os.remove(os.path.join(indiv_vcfs_tree, vcf))
            #subprocess.run(['bgzip', vcf_filtered_bname], cwd=indiv_vcfs_tree)

    #VCF MERGE - SURVIVOR

    survivor_key = get_tool_key('survivor')
    survivor_merged_vcf_bname = os.path.basename(survivor_merged_vcf_dir)

    Path(survivor_merged_vcf_dir).touch()
    ts(); print('13. Running "survivor"')
    client.containers.run(*survivor_key, f'bash -c "ls *vcf > sample_files; SURVIVOR merge sample_files 100 1 1 1 0 1 /out/{survivor_merged_vcf_bname}"', working_dir='/data/', volumes={\
        indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'}, survivor_merged_vcf_dir:{'bind':f'/out/{survivor_merged_vcf_bname}', 'mode':'rw'}})

    #VCF MERGE - SURVIVOR END

    for vcf_filt in os.listdir(indiv_vcfs_tree):
        if os.path.splitext(vcf_filt)[1] == '.vcf' and vcf_filt != SM+'_pl.vcf':
            subprocess.run(['bgzip', vcf_filt], cwd=indiv_vcfs_tree)

    for gz in os.listdir(indiv_vcfs_tree):
        if os.path.splitext(gz)[1] == '.gz':
            subprocess.run(['tabix', '-p', 'vcf', gz], cwd=indiv_vcfs_tree)

    Path(bcftools_merged_vcf_dir).touch()
    ts(); print('14. Running "bcftools merge"')

    freebayes_filtered_vcf_gz_bname = os.path.splitext(freebayes_vcf_bname)[0] + '.filter' + os.path.splitext(freebayes_vcf_bname)[1] + '.gz'
    bcf_mp_filtered_vcf_gz_bname = os.path.splitext(bcf_mp_vcf_bname)[0] + '.filter' + os.path.splitext(bcf_mp_vcf_bname)[1] + '.gz'
    gatk_hc_filtered_vcf_gz_bname = os.path.splitext(gatk_hc_vcf_bname)[0] + '.filter' + os.path.splitext(gatk_hc_vcf_bname)[1] + '.gz'

    client.containers.run(*bcftools_key, f'bash -c "bcftools merge --force-samples {freebayes_filtered_vcf_gz_bname} {bcf_mp_filtered_vcf_gz_bname} {gatk_hc_filtered_vcf_gz_bname} > /out/{bcftools_merged_vcf_bname}"', working_dir='/data/', volumes={\
        #freebayes_vcf_dir:{'bind':f'/data/{freebayes_vcf_bname}', 'mode':'rw'},\
        #bcf_mp_vcf_dir:{'bind':f'/data/{bcf_mp_vcf_bname}', 'mode':'rw'},\
        #platypus_vcf_dir:{'bind':f'/data/{platypus_vcf_bname}', 'mode':'rw'},\
        #gatk_hc_vcf_dir:{'bind':f'/data/{gatk_hc_vcf_bname}', 'mode':'rw'},\
        indiv_vcfs_tree:{'bind':f'/data/', 'mode':'rw'},\
        bcftools_merged_vcf_dir:{'bind':f'/out/{bcftools_merged_vcf_bname}', 'mode':'rw'}})

    ts(); print('15. Removing temporary files in ' + base_temp_dir)
    shutil.rmtree(base_temp_dir)

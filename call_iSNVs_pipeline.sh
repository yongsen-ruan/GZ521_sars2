metadata=./metadata ## with reference genome and primer information
refergenome=${metadata}/NC_045512.2.fasta ## reference genome
primers3end=${metadata}/primer5endrc_X.fa ## primer sequence for cutadapt
primers5end=${metadata}/primerX_5end.fa ## primer sequence for cutadapt
primerbed=${metadata}/primerbed_withstrand.txt ## primer position for iVar trim
samplelist=${metadata}/samplelist.txt ## sample list

rawfq=./exampledata ##.fastq format file
fastp_out=./fastp_out;mkdir -p ${fastp_out}
qc_afterfastp=./qc_afterfastp;mkdir -p ${qc_afterfastp}
cutadapt_out=./cutadapt_out;mkdir -p ${cutadapt_out}
bwa_out=./bwa_out;mkdir -p ${bwa_out}
logpath=./logpath;mkdir -p ${logpath}
mergefa_file=./${metadata}/mergefa_file;mkdir -p ${mergefa_file} ## file contain bam file path for different lane of the same library (optional)
merged_bam_output=./merged_bam_output;mkdir -p ${merged_bam_output}
ivartrim_output=./ivartrim_output;mkdir -p ${ivartrim_output}
ivar_out=./ivar_out;mkdir -p ${ivar_out}
mpileup_out=./mpileup_out;mkdir -p ${mpileup_out}
# step 1. quality control ##
for sp_path in $(ls ${rawfq}/*_1.fq.gz|sort); do
    sp=$(basename $sp_path _1.fq.gz)
    fastp -i ${rawfq}/${sp}_1.fq.gz -I ${rawfq}/${sp}_2.fq.gz \
		  -o ${fastp_out}/${sp}_1P.fq.gz -O ${fastp_out}/${sp}_2P.fq.gz \
		  -h ${qc_afterfastp}/${sp}.fastp.html -q 20 --thread 8 \
		  2>>${logpath}/fastp_cleaning.log 
done

## step 2. trim primer by cutadapt ##
for sp_path in $(ls ${rawfq}/*_1.fq.gz|sort); do
    sp=$(basename $sp_path _1.fq.gz)
    cutadapt -j 30 -a file:${primers3end} -A file:${primers3end} \
        -g file:${primers5end} -G file:${primers5end} \
        -o ${cutadapt_out}/${sp}_1.fq -p ${cutadapt_out}/${sp}_2.fq \
        ${fastp_out}/${sp}_1P.fq.gz ${fastp_out}/${sp}_2P.fq.gz -O 3 -m 30:30 --discard-untrimmed
done

## step 3. mapping to reference ##
##build index for reference genome
bwa index -a bwtsw ${refergenome} -p ${metadata}/NC_045512.2
##mapping fq format file to reference genome
for sp_path in $(ls ${rawfq}/*_1.fq.gz|sort); do
    sp=$(basename $sp_path _1.fq.gz)
    bwa mem -M -t 10 -R "@RG\tID:NULL\tLB:NULL\tPL:ILLUMINA\tSM:${sp}" ${metadata}/NC_045512.2 \
    ${cutadapt_out}/${sp}_1.fq ${cutadapt_out}/${sp}_2.fq \
    2>>${logpath}/bwamap.log \
    | samtools sort -O bam -o ${bwa_out}/${sp}.sort.bam --threads 4
    samtools index ${bwa_out}/${sp}.sort.bam -@ 2
done

## step 4. merge different lane of the same library and call mutation. ## 
for sample in $(sed '' ${samplelist});do
    ##merge different lane of the same library (optional)
    samtools merge -b ${mergefa_file}/${sample}_fqfile.txt -o ${merged_bam_output}/${sample}.merge.sort.bam -@ 4
    samtools index ${merged_bam_output}/${sample}.merge.sort.bam -@ 4
    ##further trim primer
    ivar trim -b ${primerbed} -p ${ivartrim_output}/${sample}.cutadapt.map.ivartrim -i ${merged_bam_output}/${sample}.merge.sort.bam -e -q 20 -m 30 -s 4
    samtools sort ${ivartrim_output}/${sample}.cutadapt.map.ivartrim.bam -o ${ivartrim_output}/${sample}.cutadapt.map.ivartrim.sorted.bam -@ 8 
    rm ${ivartrim_output}/${sample}.cutadapt.map.ivartrim.bam 
    ##call mutation and export mpileup file
    samtools mpileup -aa -d 600000 -Q 20 -f ${refergenome} ${ivartrim_output}/${sample}.cutadapt.map.ivartrim.sorted.bam \
        |ivar variants -p ${ivar_out}/${sample} -q 20 -t 0.01 -r ${refergenome} -m 10
    samtools mpileup -aa -d 600000 -Q 20 -f ${refergenome} ${ivartrim_output}/${sample}.cutadapt.map.ivartrim.sorted.bam \
        |sequenza-utils pileup2acgt -p - -o ${mpileup_out}/${sample}.mpileup_q20.txt
done




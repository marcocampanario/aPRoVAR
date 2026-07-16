#!/bin/bash

# =========================================================
# Do not change
DRAGEN_HASH_TABLE="~/hg38.alt_masked.cnv.graph.hla.rna-10-r4.0-1/"
TEMPD="~/tmp/"
NIRVANADBS="~/NirvanaData/"
ASSEMBLY="GRCh38"
MLMODEL="/opt/dragen/4.3.6/resources/ml_model/"

declare -A PROCESSED_RGSM
# =========================================================
# Change here
INPUT_FOLDER="~/"
FASTQ_LIST="${FOLDER}/fastq_list.csv"
OUPTUT_FOLDER="${FOLDER}/results"
VCTARGETBED="~/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed"
# =========================================================

{ 
    read # skip header
    while IFS=, read _ RGSM _ _ __ SEX
        do

        # Get RGSM from FASTQ list
        # RGSM=$(echo $LINHA | cut -d',' -f2)

        # Verify if the sample with RGSM was already processed
        if [[ -n "${PROCESSED_RGSM[$RGSM]}" ]];then
            echo "RGSM $RGSM already processed. Skipping..."
            continue
        fi

        # Verify if the sample has SEX declared
        
        if [[ -z "$SEX" ]]; then
        echo "Warning: SEX not defined for $RGSM, skipping..."
        continue
        fi

        # Mark RGSM as processed
        PROCESSED_RGSM["$RGSM"]=1

        # Output directory
        FOLDER_OUT="${OUPTUT_FOLDER}/${RGSM}_OUT"
        mkdir -p $FOLDER_OUT
        
        PREFIX="${RGSM}"

        # Define the input sources, select fastq list, fastq, bam, or cram
        INPUT_FASTQ_LIST="
          --fastq-list $FASTQ_LIST \
          --fastq-list-sample-id $RGSM \
        "

        INPUT_OPTIONS="
            -r $DRAGEN_HASH_TABLE \
               $INPUT_FASTQ_LIST \
        "

       #ANNOT_OPTIONS=" 
       #     --enable-variant-annotation true \
       #     --variant-annotation-data $NIRVANADBS \
       #     --variant-annotation-assembly $ASSEMBLY \
       # "

        OUTPUT_OPTIONS="
            --output-directory $FOLDER_OUT \
            --output-file-prefix $PREFIX \
        "

        MA_OPTIONS="
            --enable-map-align true \
            --enable-sort true \
            --enable-duplicate-marking true \
            --enable-hla true \
            --hla-enable-class-2 true \
            --intermediate-results-dir $TEMPD \
            --output-format BAM \
            --enable-bam-indexing true \
            --enable-map-align-output true \
            --Aligner.sec-aligns 0 \
            --Aligner.sec-aligns-hard 1 \
            --Aligner.supp-aligns 0 \
            --Aligner.supp-as-sec 1 \
        "

        CNV_OPTIONS="
            --enable-cnv false \
        "

        SEX_OPTIONS="
        --sample-sex $SEX \
        "

        SNV_OPTIONS="
            --enable-variant-caller true \
            --vc-target-bed $VCTARGETBED \
            --vc-emit-ref-confidence GVCF \
            --vc-enable-vcf-output true \
            --vc-compact-gvcf true \
            --vc-ml-dir $MLMODEL \
            --vc-ml-enable-recalibration true \
        "

        SV_OPTIONS="
            --enable-sv false \
        "

        # Construct final command line
        CMD="
        /opt/dragen/4.3.6/bin/dragen \
            $INPUT_OPTIONS \
            $OUTPUT_OPTIONS \
            $ANNOT_OPTIONS \
            $MA_OPTIONS \
            $CNV_OPTIONS \
            $SEX_OPTIONS \
            $SNV_OPTIONS \
            $SV_OPTIONS \
        "

        # Execute
        #echo $CMD
        echo ">>> $RGSM"
        echo -e "------\n"
        $CMD

        done 
} < $FASTQ_LIST

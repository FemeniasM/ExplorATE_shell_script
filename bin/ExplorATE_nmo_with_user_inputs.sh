#!/bin/bash
# Using getopt

abort()
{
    echo >&2 '
-- ABORTED --
'
    echo "An error occurred, check the input files..." >&2
    exit 1
}

trap 'abort' 0

set -e
#############################################################################
# 	This script run ExplorATE pipeline with input files provided by the user 		
##############################################################################


threads=1
kmer=31
salmon_path=="salmon"
bedtools="bedtools"
seqkit=$(dirname "${BASH_SOURCE[0]}")/seqkit
RM2Bed=$(dirname "${BASH_SOURCE[0]}")/RM2Bed.py
al=''
ovres='higher_score'
split_by='subclass'
Wlr="80,80,80"
annot_by='transcripts'
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`

print_usage () {
    echo "$0
Usage:  
 ExplorATE ${yellow}nmo ${green}[flags]${reset}

 Flags:
  ${green}-p${reset} threads [N] (1 default) 
  ${green}-k${reset} kmer [N] (31 default)
  ${green}-b${reset} bedtools binary path (bedtools default)
  ${green}-r${reset} RepeatMasker output file
  ${green}-l${reset} folder with fastq files 
  ${green}-e${reset} library format ['pe']['se']
  ${green}-o${reset} output directory path 
  ${green}-t${reset} fasta transcriptome 
  ${green}-s${reset} salmon binary path (salmon default)
  ${green}-n${reset} gene annotation file in outfmt 6 
  ${green}-d${reset} TransDecoder gff3 file 
  ${green}-w${reset} Wicker-like rule comma separated ('80,80,80' default)
  ${green}-v${reset} overlap resolution criterion ['higher_score']['longer_element']['lower_divergence'] ('higher_score' default)
  ${green}-x${reset} split by ['name']['family']['subclass']
  ${green}-q${reset} annotate target TEs by ['transcripts']['fragments'] ('transcripts' default)
  ${green}-a${reset} .align file (if use 'lower_divergence' as overlap resolution)
  ${green}-h${reset} help" 
    exit 1
}

abort() {
    echo -e >&2 '
===============
/// ABORTED ///
===============
'
    echo -e "${red}An error occurred. Exiting...${reset}" >&2
#    rm -r $outfolder/temp
    exit 1
}

trap 'abort' 0
set -e


file_extension_error ()
{
echo -e "${red}FileExtensionError: Invalid extension${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}"
return     
}

file_name_error ()
{
echo -e "${red}Filename Error: Paired end file names should contain _R1 _R2${reset}"
echo -e "${green}Example: sample_R1.fq.gz, sample_R2.fq.gz${reset}"
return 
}

echo "${green}=============================
${yellow}|||       ExplorATE       |||
${green}=============================${reset}"
while getopts ":p:k:b:r:l:e:o:t:s:n:d:w:v:x:q:a:h:" opt; do
    case $opt in

        p)
            threads="$OPTARG"
            echo "-p <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;
        b)
            bedtools=`realpath $OPTARG`
            echo "-b <bedtools binary> = $bedtools"
            ;;
        r)
            RM_file=`realpath $OPTARG`
            echo "-r <RepeatMasker output file> = $RM_file"
            ;;
        l)
            lib_folder=`realpath $OPTARG`
            echo "-l <fastq libraries path> = $lib_folder"
            ;;
        e)
            lib_format="$OPTARG"
            echo "-e <libraries format 'pe' | 'se'> = $lib_format"
            ;;
        o)
            outfolder=`realpath $OPTARG`
            echo "-o <Output files Path> = $outfolder"
            ;;
        t)
            fasta_trme=`realpath $OPTARG`
            echo "-t <fasta transcriptome> = $fasta_trme"
            ;;
        s)
            salmon_path=`realpath $OPTARG`
            echo "-s <salmon path> = $salmon_path"
            ;;
        n)
            outfmt6_file=`realpath $OPTARG`
            echo "-n <gene annotation file> = $outfmt6_file"
            ;;
        d)
            TransDecoder_file=`realpath $OPTARG`
            echo "-d <TransDecoder gff3 file> = $TransDecoder_file"
            ;;
        w)
            Wlr="$OPTARG"
            echo "-w <Wicker-like rule> = $Wlr"
            ;;
        v)
            ovres="$OPTARG"
            echo "-v <overlap resolution 'higher_score'|'longer_element'|'lower_divergence'> = $ovres"
            ;;
        x)
            split_by="$OPTARG"
            echo "-x < 'name'|'family'|'subclass'> = $split_by"
            ;;
        q)
            annot_by="$OPTARG"
            echo "-q < 'transcripts'|'fragments'> = $annot_by"
            ;;
        a)
            a="$OPTARG"
            echo "-a <*.align file> = $al"
            ;;
        h)
            print_usage && exit 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage && exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage && exit 1
            ;;
            
    esac
done

if [ -z "$outfolder" -o -z "$RM_file" -o -z "$TransDecoder_file" -o -z "$fasta_trme" -o -z "$threads" -o -z "$lib_folder" -o -z "$lib_format"  ]; then
    echo -e "${red}Error: missing required argument(s)${reset}"
    print_usage && exit 1
fi

echo "[$(printf '%(%F %T)T\n')][INFO] checking fastq files extensions.."
if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        echo "[$(printf '%(%F %T)T\n')][INFO] libraries extension format ok..."
        else 
        file_extension_error
        echo -e "${red}Error: check extension format: $fn ${reset}"
        exit 1
    fi
    done
else 
 if [[ $lib_format == 'pe' ]]; then
    for fn in $lib_folder/*R1*; do
        if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
         echo "[$(printf '%(%F %T)T\n')][INFO] $fn: extension format ok..."
            sample_name=$(basename ${fn} | sed 's/.fastq.gz\|.fq.gz\|.fastq\|.fq//g')
            if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
            echo "[$(printf '%(%F %T)T\n')][INFO] $fn: name format ok..."
            else
               file_name_error
               echo -e "${red}Error: check name format: $fn ${reset}"
               exit 1
            fi
         else
         file_extension_error
         echo -e "${red}Error: check extension format: $fn ${reset}"
         exit 1
         fi
     done
 else
    echo -e "${red}Error: Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi
fi

mkdir -p $outfolder/temp
cd $outfolder/temp

echo -e "[$(printf '%(%F %T)T\n')][INFO] recording transcript lengths"
$seqkit fx2tab --length --name -i $fasta_trme | sort -k1,1 -k2,2n > transcripts_len.bed

echo -e "[$(printf '%(%F %T)T\n')][INFO] writing UTR regions"
awk -v OFS="\t" '($3 ~ /five_prime_UTR|three_prime_UTR/) {print $1,$4,$5,$3}' $TransDecoder_file | sort -t $'\t' -k1,1 -k2,2n | $bedtools merge -c 4 -o first > UTRs.bed
awk -F".p" '{print $1}' $outfmt6_file | sort | uniq > IDGenes.txt
awk 'FNR==NR { a[$1]=$1; next } ($1 in a)' IDGenes.txt UTRs.bed >UTR_genes.bed

#Overlap resolution
if [[ $ovres == 'higher_score' || $ovres == 'longer_element' || $ovres == 'lower_divergence' ]]; then
echo -e "${yellow}
Resolving overlaps with RM2Bed
==============================${reset}"

python3 $RM2Bed --out_dir . --out_prefix RM_ovres_RM2Bed --sort_criterion 'size' --split 'subclass' --ovlp_resolution $ovres $RM_file $al

echo -e "[$(printf '%(%F %T)T\n')][INFO] filtering TE sequences non-overlapped with genes"
awk -v OFS="\t" '{print $1, $2, $3, $4":"$7":"$8 }' RM_ovres_RM2Bed_rm.bed | sort -k1,1 -k2,2n | $bedtools merge -c 4 -o first > RM_overes_merge_sort.bed
$bedtools intersect -wb -wa -a UTR_genes.bed -b RM_overes_merge_sort.bed -sorted | awk '!($5="")' > repeats_in_UTRs.txt
awk '{print $1}' repeats_in_UTRs.txt | sort | uniq > seqIDGenesWithRepinUTRs.txt
$seqkit grep $fasta_trme -v -w 0 -f seqIDGenesWithRepinUTRs.txt > NoGenes.fa

echo -e "[$(printf '%(%F %T)T\n')][INFO] creating references to apply Wicker-like rule"
grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" RM_ovres_RM2Bed_rm.bed \
| awk '{print $1, $4, $7, $8, $5}' > BED_wout_ovlp.bed
grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" $RM_file >RMclean.out

    if [[ $split_by == 'name' ]];then
    echo -e "[$(printf '%(%F %T)T\n')][INFO] making names-based references for Wicker-like rule"
    awk 'FNR==NR{a[$1]=$2;next}{ print $1"::"$2,$3, a[$1], $5/a[$1]*100}' transcripts_len.bed BED_wout_ovlp.bed > trLen_perRep.txt
    awk 'NR>3 {print $5"::"$10, $2}' RMclean.out | awk '{a[$1]++}{b[$1]=b[$1]+$2}END{for (i in a) { print i, 100-(b[i]/a[i])} }' > PerIdentity_ok.txt
    awk 'NR>3 {print $5, $6, $7,$5"::"$10}' RMclean.out > RMclean.bed
    fi
    if [[ $split_by == 'subclass' ]];then
    echo -e "[$(printf '%(%F %T)T\n')][INFO] making subclass-based references for Wicker-like rule"
    awk 'FNR==NR{a[$1]=$2;next}{ print $1"::"$3,$3, a[$1], $5/a[$1]*100}' transcripts_len.bed BED_wout_ovlp.bed > trLen_perRep.txt
    awk 'NR>3{gsub("/","\t",$11);print $0}' RMclean.out | awk 'NR>3 {print $5"::"$11, $2}' | awk '{a[$1]++}{b[$1]=b[$1]+$2}END{for (i in a) { print i, 100-(b[i]/a[i])} }' > PerIdentity_ok.txt
    awk 'NR>3{gsub("/","\t",$11);print $0}' RMclean.out | awk 'NR>3 {print $5, $6, $7, $5"::"$11}' > RMclean.bed
    fi
    if [[ $split_by == 'family' ]];then
    echo -e "[$(printf '%(%F %T)T\n')][INFO] making families-based references for Wicker-like rule"
    awk 'FNR==NR{a[$1]=$2;next}{ print $1"::"$4,$3, a[$1], $5/a[$1]*100}' transcripts_len.bed BED_wout_ovlp.bed > trLen_perRep.txt
    awk 'NR>3{gsub("/","\t",$11);print $0}' RMclean.out | awk 'NR>3 {print $5"::"$12, $2}' | awk '{a[$1]++}{b[$1]=b[$1]+$2}END{for (i in a) { print i, 100-(b[i]/a[i])} }' > PerIdentity_ok.txt
    awk 'NR>3{gsub("/","\t",$11);print $0}' RMclean.out | awk 'NR>3 {print $5, $6, $7, $5"::"$12}' > RMclean.bed
    fi

awk 'FNR==NR{a[$1]=$2 ;next}{ print $0, a[$1]}' PerIdentity_ok.txt trLen_perRep.txt >WickerLike_tbl.txt
echo -e "[$(printf '%(%F %T)T\n')][INFO] applying Wicker-like rule"
awk -v perId=$(echo $Wlr | awk -F',' '{print $1}') -v perLen=$(echo $Wlr | awk -F',' '{print $2}') -v trLen=$(echo $Wlr | awk -F',' '{print $3}') \
'$5>=PerId && $4>=perLen && $3>=trLen' WickerLike_tbl.txt > TranscriptsID.txt
awk -F'::' '{print $1}' TranscriptsID.txt > targetTEs_ids.txt

    if [[ $annot_by == 'transcripts' ]];then
    echo -e "[$(printf '%(%F %T)T\n')][INFO] target TEs annotation by transcripts"
    seqkit grep NoGenes.fa -w 0 -f  targetTEs_ids.txt > targetTEs.fa
    grep '>' targetTEs.fa | sed 's/>//g' | awk '{print $1}' >targetTEs_tr_ids_noGenes.txt
    grep -w -f targetTEs_tr_ids_noGenes.txt TranscriptsID.txt | awk '{gsub("::",";", $1); print $1}' > references.csv
    else
        if [[ $annot_by == 'fragments' ]];then
        echo -e "[$(printf '%(%F %T)T\n')][INFO] target TEs annotation by fragments"
        awk '{print $1}' TranscriptsID.txt >fragmentsIDs.txt
        grep -v -w -f seqIDGenesWithRepinUTRs.txt RMclean.bed | grep -w -f fragmentsIDs.txt | awk '{print $1, $2, $3, $4":"$2$3}' \
        | sort -k1,1 -k2,2n | $bedtools merge -c 4 -o first >  targetTEs.bed
        $bedtools getfasta -fi $fasta_trme -bed targetTEs.bed -nameOnly -fo  targetTEs.fa
        grep '>' targetTEs.fa | sed 's/>//g' | awk -F':' '{print $0";"$3}' > references.csv
        grep '>' targetTEs.fa | sed 's/>//g' | awk -F':' '{print $1}' > targetTEs_tr_ids_noGenes.txt
        fi
    fi
fi

echo -e "[$(printf '%(%F %T)T\n')][INFO] removing redundant sequences"
$seqkit rmdup -s targetTEs.fa | $seqkit seq -m $kmer -g > targetTEs_clean.fa
$seqkit grep $fasta_trme -v -w 0 -f targetTEs_tr_ids_noGenes.txt | $seqkit rmdup -s | $seqkit seq -m $kmer -g > decoys.fa
grep '>' decoys.fa | sed 's/>//g' | awk '{print $1}' > decoys.txt
echo -e "[$(printf '%(%F %T)T\n')][INFO] making references for Salmon"
cat targetTEs_clean.fa decoys.fa > trmeSalmon.fa
rm targetTEs.fa targetTEs_clean.fa decoys.fa 

echo -e "${yellow}
Indexing in Salmon
==================${reset}"

$salmon_path index -p $threads -t trmeSalmon.fa -k $kmer -i SALMON_INDEX --decoys decoys.txt

echo -e "${yellow}
Quantifying with Salmon
=======================${reset}"

if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        sample_name=$(basename ${fn} | sed 's/.fastq.gz\|.fq.gz\|.fastq\|.fq//g')
        $salmon_path quant -i SALMON_INDEX -l A --gcBias --useVBOpt -r ${fn} -p $threads --validateMappings -o $outfolder/quant_out/$sample_name
        else 
        file_extension_error
        exit 1
    fi
    done
else 
 if [[ $lib_format == 'pe' ]]; then
    for fn in $lib_folder/*R1*; do
        if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
         sample_name=$(basename ${fn} | sed 's/_R1.fastq.gz\|_R1.fq.gz\|_R1.fastq\|_R1.fq//g')
            if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
             ext=$(basename ${fn} | awk -F'_R1' '{print $2}')
             R1=${sample_name}_R1${ext}
             R2=${sample_name}_R2${ext}
             $salmon_path quant -i SALMON_INDEX -l A --gcBias --useVBOpt -1 $lib_folder/$R1 -2 $lib_folder/$R2 -p $threads --validateMappings -o $outfolder/quant_out/$sample_name
            else
               file_name_error
               exit 1
            fi
         else
         file_extension_error
         exit 1
         fi
     done
 else
    echo -e "\n${red} Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi
fi

mv references.csv $outfolder
mv repeats_in_UTRs.txt $outfolder
echo "[$(printf '%(%F %T)T\n')][INFO] removing temporary files"
rm -r $outfolder/temp


trap : 0
echo >&2 "${green}
======================================================
>         ExplorATE finished successfully            <
======================================================
"


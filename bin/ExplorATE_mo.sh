#!/bin/bash
# Using getopt

#########################################################################################
# This script is based on Selective Alignment of sequences, 
# (Srivastava, A., Malik, L., Sarkar, H. et al. Alignment 
# and mapping methodology influence transcript abundance 
# estimation. Genome Biol 21, 239 (2020). 
# https://doi.org/10.1186/s13059-020-02151-8)
#----------------------------------------------------
# It assumes awk, bedtools and Salmon is 
# available.
# We have tested this script with awk 4.1.4, 
# bedtools v2.30.0 on an Ubuntu system. 
# 
# This script uses SeqKit (Wei Shen, Copyright ©2016-2019 Oxford Nanopore Technologies.
# W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for 
# FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962)
#
# This script uses RM2Bed.py from RepeatMaster distribution (Authors: David Ray 
# and Robert Hubley) for overlaps resolution.
########################################################################################

threads=1
kmer=31
salmon_path="salmon"
bedtools="bedtools"
seqkit=$(dirname "${BASH_SOURCE[0]}")/seqkit
RM2Bed=$(dirname "${BASH_SOURCE[0]}")/RM2Bed.py
al=''
ovres=''
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`

# Argument Parsing
print_usage () {
    echo "$0
Usage:  
ExplorATE ${yellow}mo${green} [flags]${reset}

Flags:
   ${green}-p${reset} threads [N] (1 default) 
   ${green}-k${reset} kmer [N] (31 default)
   ${green}-c${reset} chromosome alias file
   ${green}-b${reset} bedtools binary path (bedtools default)
   ${green}-s${reset} salmon binary path (salmon default)
   ${green}-f${reset} fasta genome 
   ${green}-g${reset} gtf file 
   ${green}-r${reset} RepeatMasker .out file 
   ${green}-e${reset} library format ['pe']['se'] 
   ${green}-l${reset} folder with fastq files
   ${green}-o${reset} output directory path 
   ${green}-t${reset} fasta transcriptome ** 
   ${green}-u${reset} RepeatMasker .out file from transcriptome ** 
   ${green}-v${reset} overlap resolution ['higher_score']['longer_element']['lower_divergence'] 
   ${green}-a${reset} .align file (if use 'lower_divergence' as overlap resolution) 
   ${green}-h${reset} help
** if the TEs targets are based on a de novo transcriptome" 
    exit 1
}

abort() 
{
    echo -e >&2 '
===============
/// ABORTED ///
===============
'
    echo -e "${red}An error occurred. Exiting...${reset}" >&2
   [ -f $outfolder/temp ] && rm -r $outfolder/temp
    exit 1
}

trap 'abort' 0
set -e

file_extension_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] FileExtensionError: Invalid extension${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}"
return     
}

file_name_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] Filename Error: Paired end file names should contain _R1 _R2${reset}"
echo -e "${green}Example: sample_R1.fq.gz, sample_R2.fq.gz${reset}"
return 
}

echo "${green}=============================
${yellow}|||       ExplorATE       |||
${green}=============================${reset}"
while getopts ":a:p:b:c:s:f:g:r:o:k:l:t:u:v:e:h:" opt; do
    case $opt in

        b)
            bedtools=`realpath $OPTARG`
            echo "-b <bedtools binary> = $bedtools"
            ;;
        f)
            fasta_genome=`realpath $OPTARG`
            echo "-f <fasta genome> = $fasta_genome"
            ;;
        c)
            chrAl=`realpath $OPTARG`
            echo "-c <chromosome alias file> = $chrAl"
            ;;
        g)
            gtf_genome=`realpath $OPTARG`
            echo "-g <gtf genome> = $gtf_genome"
            ;;
        r)
            RM_genome=`realpath $OPTARG`
            echo "-r <RepeatMasker file> = $RM_genome"
            ;;
        s)
            salmon_path=`realpath $OPTARG`
            echo "-s <salmon path> = $salmon_path"
            ;;
        l)
            lib_folder=`realpath $OPTARG`
            echo "-l <fastq libraries path> = $lib_folder"
            ;;
        e)
            lib_format="$OPTARG"
            echo "-e <libraries format 'pe' | 'se'> = $lib_format"
            ;;

        p)
            threads="$OPTARG"
            echo "-p <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;

        o)
            outfolder=`realpath $OPTARG`
            echo "-o <Output files Path> = $outfolder"
            ;;
        t)
            fasta_trme=`realpath $OPTARG`
            echo "-t <fasta transcriptome> = $fasta_trme"
            ;;
        u)
            RM_trme=`realpath $OPTARG`
            echo "-u <gtf RepeatMasker from transcriptome> = $RM_trme"
            ;;
        v)
            ovres="$OPTARG"
            echo "-v <overlap resolution 'higher_score'|'longer_element'|'lower_divergence'> = $ovres"
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

if [ -z "$bedtools" -o -z "$fasta_genome" -o -z "$gtf_genome" -o -z "$RM_genome" -o -z "$salmon_path" -o -z "$lib_format" -o -z "$lib_folder" -o -z "$threads" -o -z "$kmer" -o -z "$outfolder" ]
then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] missing required argument(s)${reset}"
    print_usage && exit 1
fi

bedtools_ver=$($bedtools --version | awk -F" v" '{print $2}')
salmon_ver=$($salmon_path --version | awk -F" " '{print $2}')

if { echo "$bedtools_ver"; echo "2.29.0"; } | sort --version-sort --check=silent; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] ExplorATE requires bedtools version 2.29.1 or later. 
    You are using betools version $bedtools_ver ${reset}"
    print_usage && exit 1
fi

if { echo "$salmon_ver"; echo "1.3.0"; } | sort --version-sort --check=silent; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] ExplorATE requires salmon version 1.4.0 or later. 
    You are using Salmon version $salmon_ver ${reset}"
    print_usage && exit 1
fi

mkdir -p $outfolder/temp
cd $outfolder/temp

echo -e "${yellow}
Making TEs references
=====================${reset}"
echo "---------------------------
Genomic features processing
---------------------------"
echo "[$(printf '%(%F %T)T\n')][INFO] recording chromosome lengths..."

$seqkit fx2tab --length --name -i $fasta_genome | sort -k1,1 -k2,2n > genome_len.bed

if [[ -e $fasta_trme ]]; then
echo "
-------------------------------------------
Searching TE targets based in transcriptome
-------------------------------------------"
        if [[ $ovres == 'higher_score' || $ovres == 'longer_element' || $ovres == 'lower_divergence' ]]; then
        
echo -e "${yellow}
Resolving overlaps with RM2Bed
==============================${reset}"
        
        python3 $RM2Bed --out_dir . --out_prefix RM_ovres --sort_criterion 'size' --ovlp_resolution $ovres $RM_trme $al
        
        grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" RM_ovres_rm.bed \
        | awk -v OFS='\t' '{print $1, $2, $3, $1":"$2":"$3":"$4"/"$7"/"$8}' | sort -k1,1 -k2,2n | $bedtools merge -c 4 -o first > RMtrme.bed
        
        else
        
        grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" $RM_trme \
        | awk -v OFS='\t' 'NR>3 {print $5, $6, $7, $5":"$6":"$7":"$10"/"$11}' | sort -k1,1 -k2,2n | $bedtools merge -c 4 -o first > RMtrme.bed
        
        fi

$bedtools getfasta -fi $fasta_trme -bed RMtrme.bed -nameOnly -fo RM_trme.fa

$seqkit rmdup -s RM_trme.fa | $seqkit seq -m $kmer -g > tarTEs.fa

grep '>' tarTEs.fa | sed 's/>//g' |awk -F":" '{print $0";"$4}' > references.csv

echo -e "[$(printf '%(%F %T)T\n')][INFO] ${green} $(grep -c '>' tarTEs.fa) target TEs found based in transcriptome${reset}"

else

echo "
------------------------------------
Searching TE targets based in genome
------------------------------------"

    if [[ $ovres == 'higher_score' || $ovres == 'longer_element' || $ovres == 'lower_divergence' ]]; then
    
echo -e "${yellow}
Resolving overlaps with RM2Bed$
==============================${reset}"
    
    python3 $RM2Bed --out_dir . --out_prefix RM_ovres --sort_criterion 'size' --ovlp_resolution $ovres $RM_genome $al
    
    echo "[$(printf '%(%F %T)T\n')][INFO] Removing short repeats (small RNAs, Satellites, Low complexity and Simple repeats, etc) from reference"
    
    grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" RM_ovres_rm.bed \
    | awk -v OFS='\t' '{print $1, $2, $3, $1":"$2":"$3":"$4"/"$7"/"$8}' | sort -k1,1 -k2,2n > RMgen.bed
    
    else
    
    echo "[$(printf '%(%F %T)T\n')][INFO] Removing short repeats (small RNAs, Satellites, Low complexity and Simple repeats, etc) from reference"
    
    grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" $RM_genome \
    | awk -v OFS='\t' 'NR>3 {print $5, $6, $7, $5":"$6":"$7":"$10"/"$11}' | sort -k1,1 -k2,2n > RMgen.bed
    
    fi
    
    if [[ -e $chrAl ]]; then
    
    echo "[$(printf '%(%F %T)T\n')][INFO] homogenizing chromosome names from chromosome alias file"
    
    awk '{print $1}' genome_len.bed > ChrNamGen.txt
    awk '{print $1}' RMgen.bed | sort | uniq > ChrNamRM.txt
    grep -v -w -f ChrNamRM.txt ChrNamGen.txt > NamesNoMatch.txt
    grep -w -f NamesNoMatch.txt $chrAl | awk '{print $2, $1}' > chromAlias_noMatch.txt
    cat RMgen.bed | sed -f <(printf 's/%s/%s/g\n' $(<chromAlias_noMatch.txt)) > RM_gen.bed
    mv RM_gen.bed RMgen.bed
    rm *.txt
        
    fi
    
    echo -e "[$(printf '%(%F %T)T\n')][INFO] Making references.csv file"
    
    $bedtools merge -c 4 -o first -i RMgen.bed | sort -k1,1 -k2,2n > RMgen_merged.bed
    
    awk '{print $4}' RMgen.bed | awk -F":" -v OFS=";" '{print $1":"$2":"$3":"$4, $4}' > references.csv

    echo -e "[$(printf '%(%F %T)T\n')][INFO] Extract target sequences"

    $bedtools getfasta -fi $fasta_genome -bed RMgen.bed -nameOnly -fo target1.fa
   
    echo -e "[$(printf '%(%F %T)T\n')][INFO] Removing duplicated and short sequences from target file"

    $seqkit rmdup -s target1.fa | $seqkit seq -m $kmer -g > target.fa

    rm target1.fa 
    mv RMgen.bed references.bed
    
    echo -e "[$(printf '%(%F %T)T\n')][INFO] ${green}$(grep -c '>' target.fa) ${reset}target TEs found based in genome"

    echo -e "[$(printf '%(%F %T)T\n')][INFO] Extract decoys reference"

    $bedtools complement -i RMgen_merged.bed -g  genome_len.bed -L | sort -k1,1 -k2,2n | $bedtools merge \
    | awk -v OFS='\t' '{print $1, $2, $3, "d_"$1"_"$2"_"$3}' > decoys.bed

    $bedtools getfasta -fi $fasta_genome -bed decoys.bed -nameOnly -fo decoys1.fa

    echo -e "[$(printf '%(%F %T)T\n')][INFO] Removing duplicated and short sequences from decoy file"

    $seqkit rmdup -s decoys1.fa | $seqkit -is replace -p "n" -r "" | $seqkit seq -g -m $((2 * $kmer)) > decoys.fa

    rm decoys1.fa RMgen_merged.bed genome_len.bed decoys.bed
    
    awk '/^>/ {print $0; next}' decoys.fa | sed 's/^>//g' > decoys.txt
    
    echo -e "[$(printf '%(%F %T)T\n')][INFO] ${green}$(wc -l decoys.txt | awk '{print $1}') ${reset}decoys sequences generated"
    
    echo "[$(printf '%(%F %T)T\n')][INFO] Making reference gentrome"

    cat target.fa decoys.fa > trmeSalmon.fa
    
    rm decoys.fa target.fa 

fi

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
    for fn in $lib_folder/*_R1*; do
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
echo "[$(printf '%(%F %T)T\n')][INFO] removing temporary files"
rm -r $outfolder/temp/

trap : 0
echo >&2 "${green}
======================================================
          ExplorATE finished successfully            
======================================================
"

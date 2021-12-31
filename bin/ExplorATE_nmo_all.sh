#!/bin/bash
# Using getopt

######################################################################
# 	This script generates input files for ExplorATE. 		
# __________________________________________________________________
# TransDecoder, BLAST, HMMER and RepeatMasker must be installed.
# In addition, you must provide the UniProt/SwissProt and Pfam and
# databases, and a repeats library (such as Dfam, RepBase or specie-
# specific library).
# ------------------------------------------------------------------
# Tested with BLAST 2.2.29, HMMER 3.2.1, TransDecoder 5.5.0 and
# RepeatMasker 4.1.0
######################################################################



BLASTP_path=blastp
HMMER_path=hmmscan
RM_path=RepeatMasker
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
  ExplorATE ${yellow}nmo_all ${green}[flags]${reset}
  
  Flags:
    ${green}-p${reset} threads [N] (1 default)
    ${green}-k${reset} kmer [N] (31 default)
    ${green}-b${reset} bedtools binary path (bedtools default)
    ${green}-r${reset} RepeatMasker binary path (RepeatMasker default)
    ${green}-l${reset} folder with fastq files
    ${green}-e${reset} library format ['pe']['se']  
    ${green}-o${reset} output directory path 
    ${green}-t${reset} fasta transcriptome
    ${green}-s${reset} salmon binary path (salmon default)
    ${green}-n${reset} blastp binary path (blastp default)
    ${green}-m${reset} hmmer binary path (hmmscan default)
    ${green}-d${reset} TransDecoder directory path
    ${green}-u${reset} swissprot database
    ${green}-f${reset} Pfam database
    ${green}-i${reset} RepeatMasker -lib argument (custom library)
    ${green}-j${reset} RepeatMasker -species argument (species library)
    ${green}-w${reset} Wicker-like rule comma separated ('80,80,80' default)
    ${green}-v${reset} overlap resolution criterion ['higher_score']['longer_element']['lower_divergence']
    ${green}-x${reset} split by ['name']['family']['subclass'] ('subclass' default)
    ${green}-q${reset} annotate target TEs by ['transcripts']['fragments'] ('transcripts' default)
    ${green}-a${reset} .align file (if use 'lower_divergence' as overlap resolution)
    ${green}-h${reset} help " 
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
RepeatMasker_lib_error ()
{
echo -e "${red}Error: A repeats library is required for RepeatMasker.${reset}"
echo -e "${green}You can enter a custom library ('-i' flag on ExplorATE), or use a RepeatMasker '-species' argument (-j flag on ExplorATE). 
See the RepeatMasker documentation for more details.${reset}"
return 
}

echo "${green}=============================
${yellow}|||       ExplorATE       |||
${green}=============================${reset}"
while getopts ":p:k:b:r:l:e:i:j:o:t:s:n:m:d:w:f:u:q:v:x:a:h:" opt; do
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
            RM_path=`realpath $OPTARG`
            echo "-r <RepeatMasker binary path> = $RM_path"
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
            BLASTP_path=`realpath $OPTARG`
            echo "-n <blastp binary path> = $BLASTP_path"
            ;;
        m)
            HMMER_path=`realpath $OPTARG`
            echo "-m <hmmer binary path> = $HMMER_path"
            ;;
        d)
            TransDecoder_path=`realpath $OPTARG`
            echo "-d <TransDecoder directory path> = $TransDecoder_path"
            ;;
        u)
            sprot_db=`realpath $OPTARG`
            echo "-u <swissprot database> = $sprot_db"
            ;;
        f)
            pfam_db=`realpath $OPTARG`
            echo "-f <Pfam database> = $pfam_db"
            ;;
        i)
            TEs_lib=`realpath $OPTARG`
            echo "-i <Repeats (custom library)> = $TEs_lib"
            ;;
        j)
            TEs_sp="$OPTARG"
            echo "-j <Repeats (species library)> = $TEs_sp"
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
            echo "-x <'name'|'family'|'subclass'> = $split_by"
            ;;
        q)
            annot_by="$OPTARG"
            echo "-q <'transcripts'|'fragments'> = $annot_by"
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

if [ -z "$outfolder" -o -z "$RM_path" -o -z "$TransDecoder_path" -o -z "$fasta_trme" -o -z "$threads" -o -z "$lib_folder" -o -z "$lib_format" ]
then
    echo -e "\n${red}Error: missing required argument(s)${reset}\n"
    print_usage && exit 1
fi

if [ -z "$sprot_db" -o -z "$pfam_db" ]; then
    echo -e "\n${red}Error: missing required database(s)${reset}\n"
    print_usage && exit 1
fi

if [ -z "$TEs_lib" ]; then
    if [ -z "$TEs_sp" ]; then
    RepeatMasker_lib_error && exit 1
    fi
fi

if [ -z "$ovres" ]; then
  echo -e "\n${red}Error: missing overlap resolution argument${reset}\n"
  print_usage && exit 1
  else
    if [[ $ovres == 'higher_score' || $ovres == 'longer_element' || $ovres == 'lower_divergence' ]]; then
      echo -e "[$(printf '%(%F %T)T\n')][INFO] Overlap resolution argument ok: $ovres"
      else
      echo -e "\n${red}Error: Wrong overlap resolution argument${reset} : $ovres\n"
      print_usage && exit 1
    fi
fi

echo "[$(printf '%(%F %T)T\n')][INFO] checking fastq files extensions.."
if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        echo "[$(printf '%(%F %T)T\n')][INFO] libraries extension format ok..."
        else 
        file_extension_error
        echo -e "\n${red}Error: check extension format :${reset} $fn \n"
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
               echo -e "\n${red}Error: check name format: ${reset} $fn\n"
               exit 1
            fi
         else
         file_extension_error
         echo -e "\n${red}Error: check extension format: ${reset} $fn\n"
         exit 1
         fi
     done
 else
    echo -e "\n${red}Error: Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi
fi

mkdir -p $outfolder/temp
cd $outfolder/temp

#TransDecoder
echo -e "${yellow} 
Identifying candidate coding regions with TransDecoder
======================================================${reset}"

echo -e "[$(printf '%(%F %T)T\n')][INFO] TransDecoder Step 1/3: extract the long open reading frames"
$TransDecoder_path/TransDecoder.LongOrfs -S -t $fasta_trme -O $outfolder/temp -S 

echo -e "[$(printf '%(%F %T)T\n')][INFO] TransDecoder Step 2/3: identifying ORFs with homology to known proteins via blast and pfam searches";
echo -e "[$(printf '%(%F %T)T\n')][INFO] runing BLAST"
$BLASTP_path -query $outfolder/temp/longest_orfs.pep  \
-db $sprot_db  -max_target_seqs 1 \
-outfmt 6 -evalue 1e-8 -num_threads $threads > $outfolder/temp/blastp.outfmt6

echo -e "[$(printf '%(%F %T)T\n')][INFO] runing HMMER"
$HMMER_path --cpu $threads --domtblout $outfolder/temp/pfam.domtblout $pfam_db $outfolder/temp/longest_orfs.pep 

echo -e "[$(printf '%(%F %T)T\n')][INFO] TransDecoder Step 3/3: predict the likely coding regions"
$TransDecoder_path/TransDecoder.Predict -t $fasta_trme --retain_pfam_hits $outfolder/temp/pfam.domtblout --retain_blastp_hits \
$outfolder/temp/blastp.outfmt6 --single_best_only -O $outfolder/temp 

#RepeatMasker
echo -e "${yellow}
Masking transcriptome with RepeatMasker
=======================================${reset}"
if [ $TEs_lib ]; then
echo -e "[$(printf '%(%F %T)T\n')][INFO] starting masking with custom library..."
$RM_path -a -e rmblast -pa $threads -nolow -xm -u -gff -dir $outfolder/temp -lib $TEs_lib $fasta_trme  
else
    if [ $TEs_sp ]; then
    echo -e "[$(printf '%(%F %T)T\n')][INFO] starting masking with $TEs_sp library..."
    $RM_path -a -e rmblast -pa $threads -nolow -xm -u -gff -dir $outfolder/temp -species $TEs_sp $fasta_trme
    fi
fi

echo -e "[$(printf '%(%F %T)T\n')][INFO] masking finished, moving files"
mkdir $outfolder/ExplorATE_inputs
mv -t $outfolder/ExplorATE_inputs blastp.outfmt6 $(basename $fasta_trme).transdecoder.gff3 $(basename $fasta_trme).out

echo -e "[$(printf '%(%F %T)T\n')][INFO] recording transcript lengths"
$seqkit fx2tab --length --name -i $fasta_trme | sort -k1,1 -k2,2n > transcripts_len.bed

echo -e "[$(printf '%(%F %T)T\n')][INFO] writing UTR regions"
awk -v OFS="\t" '($3 ~ /five_prime_UTR|three_prime_UTR/) {print $1,$4,$5,$3}' $outfolder/ExplorATE_inputs/$(basename $fasta_trme).transdecoder.gff3 \
| sort -t $'\t' -k1,1 -k2,2n | $bedtools merge -c 4 -o first > UTRs.bed
awk -F".p" '{print $1}' blastp.outfmt6 | sort | uniq > IDGenes.txt
awk 'FNR==NR { a[$1]=$1; next } ($1 in a)' IDGenes.txt UTRs.bed >UTR_genes.bed

#Overlap resolution

echo -e "${yellow}
Resolving overlaps with RM2Bed
==============================${reset}"

python3 $RM2Bed --out_dir . --out_prefix RM_ovres_RM2Bed --sort_criterion 'size' --split 'subclass' --ovlp_resolution $ovres $outfolder/ExplorATE_inputs/$(basename $fasta_trme).out $al

echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] filtering non-gene-overlapping sequences"
awk -v OFS="\t" '{print $1, $2, $3, $4":"$7":"$8 }' RM_ovres_RM2Bed_rm.bed | sort -k1,1 -k2,2n | $bedtools merge -c 4 -o first > RM_overes_merge_sort.bed
$bedtools intersect -wb -wa -a UTR_genes.bed -b RM_overes_merge_sort.bed -sorted | awk '!($5="")' > repeats_in_UTRs.txt
awk '{print $1}' repeats_in_UTRs.txt | sort | uniq > seqIDGenesWithRepinUTRs.txt
$seqkit grep $fasta_trme -v -w 0 -f seqIDGenesWithRepinUTRs.txt > NoGenes.fa

echo -e "[$(printf '%(%F %T)T\n')][INFO] creating references to apply Wicker-like rule"
grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" RM_ovres_RM2Bed_rm.bed \
| awk '{print $1, $4, $7, $8, $5}' > BED_wout_ovlp.bed
grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" $outfolder/ExplorATE_inputs/$(basename $fasta_trme).out >RMclean.out

    if [[ $split_by == 'name' ]];then
    echo -e "$[$(printf '%(%F %T)T\n')][INFO] making names-based references for Wicker-like rule"
    awk 'FNR==NR{a[$1]=$2;next}{ print $1"::"$2,$3, a[$1], $5/a[$1]*100}' transcripts_len.bed BED_wout_ovlp.bed > trLen_perRep.txt
    awk 'NR>3 {print $5"::"$10, $2}' RMclean.out | awk '{a[$1]++}{b[$1]=b[$1]+$2}END{for (i in a) { print i, 100-(b[i]/a[i])} }' > PerIdentity_ok.txt
    awk 'NR>3 {print $5, $6, $7,$5"::"$10}' RMclean.out > RMclean.bed
    fi
    if [[ $split_by == 'subclass' ]];then
    echo -e "$[$(printf '%(%F %T)T\n')][INFO] making subclass-based references for Wicker-like rule"
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
echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] applying Wicker-like rule${reset} "
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
        grep -v -w -f seqIDGenesWithRepinUTRs.txt RMclean.bed | grep -w -f fragmentsIDs.txt | awk '{print $1, $2, $3, $4":"$2$3}' | sort -k1,1 -k2,2n \
        | $bedtools merge -c 4 -o first >  targetTEs.bed
        $bedtools getfasta -fi $fasta_trme -bed targetTEs.bed -nameOnly -fo  targetTEs.fa
        grep '>' targetTEs.fa | sed 's/>//g' | awk -F':' '{print $0";"$3}' > references.csv
        grep '>' targetTEs.fa | sed 's/>//g' | awk -F':' '{print $1}' > targetTEs_tr_ids_noGenes.txt
        fi
    fi

echo -e "[$(printf '%(%F %T)T\n')][INFO] removing redundant sequences"
$seqkit rmdup -s targetTEs.fa | $seqkit seq -m $kmer -g > targetTEs_clean.fa
$seqkit grep $fasta_trme -v -w 0 -f targetTEs_tr_ids_noGenes.txt | $seqkit rmdup -s | $seqkit seq -m $kmer -g > decoys.fa
grep '>' decoys.fa | sed 's/>//g' | awk '{print $1}' > decoys.txt
echo -e "[$(printf '%(%F %T)T\n')][INFO] making references for Salmon"
cat targetTEs_clean.fa decoys.fa > trmeSalmon.fa
mv repeats_in_UTRs.txt $outfolder
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
echo "[$(printf '%(%F %T)T\n')][INFO] removing temporary files"
rm -r $outfolder/temp*


trap : 0
echo >&2 "${green}
======================================================
>         ExplorATE finished successfully            <
======================================================
"


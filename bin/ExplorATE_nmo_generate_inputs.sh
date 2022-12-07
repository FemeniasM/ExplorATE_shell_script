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

BLASTP_path="blastp"
HMMER_path="hmmscan"
RM_path="RepeatMasker"
threads=1
kmer=31
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
    ${green}-r${reset} RepeatMasker binary path (RepeatMasker default)
    ${green}-o${reset} output directory path
    ${green}-t${reset} fasta transcriptome
    ${green}-n${reset} blastp binary path (blastp default)
    ${green}-m${reset} hmmer binary path (hmmscan default)
    ${green}-d${reset} TransDecoder directory path 
    ${green}-u${reset} swissprot database
    ${green}-f${reset} Pfam database
    ${green}-i${reset} RepeatMasker -lib argument (custom library)
    ${green}-j${reset} RepeatMasker -species argument (species library)
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
while getopts ":p:k:r:o:t:n:m:d:u:f:i:j:h:" opt; do
    case $opt in

        p)
            threads="$OPTARG"
            echo "-p <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;
        r)
            RM_path=`realpath $OPTARG`
            echo "-r <RepeatMasker binary path> = $RM_path"
            ;;
        o)
            outfolder=`realpath $OPTARG`
            echo "-o <Output files Path> = $outfolder"
            ;;
        t)
            fasta_trme=`realpath $OPTARG`
            echo "-t <fasta transcriptome> = $fasta_trme"
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


# Required arguments
if [ -z "$outfolder" -o -z "$RM_path" -o -z "$TransDecoder_path" -o -z "$fasta_trme" -o -z "$threads" ]; then
    echo -e "${red}Error: missing required argument(s)${reset}"
    print_usage && exit 1
fi

if [ -z "$sprot_db" -o -z "$pfam_db" ]; then
    echo -e "${red}Error: missing required database(s)${reset}"
    print_usage && exit 1
fi

if [ -z "$TEs_lib" ]; then
    if [ -z "$TEs_sp" ]; then
    RepeatMasker_lib_error && exit 1
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
mv -t $outfolder/ExplorATE_inputs $outfolder/temp/blastp.outfmt6 $outfolder/temp/$(basename $fasta_trme).transdecoder.gff3 $(basename $fasta_trme).out
echo "[$(printf '%(%F %T)T\n')][INFO] removing temporary files"
rm -r $outfolder/temp*
echo -e "${green}[$(printf '%(%F %T)T\n')][INFO] input files are ready in directory $outfolder/ExplorATE_inputs ${reset}\n"

trap : 0
echo >&2 "${green}
======================================================
          ExplorATE finished successfully            
======================================================
"


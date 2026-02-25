#!/bin/bash

#Exit on error
set -e

vcf_input="$1"
vcf_names="$2"

genepop_dir="genepop_files_directory"

##### Create output directory #####
mkdir -p "${genepop_dir}"

if [[ ! -d "${genepop_dir}" ]]; then
    echo "ERROR: Failed to create output genepop directory: ${genepop_dir}" >&2
fi

##############################################################
#Function: spid_file
#Description: Creation of the spid file required for PGDSpider
###############################################################

spid_file(){
    cat > "spid_file.spid" << EOF
# spid-file generated: $(date '+%a %b %d %H:%M:%S %Z %Y')

# VCF Parser questions
PARSER_FORMAT=VCF
    
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=false
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false

# GENEPOP Writer questions
WRITER_FORMAT=GENEPOP

# Specify which data type should be included in the GENEPOP file  (GENEPOP can only analyze one data type per file):
GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP
# Specify the locus/locus combination you want to write to the GENEPOP file:
GENEPOP_WRITER_LOCUS_COMBINATION_QUESTION=
EOF
}

##################################################################
#Function: vcf_2_genepop
#Description: Function to convert VCF into genepop with PGDSpider
##################################################################

vcf_2_genepop(){

    ##### Parameters #####
    local vcf_input="$1"    #Txt file containing the path to each VCF (1 per line)
    local vcf_names="$2"
    local spid_file="spid_file.spid" 

    # Convert comma-separated lists to arrays
    IFS=',' read -ra vcf_array <<< "$vcf_input"
    IFS=',' read -ra name_array <<< "$vcf_names"

    ##### Process each VCF #####
    for i in "${!vcf_array[@]}"; do
        local vcf="${vcf_array[$i]}"
        local original_name="${name_array[$i]}"

        ##### Check if file exists #####
        if [[ ! -f "$vcf" ]]; then
            echo "File not found, ignored: $vcf"
            continue
        fi

        # Extract base name (handle .vcf)
        local base_name
        local regex='\(([^)]+)\)[[:space:]]*$'
        if [[ "$original_name" =~ $regex ]]; then
            #Extract content between last parentheses
            base_name="${BASH_REMATCH[1]}"
        else
            # No parentheses, use original name
            base_name=$(basename "$original_name")
        fi
        
        base_name=${base_name%.vcf}

        #Output GENEPOP file
        output_genfile="${genepop_dir}/${base_name}_genepop.txt"

        #Check if already processed
        if [[ -f "$output_genfile" ]]; then
            continue
        fi

        #Conversion VCF to GENEPOP
        PGDSpider2-cli -inputfile "$vcf" -inputformat VCF -outputfile $output_genfile -outputformat GENEPOP -spid "$spid_file"

    done < <(printf "%s\n" "$vcf_input" | tr ',' '\n')

}

##########################
# Main execution 
##########################

spid_file
vcf_2_genepop "$vcf_input" "$vcf_names"
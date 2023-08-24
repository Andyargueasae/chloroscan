# This script could make use of snakemake parameters to give cds extraction for those MAGs.

# make the directory first.
mkdir -p ${snakemake_output}
echo ${snakemake_output}

for i in $(ls ${snakemake_input});
do
    gffread -C ${snakemake_input}/$i -o $(basename $i .gff).cds 2> ${snakemake_log};
    mv $(basename $i .gff).cds ${snakemake_output}
done

# So basically, every script uses a for loop to do batch jobs for each of the files.
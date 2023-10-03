import os
import click

@click.command()

@click.option("--assembly", type=click.Path(exists=True))
@click.option("--accession", type=click.Path(exists=True))
@click.option("--fasta_file_name", type=click.STRING)

# sample_assembly = "mediaflux_BATCHES/MGYS00002105_baltic_ocean_surface/assembly/merged.megahit.contigs.fa"
# sample_plastid_accession = "MGYS00002105-plastid_acc.txt"
# fasta_file_name = "mediaflux_BATCHES/MGYS00002105_baltic_ocean_surface/plastid.fasta"

def main(assembly, accession, fasta_file_name):
    plastid_accession_id = []
    with open(accession, "r") as plastida:
        for lines in plastida:
            plastid_accession = lines.replace("\n","")
            # print(plastid_accession)
            plastid_accession_id.append(plastid_accession)
    print(len(plastid_accession_id)) 

    with open(assembly, "r") as assemblyf:
        assembly_lines = assemblyf.readlines()

    # print(assembly_lines[:20])

    # make a judgement, whether there are files already there.
    if os.path.exists(fasta_file_name):
        os.remove(fasta_file_name)
    CONTIG_IDENTIFER = ">"
    index = 0
    while index < len(assembly_lines):
        accession_id = assembly_lines[index].replace(">","").split(" ")[0]
        # print(accession_id)
        if accession_id in plastid_accession_id:
            sequence = assembly_lines[index+1]
            with open(fasta_file_name, "a") as fasta:
                fasta.write(CONTIG_IDENTIFER)
                fasta.write(accession_id)
                fasta.write("\n")
                fasta.write(sequence)
        index += 2
    return

if __name__ == "__main__":
    main()
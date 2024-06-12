import os
import click

@click.command()

@click.option("--assembly", type=click.Path(exists=True))
@click.option("--accession", type=click.Path(exists=True))
@click.option("--fasta_file_name", type=click.STRING)

def main(assembly, accession, fasta_file_name):
    plastid_accession_id = []
    with open(accession, "r") as plastida:
        for lines in plastida:
            plastid_accession = lines.replace("\n","")
            # print(plastid_accession)
            plastid_accession_id.append(plastid_accession)
    # print(len(plastid_accession_id)) 

    with open(assembly, "r") as assemblyf:
        assembly_lines = assemblyf.readlines()

    # print(assembly_lines[:20])

    # make a judgement, whether there are files already there.
    if os.path.exists(fasta_file_name):
        os.remove(fasta_file_name)
    CONTIG_IDENTIFER = ">"
    index = 0

    while index < len(assembly_lines):
        accession_id = assembly_lines[index].replace(">","").split(" ")[0].replace("\n", "")
        # print(accession_id)
        if accession_id in plastid_accession_id:
            # print("sequence classified as plastidial.\n")
            sequence = assembly_lines[index+1]
            # print(sequence)
            with open(fasta_file_name, "a") as fasta:
                fasta.write(CONTIG_IDENTIFER)
                fasta.write(accession_id)
                fasta.write("\n")
                fasta.write(sequence)
        index += 2
    
    # print(os.path.exists(fasta_file_name))
    return

if __name__ == "__main__":
    main()
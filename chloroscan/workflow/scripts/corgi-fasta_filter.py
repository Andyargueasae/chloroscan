import os
import click

CONTIG_IDENTIFER = ">"

@click.command()

@click.option("--assembly", type=click.Path(exists=True))
@click.option("--accession", type=click.Path(exists=True))
@click.option("--fasta_file_name", type=click.STRING)


def main(assembly, accession, fasta_file_name):
    plastid_accession_id = []
    with open(accession, "r") as plastida:
        for lines in plastida:
            plastid_accession = lines.strip()
            # print(plastid_accession)
            plastid_accession_id.append(plastid_accession)
    # print(len(plastid_accession_id)) 

    n_contig = len(plastid_accession_id)
    with open(assembly, "r") as assemblyf:
        assembly_lines = assemblyf.readlines()

    # print(assembly_lines[:20])

    # make a judgement, whether there are files already there.
    if os.path.exists(fasta_file_name):
        os.remove(fasta_file_name)
    
    index = 0
    n_count = 0
    # Here, adjust the code to jump across different contig accession line, and store the sequences under the desired accession ID.
    while index < len(assembly_lines):
        if n_count == n_contig:
            break
        # By default, the first line is always an accession for the first sequence.
        accession_id = assembly_lines[index].replace(">","").split(" ")[0].strip()
        # Now, find the next index id. 
        next_accession_index = get_next_accession(assembly_lines, index)
        if accession_id in plastid_accession_id:
            sequence = "".join([seq.strip() for seq in assembly_lines[index+1:index+next_accession_index]])+"\n"
            with open(fasta_file_name, "a") as fasta:
                fasta.write(CONTIG_IDENTIFER)
                fasta.write(accession_id)
                fasta.write("\n")
                fasta.write(sequence)
                n_count += 1
        index += next_accession_index
    
    # print(os.path.exists(fasta_file_name))
    return

def get_next_accession(assembly_lines, index):
    line_skip = 0
    for i in range(index+1, len(assembly_lines)):
        if CONTIG_IDENTIFER in assembly_lines[i]:
            line_skip = i - index
            return line_skip
        else: # Make sure the sequence is already at the end.
            line_skip = len(assembly_lines)-index

    return line_skip

if __name__ == "__main__":
    main()
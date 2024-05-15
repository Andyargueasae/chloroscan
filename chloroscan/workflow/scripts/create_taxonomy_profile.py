import pandas as pd
import click

@click.command()
@click.option("--cat_prediction_txt", type=click.Path(exists=True))
@click.option("--output", type=click.STRING)

def main(cat_prediction_txt, output):
    # transfer it into a pandas dataframe.
    cat_df=pd.read_csv(cat_prediction_txt, sep="\t")
    copy_df = cat_df.loc[cat_df['lineage'].notna()]
    lineage_hierarchy = copy_df['lineage']
    finest_lineage = []
    for i in lineage_hierarchy:
        indi_hierarchy = i.replace("*","").split(";")
        finest_classification = indi_hierarchy[-1]
        print(finest_classification)
        finest_lineage.append(finest_classification)
    
    copy_df['intake for krona'] = finest_lineage
    intake_dict = dict()
    for i in finest_lineage:
        if i in intake_dict.keys():
            intake_dict[i]+=1
        else:
            intake_dict[i]=1
    

# Here, only counts of contigs were used as magnitude, but that is arbitrarily trivial
# So what we need to do is try to get the real abundance of the sequences/taxonomy to work. 
    counts = intake_dict.values()
    taxon_id = intake_dict.keys()
    krona_intake = pd.DataFrame()
    krona_intake['count'] = counts
    krona_intake['taxon id'] = taxon_id
    krona_intake.to_csv(output, sep="\t", index=False)
    ##firstly, filter those without a ORF hit, retain those having it.
    return

if __name__ == "__main__":
    main()

    # Tips for changing the output for our krona tools:
    # 1. Shifting the count of taxonomy to the count of contig magnitudes + the official format of krona tools to work on.
    # 2. try to get the magnitude for each contig and then go for it. 
    # 3. Here is the problem: If you only use count of taxa found, it gives you how precise corgi did to classify each contig.
    #    But, adding the magnitude of contigs by its depth, can adjust the abundance to make it precise. 
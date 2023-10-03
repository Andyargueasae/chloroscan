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
    

    counts = intake_dict.values()
    taxon_id = intake_dict.keys()
    krona_intake = pd.DataFrame()
    krona_intake['count'] = counts
    krona_intake['taxon id'] = taxon_id
    krona_intake.to_csv(output, sep="\t", index=False)
    #firstly, filter those without a ORF hit, retain those having it.
    return

if __name__ == "__main__":
    main()
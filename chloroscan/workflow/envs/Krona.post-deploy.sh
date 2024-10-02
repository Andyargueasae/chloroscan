# Determine if ktUpdateTaxonomy.sh has to be run.
TAXONOMY_DIR=$CONDA_PREFIX/opt/krona/taxonomy

if [ -n $(ls -l $TAXONOMY_DIR | grep "taxonomy.tab") ]; then
    echo "Taxonomy data found, no need to run UpdateTaxonomy.sh."
else
    ktUpdateTaxonomy.sh
fi

echo "The Krona environment is all set."
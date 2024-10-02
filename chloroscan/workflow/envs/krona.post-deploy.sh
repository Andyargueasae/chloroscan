# Determine if ktUpdateTaxonomy.sh has to be run.
TAXONOMY_DIR_FILE=$CONDA_PREFIX/opt/krona/taxonomy/taxonomy.tab

if [ -f "$TAXONOMY_DIR_FILE" ]; then
    echo "Taxonomy data found, no need to run UpdateTaxonomy.sh."
else
    ktUpdateTaxonomy.sh
fi

echo "The Krona environment is all set."
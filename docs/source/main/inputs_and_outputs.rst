===================
Inputs and Outputs
===================

Overall ChloroScan takes the following inputs:
    - Metagenomic assembly in fasta format (required): configured by ``--Inputs-assembly``
    - Depth profile in tabular format (required ``tsv``, mutually exclusive with the next): configured by ``--Inputs-depth-profile``
    - directory storing alignment in sorted bam format (required): configured by ``--Inputs-alignment``
    - Batch name (required, just to distinguish different jobs when running in parallel to avoid overwriting): configured by ``--Inputs-batch-name``

Commonly for a **de novo** run we recommend using bam directory, but if you want to rerun your analysis (perhaps adjusting parameters) with present depth profile, you can use it.

The output directory structure of ChloroScan is as follows:
(This is a typical example showing only directories)


.. code-block:: text

    chloroscan_output
    ├── logging_info
    └── working
        ├── binny
        │   ├── bins
        │   ├── intermediary
        │   │   └── mantis_out
        │   ├── job.errs.outs
        │   ├── logs
        │   └── tmp
        ├── call_contig_depth
        ├── CAT
        ├── cds-extraction
        │   ├── cds
        │   └── faa
        ├── corgi
        ├── FragGeneScanRs
        │   └── GFFs
        ├── refined_bins
        ├── summary
        └── visualizations



Here is a brief description of each step regarding their inputs and outputs.

1. Corgi predictions
============================

The input of the first step: contig classification by Corgi requires your ``input-contigs.fasta`` only.

The program firstly classifies each contig into five RefSeq-based categories: chloroplast, mitochondria, prokaryote, eukaryote, and unknown.

Then the outputs: classification results: ``corgi-prediction.csv``, and the isolated ``plastids.fasta`` file will be generated in the ``output/working/corgi`` directory.

To configure this step:
    - ``--corgi-pthreshold``: the posterior probability threshold to classify a contig into one of the five categories. Default is 0.5.
    - ``--corgi-min-length``: minimum length of contigs to be considered. Default is 1000bp.
    - ``--corgi-save-filter``: run corgi's default filtering step to isolate plastid contigs, **not recommended** due to long runtime.
    - ``--corgi-batch-size``: batch size when inferring contigs, default is 1.

2. binning results
============================

The inputs of this step is the ``plastids.fasta`` from Corgi and the tabular depth profile with the first column the contig ids and the rest columns the average depths in each sample.

 - ``Note``: if users did not specify the depth table, ChloroScan directly incorporates the code from binny to generate depth profile. It's path is ``output/working/depth_profile.txt``.

The output directory of binny is in ``output/working/binny``. The directory ``bins`` contains the final plastid MAGs in fasta format.

Meanwhile because binny is a snakemake workflow, other intermediary results are also retained in ChloroScan.

Here is a typical example of it:

 .. code-block:: text

    intermediary
    ├── annotation_CDS_RNA_hmms_checkm.gff
    ├── annotation.filt.gff
    ├── assembly.contig_depth.txt
    ├── assembly.formatted.fa
    ├── iteration_1_contigs2clusters.tsv
    ├── iteration_2_contigs2clusters.tsv
    ├── iteration_3_contigs2clusters.tsv
    ├── mantis_out
    │   ├── consensus_annotation.tsv
    │   ├── integrated_annotation.tsv
    │   ├── Mantis.out
    │   └── output_annotation.tsv
    ├── prokka.faa
    ├── prokka.ffn
    ├── prokka.fna
    ├── prokka.fsa
    ├── prokka.gff
    ├── prokka.log
    ├── prokka.tbl
    ├── prokka.tsv
    └── prokka.txt

To simply, here are the things you need to know:
 - The input contigs and depth profile in tabular format are in ``assembly.formatted.fa`` and ``assembly.contig_depth.txt`` respectively. 
 - The predicted open reading frames from contigs are dumped in the file ``annotation_CDS_RNA_hmms_checkm.gff``.

To configure this step:
    - ``--binning-clustering-hdbscan-min-sample-range``: the range of minimum samples parameter when running HDBSCAN clustering, default is 1,5,10.
    - ``--binning-bin-quality-purity``: the minimum purity (i.e. max contamination) threshold for a bin to be retained, default is 90.
    - ``--binning-bin-quality-min-completeness``: the minimum completeness threshold for a bin to be retained, default is 50.
    - ``--binning-outputdir``: the output directory to store binny results before moving it to ChloroScan output, default is ``output/working/binny``.
    - ``--binning-universal-length-cutoff``: the length cutoff applied to filter contigs for marker gene prediction and clustering, default is 1500bp.
    - ``--binning-clustering-epsilon-range``: the epsilon range for HDBSCAN. We do not recommend changing this, as it shows minimal effect, default is 0.250,0.000.

For further information on HDBSCAN, see the sklearn docs: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html.

3. CAT
============================
The input is the ``assembly.formatted.fasta`` from the intermediary directory. It stores putative plastid contigs longer than 500bp. 

Overall the only useful information is in ``out.CAT.contig2classification.txt``, it is a tabular text file storing the classification outcome, putative lineage and the corresponding scores for each taxon for each contig.

It looks like this:

.. code-block:: text

        # contig	classification	reason	lineage	lineage scores
    k141_1000298	no taxid assigned	no hits to database
    k141_1000835	taxid assigned	based on 2/2 ORFs	1;131567;2;1224;28211;54526	1.00;1.00;1.00;1.00;1.00;1.00
    k141_1001040	taxid assigned	based on 2/2 ORFs	1;131567;2;1224;28211;766;1699067;2026788	1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00
    k141_1002229	taxid assigned	based on 2/3 ORFs	1;131567;2	1.00;1.00;1.00
    ...

To configure this step:
 - ``--cat-database``: path to CAT database (the ``db`` directory), we refer to the directory structure of the most recent CAT db release.
 - ``--cat-taxonomy``: path to the taxonomy file (the ``taxonomy`` / ``tax`` directory).

4. Summary files
============================

This step is to summarize the metadata including binning and CAT results. It takes inputs from previous jobs to summarize:
 - basic information of contigs like GC contents, length, average depths in each sample.
 - binning results: bin id for each contigs. Unbinned contigs are given as NaN.
 - CAT results: taxonomic lineage for each contig.
 - raw sequence.

The output ``summary_table.tsv`` is a tabular text file storing all the above information for each contig.

5. Refinement module
============================

Though in binning step, binny already takes serious efforts in maintaining the purity of bins. We still found it vulnerable to fragmental contigs by appending them into bins.

So in order to give users a notification, we here reports the contigs with no marker genes predicted by binny and has ambiguous taxonomic classification by CAT (i.e. not eukaryotic or unclassified).

Another way to assist users is to run anvi'o interactive mode to see if contigs form uniform clusters in the phylogram. Further guides: https://astrobiomike.github.io/genomics/metagen_anvio.

The input is the ``summary_table.tsv`` from the previous step. The output directory is ``output/working/refined_bins``. We also summarized the information into ``output/working/refinement_contig_summary.txt`` to report which contig in which bin has suspicious identity.

6. cds-extraction
============================

For each bin we predict their coding nucleotide sequences and proteins via FragGeneScanRs.

The output is in ``output/working/cds-extraction``. The directory ``cds`` contains nucleotide sequences in fna format while ``faa`` contains protein sequences in faa format.

We configure FragGeneScanRs with the following parameters:
 - ``-t illumina_5``: to use the Illumina error model with 0.5% error rate, which is suitable for metagenomic data.

7. Visualizations
============================
This step generates some visualizations for users to have a glance at their data.

For all MAGs in the sample, we plot their GC contents against their log-transformed pooled average depths. Pooled means summing up the depths in each sample for each contig. The output is in ``output/working/visualizations``. The file ``GC_vs_depth.png`` is a scatter plot showing the distribution of all MAGs in the sample.
We also generate a violin plot showing the distribution of pooled average depths for all MAGs in the sample. The output is in ``output/working/visualizations``. The file ``depth_distribution.png`` is a violin plot showing the distribution of pooled average depths for all MAGs in the sample.
For each bin, we plot their contig-level taxonomy classification in pie charts. The are named as ``{Batch_name}_{bin_id}_taxonomy_composition.png``.

8. Some notes.
============================
Parameters used for this analysis are stored in ``output/arguments.txt`` for users to check.
"Plastid" contigs count out of total contigs is reported in ``corgi.summary.txt``.

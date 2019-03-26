# RSAT analysis

This protocol was modified based on [Castro-Mpndragon et al, 2019. RSAT::Plants: Motif Discovery in ChIP-seq Peaks of Plant Genomes](https://link.springer.com/protocol/10.1007%2F978-1-4939-6396-6_19), presenting a step-by-step instruction for the task of discovering and annotating DNA motifs in MNase-sensitive footprints (MSFs) for cotton using the [RSAT::plants](http://rsat.eead.csic.es/plants/) web server. Because the cotton genome used is not supported (*G. raimondii* genome was there, but probably different version from mine, I got errors using the "Sequence tools"), input sequence files need to be prepared ahead. 

Research questions and planned analyses:

1. What motifs are enriched in MSF regions of each cotton genome studied? 
    * 1) ab initio discovery of MSF sequences using whole reference genome as control, 2) compare discovered motfis with plant motif databases, and 3) clustering the discovered motifs to obtain a nonredundant collection. Optionally, repeat for 10 sets of random sequences as negative control.
2. What motifs are differentially enriched between cotton genomes (A2, D5, AD1, A2xD5)?
    * perform motif clustering for motifs discovered from all genomes

Analytic questions:

1. Besides motif discovery of MSF against whole genome, what about MRFs against whole genome, and MSFs against MRFs? 
2. RSAT online submission allows the “Number of motifs per algorithm” is set to up to 10. Is it possible to use larger number online?
3. It has been over 24 hrs since submission of my smallest input (1% MSFs of 800Mb genome), still running. My largest input will be about maize size. How long does it take to analyze the maize MOA-seq data with local installation to find 500 motfis?
 
    

## Input file preparation

Using the MSF coordinates (BED) and refernce genome to retrieve footprint sequences. The default setting in RSAT is to use masked reference sequence with repeats as N's. The reference will be used as control genome.

### R srcipt

     cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/isegv1.3.2
     
     # get MSFs BED
     grep '20,20,255' DcD_bc6.0.Fus.bed >DcD_bc6.0.Fus.MSF.bed
     
     # get fasta
     module load bedtools2
     bedtools getfasta -fi ../refGenomes/Dgenome2_13.fasta -bed DcD_bc6.0.Fus.MSF.bed -fo DcD_bc6.0.Fus.MSF.fasta
     
     # Or use repeat masked reference, N's 
     # https://www.cottongen.org/species/Gossypium_raimondii/jgi_genome_221
     bedtools getfasta -fi ../refGenomes/G.raimondii_JGI_221_v2.0.assembly-hardmasked.fasta -bed DcD_bc6.0.Fus.MSF.bed -fo DcD_bc6.0.Fus.MSF.ms.fasta
           

## Step-by-step instruction for online submission

### Ab Initio Motif Discovery by peak-motifs section
1. Open a connection to the RSAT::Plants server. It can be reached at [http://plants.rsat.eu](http://plants.rsat.eu) and also at [http://floresta.eead.csic.es/rsat](http://floresta.eead.csic.es/rsat).
2. In the left menu of the RSAT server, click on “NGS - ChIP-seq” and next "peak-motifs" to open the pipeline form.
3. Add a title for this job, such as ```D_BC6_ms```. Upload [DcD_bc6.0.Fus.MSF.ms.fasta]() as **Peak sequences**, and upload [G.raimondii_JGI_221_v2.0.assembly-hardmasked.fasta]() as **Control sequences**.
4. Click on "Reduce peak sequences" and leave both fields blank ("number of top sequences to retain" and "cut peak sequences") to avoid having the sequences clipped.
5. Click on "Motif discovery parameters". Select "Discover over-represented words" [oligo - analysis] and "Discover over-represented spaced word pairs" [dyad - analysis] typically bound by dimeric transcription factors. The other two algorithms - "Discover words with a positional bias" [position - analysis] and "Discover words with local over-representation" [local-word-analysis] for detecting k-mers with local overrepresentation, will not run when control set is provided. If without control set, the two default selected [oligo - analysis] and [position - analysis] was considered to be complementary criteria (over- representation and positional distribution relative to peak centers) offering a good trade-off between computing time and sensitivity; other algorithms can be selected are andMore selections may search more thoroughly but often led to redundancy in the result.
6. Under "Motifs discovery" activate oligomer lengths 6&7 (not motif length, like k-mer length used for building position-specific scoring matrices) and set Markov order is to automatic. **The website allows the “Number of motifs per algorithm” is set to up to 10. Is it possible to define larger number?**
7. Click on "Compare discovered motifs with databases" and
select appropriate databases which will be used to annotate any found motifs: "footprintDB-plants","RSAT non-redundant plants","JASPAR core nonredundant plants". **About "Add your own motif database" and "add known reference motifs for this experiment", will it be useful to provide a species-specific motif set predicted by plantTFDB (http://planttfdb.cbi.pku.edu.cn/download.php#bind_motif)?**
8. Click on "Locate motifs and export predicted sites as custon USCS tracks", check the option "Search putative binding sites in the peak sequences" and "Peak coordinates specified in fasta headers in bedtools getfasta format (also for retrieve-seq-bedoutput)".
9. Select output type (display or email) and press “GO”.
10. Once the job is complete click on the link [Download all results (peak -motifs_archive.zip)] to save the results on your computer, and click on the link [Download all matrices (transfac format)] to download the full set of discovered motifs to a local file named ```D_BC6_ms.motifs.tf```. This file contains all motifs in the form of position- weight matrices (PWMs) in TRANSFAC format. This file will be used for matrix clustering next
11. In the Sequence composition (test sequences) section,
right-click on the link “[coordinates: UCSC BED track]”
(right panel) and save the BED file as peak-motifs_test_seqcoord.bed . This fi le contains the peaks used for the peakmotifs analysis.
12. At the bottom of the Web page, look for section Motif
locations (sites) , then Predicted sites on test peaks (all
motifs) . Right-click on the “[bed]” link to download the
corresponding fi le peak-motifs_all_motifs_seqcoord.bed . This file can be loaded in a genome browser such as IGV.

The execution time was .....

#### Sequence composition
* Sequence length distribution showed that most MSFs detected have a mean length around 77 bp.
* Nucleotide composition showed an enrichment of A and T around segment center
* Dinucleotide composition showed enrichment of AA/TT followed by AT/TA.

#### Motif discovery
check significance score

#### Comparison with known motifs
check correlations, positional distribution and number of binding sites per sequence


### Motif clustering

1. Open the toolbox "Matrix tools" and click "matrix-clustering", next click "Mandatory inputs".
2. On the title box give a title to the analysis ```D6ms_cluster```.
3. Upload the motif file obtained from peak-motifs and select the TRANSFAC format. Up to three input files can be uoloaded by "Optional inputs".
4. Click "Advanced options" to fine-tune the thresholds. Set w to 5, cor to 0.75, and Ncor to 0.55, leaving others as default (e.g., Ncor (Normalized Pearson Correlation) as a Metric to build the trees and average as the Agglomeration rule).
5. Click "Ran analysis", select "Heatmap", and leave "Export Radial Tree (motif browser)" and "Negative control: the input motifs columns are randomly permuted" unselected, select email output and fill up address and click GO.
6. Optionally, repeat step 1-5 except selecting "Negative control" to estimate how many clusters can be found when motifs are permuted, 100 times followed by boxplot could be a good test.

The execution time was .....

#### Interpretation

### Negative control with random genomic regions

The idea is to first generate a set of random sequences with the same lengths as the input MSFs from the same reference genome, and then submit this random set as input to "peak-motifs" with the same options. The motifs returned (much better no enrichment)from radom sets are expected to show much lower significance score or have low complexity. Next, through matrix clustering, these control motifs are expected NOT to be clustered with discovered motifs.

 
Server Command:
 
    # motif discovery
    peak-motifs  -v 1 -title D_BC6_ms -i $RSAT/public_html/tmp/apache/2019/03/04/peak-motifs.2019-03-04.222721_2019-03-04.222721_8HcKCq/peak-motifspeak_seq -ctrl $RSAT/public_html/tmp/apache/2019/03/04/peak-motifs.2019-03-04.222721_2019-03-04.222721_8HcKCq/peak-motifscontrol_seq.gz -markov auto -disco oligos,dyads,positions -nmotifs 10 -minol 6 -maxol 7 -no_merge_lengths -2str -origin center -motif_db RSAT_nonredundant_plants tf $RSAT/public_html/motif_databases/RSAT_nonredundant_plants_2017.tf -motif_db footprintDB-plants tf $RSAT/public_html/motif_databases/footprintDB/footprintDB.plants.motif.tf -motif_db jaspar_core_nonredundant_plants tf $RSAT/public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_plants_non-redundant_pfms_transfac.tf -scan_markov 1 -source getfasta -task purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan -prefix peak-motifs -noov -img_format png -outdir $RSAT/public_html/tmp/apache/2019/03/04/peak-motifs.2019-03-04.222721_2019-03-04.222721_8HcKCq
    # clustering
    $RSAT/perl-scripts/matrix-clustering  -v 1 -max_matrices 300 -matrix d6ms_motifs $RSAT/public_html/tmp/apache/2019/03/11/matrix-clustering_2019-03-11.180237_qSyrTu/matrix-clustering_query_matrices.transfac transfac -hclust_method average -calc sum -title 'd5' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.75 -lth Ncor 0.55 -quick -label_in_tree name -return json,heatmap  -o $RSAT/public_html/tmp/apache/2019/03/11/matrix-clustering_2019-03-11.180237_qSyrTu/matrix-clustering 2> $RSAT/public_html/tmp/apache/2019/03/11/matrix-clustering_2019-03-11.180237_qSyrTu/matrix-clustering_err.txt


### Result interpretation

The execution time was .....

#### Sequence composition
* Sequence length distribution showed that most MSFs detected have a length around ??? bp.
* Nucleotide composition showed an enrichment of ?? around segment center
* Dinucleotide composition showed

#### Motif discovery
check significance score

#### Comparison with known motifs
check correlations, positional distribution and number of binding sites per sequence

### Negative Control using Random Groups of Genes

1. Repeat above steps 10 times, by using random promoters as input.
2. The top motifs found by oligo-analysis and
dyad-analysis are used as negative control, in comparison with motifs found from *trans* affected genes.

## Local RSAT Installation

[http://download.rsat.eu](http://download.rsat.eu)

    cd /Users/Jing/rsat
    source RSAT_config.bashrc
    echo $RSAT
    
Something not right with `matrix-clustering`, and I had to download and reinstall it from <https://github.com/jaimicore/matrix_clustering>, otherwise error with the TFBMclust package.

    
### Test Run 1 ### 

Default discovery without control, 5 motifs defined. About 10 minutes run, 5 motifs each reported by algorithms `oligos_6nt_mkv4`,`positions_6nt`,`oligos_7nt_mkv5`,`positions_7nt`.

    # motif discovery, no control, 5 motifs. Pretty quick. Oligos and Position analyses by default.
    peak-motifs -v 3 -title D_BC6_ms -i ~/Downloads/DcD_bc6.0.Fus.MSF.ms.fasta -outdir ~/Downloads/D6ms/ -prefix D6ms -motif_db RSAT_nonredundant_plants tf $RSAT/public_html/motif_databases/RSAT_nonredundant_plants_2017.tf -motif_db footprintDB-plants tf $RSAT/public_html/motif_databases/footprintDB/footprintDB.plants.motif.tf -motif_db jaspar_core_nonredundant_plants tf $RSAT/public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_plants_non-redundant_pfms_transfac.tf -img_format png -nmotifs 5 -markov auto -minol 6 -maxol 7 -no_merge_lengths -2str -origin center -scan_markov 1 -noov 2>~/Downloads/D6ms/log.txt &
    
    # matrix-clustering. Use Firefox to view output.
    matrix-clustering  -v 1 -max_matrices 300 -matrix 'd6ms_motifs' ~/Downloads/D6ms/results/discovered_motifs/D6ms_motifs_discovered.tf transfac -hclust_method average -calc sum -title ''d5'' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.75 -lth Ncor 0.55 -quick -label_in_tree name -return json,heatmap  -o ~/Downloads/D6ms/clustering/D6ms 

[(Castro-Mondragon et al 2017)](https://academic.oup.com/nar/article/45/13/e119/3862068#119463895) `-radial_tree_only` will list only a single rooted tree in summary html. `-rand` When this option is selected, the columns of the input motifs are randomly permuted (conserving thus the Information Content), the new motifs are used as input for the pairwise-comparison and clustering. The resulted number of clusters are probably different every time, but povide a rough estimation of the number of clusters that can be randomly generated. 


### Test Run 2 ### 

Default discovery without control, 100 motifs defined. About ?? minutes run time, 5 motifs each reported by algorithms `oligos_6nt_mkv4`,`positions_6nt`,`oligos_7nt_mkv5`,`positions_7nt`.

    # motif discovery: job 20615
    peak-motifs -v 3 -title D_BC6_ms -i ~/Downloads/DcD_bc6.0.Fus.MSF.ms.fasta -outdir ~/Downloads/D6ms_100/ -prefix D6ms -motif_db RSAT_nonredundant_plants tf $RSAT/public_html/motif_databases/RSAT_nonredundant_plants_2017.tf -motif_db footprintDB-plants tf $RSAT/public_html/motif_databases/footprintDB/footprintDB.plants.motif.tf -motif_db jaspar_core_nonredundant_plants tf $RSAT/public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_plants_non-redundant_pfms_transfac.tf -img_format png -nmotifs 100 -markov auto -minol 6 -maxol 7 -no_merge_lengths -2str -origin center -scan_markov 1 -noov 2>~/Downloads/D6ms_100/log.txt &
    ---book
    # matrix-clustering. Use Firefox to view output.
    matrix-clustering  -v 1 -max_matrices 300 -matrix 'd6ms_motifs' ~/Downloads/D6ms/results/discovered_motifs/D6ms_motifs_discovered.tf transfac -hclust_method average -calc sum -title ''d5'' -metric_build_tree 'Ncor' -lth w 5 -lth cor 0.75 -lth Ncor 0.55 -quick -label_in_tree name -return json,heatmap  -o ~/Downloads/D6ms/clustering/D6ms
    
    
## Reference

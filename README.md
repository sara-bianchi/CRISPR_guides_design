# CRISPR_guides_design
Automated pipeline for customizable CRISPRa and CRISPRi guides design

# environment set up and installation
All the installations are provided in the docker published on our repository as hedgelab/crispr:image1.
To use this pipeline is highly recommended to use it through the docker image to ensure reproducibility, alernatively it is possible to manually repeat the installation on your computer following the steps in the Dockerfile (optimised for a linux bash shell).

If you do not have Docker installed on your computer yet, you can follow the instructions in the following link: https://docs.docker.com/get-started/get-docker/.
Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run the analysis on an already indexed genome, 10 GB of RAM will be enough. Else, if you will need to run the indexing step, at least 30 GB are needed. To set the proper RAM limits for your analysis, some versions of Docker Desktop allow you to directly edit these parameters in the settings, while others require the use of an external settings file. In the latter case, create a .wslconfig file in the User/user_name folder of your computer. You can follow the example template provided. After making these changes, restart the Docker engine to apply them effectively.

Once Docker is installed and properly configured, activate the Docker engine directly opening the Docker Desktop application or by running on your terminal: start docker

# input files preparation
To run the analysis it is required to select a working directory on your computer. All the input files must be located in this directory, either directly or organized in subfolders, and the output will be written in this directory as well. This is important as the docker container which will be created needs to be connected with a specific folder which will be shared between your computer and the container itself. Since all the data within this folder will be accessed through the container it is important to use the prefix /home/shared_folder/ at the beginning of each path required. This will be the path under which your folder will be located in your container. For example, if your folder is D:/home/data/CRISPR_design/index, and my working directory is CRISPR_design, to indicate the index folder I should write /home/shared_folder/index.

Within this folder it is necessary to locate:
- the settings.xlsx file: an excel file with all the information needed to personalize your design. This file needs to be saved directly in the shared folder (not in subdirectories) to allow direct recognition from the pipeline as /home/shared_folder/settigs.xlsx. This table is made up of three columns:
  - 1^st columns: the setting parameter name, which must remain unchanged
  - 2^nd column: a short description of the paramter and instructions to compile it
  - 3^rd and further columns: the settings value assigned to each paramter. This are the columns which the user is required to edit. If more than one single value is needed, further values must be added in subsequent columns. For example the target_window, which indicates the range of nucleotides in which the guide is designed, requires an initial and a final coordinate, which must be added in the first two columns.
  It is possible to find two examples of the structure of this file.
- genome information: files needed for the indexing. If the Bowtie index has already been created, the genome folder is needed, else the FASTA file is necessary. For hg38 it is possible to use the https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz, and for hg19 https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz. At the moment the pipeline can only support these genomes but it will be updated with other versions and species.
- chromatin accessibility information: to evaluate if the region in which the guide will be alligned is accessible it is important to compare the alignment results with a chromatin accessibility ATAC-Seq track in BedGraph format. If your file is provided in another format like a BigWig it is necessary to convert it, for example with the UCSC bigWigToBedGraph in.bigWig out.bedGraph command (https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-bigwigtobedgraph) This track must be chosen dependently on the cell type you are interested in. For example, for hESC and subsequent steps of mesoderm to cardiac differentiation it is possible to use the GSE106689 accession number.

# running the pipeline
Once your input files are correctly set up it is possible to easily run the pipeline in two steps: open your terminal and write

docker run -d -v  [Path/To/Your/Folder]:/home/shared_folder --name [container_name] hedgelab/crispr:image1

with this first step the connection between your folder and the container is created. If you have ot previously downloaded the docker image required it will be downloaded with this first command. The only parts that you need to edit in this command are the ones in squared brackets with the path of your working directory and a unique identifier name that you want to assign to your container. Once it is done you can proceed with the next step:

docker exec -it [container_name] Rscript /home/shared_folder/Rscript.R

This command will run the Rscript.R in the container of the provided name (which must be the same used previously), and the output will be written in the Output_yymmdd directory.

# rpipeline description and output interpretation
The analysis consists of three main steps plus an optional one:
- optional step: if the needs to be generated the pipeline will create it from the provided fasta file, and locate it in the output folder into the directory called genome_hg38 or genome_hg19 dependening on the genome type used.
- 1^st step: the CrisprVerse pipeline https://github.com/crisprVerse is used for the identification of the potential guides in the selected window, with the selected length, with a PAM site dependent on the nuclease on the selected indexed genome.
- 2^nd step: guides are realigned on the reference genome in order to detect potentially off-targets allignment. All the allignments from 0 to 3 mismatches are taken nto consideration. Guides with one or more alignment with 0 or 1 mismatch out of the target window are automatically discarded. Guides with 2 or 3 mismatches in off-target alignments are scored based on an aligment score caluclated as follows:

  alignment_score = prom_3 * 2 + (al_3 - prom_3) + prom_2 * 3 + (al_2 - prom_2) * 2

  where al_3 and al_2 are the total numebr of alignment with either 2 or 3 mismatches, while the prom_2 and prom_3 values are, for each category of mismatches, the number of alignments within the region surrounding another promoter. This is calculated considering a window between -2000 and 500 bp from the TSS.
- 3^rd step: for each guide it is calculated the area under the curve of the ATAC track, meaning that for each of the x bases in the gRNA it is calculated the sum of the correpsonding ATAC scores. If this value is 0, when the gRNA is not in correspondence of an open chromatin region, for each guide it is calculated the minimal distance from another ATAC peak both at the 3' or 5'.

The output for each gene will contain:
- a spacer_ID, unique for each gRNA
- genomic information on the location of each guide, including chr, the start and end positions, the strand, the gene name, the ENSEMBL_ID for that gene, and if more isoforms are present it is insteresting to look also at the transcript_ID, promoter_ID and TSS_ID.
- the position of some functional elements, like the PAM, the TSS_position and distance_to_TSS of the guide
- sequence information with the spacer and protospacer sequences
- alignment score information with the total alignment with 2 or 3 mismatches (alignment_2_mismatches, alignment_3_mistmatches), and the list of the other off targets promoters (other_promoters_2, other_promoters_3) and the alignment_score
- ATAC score information, with the atac_score and the distance of the closer ATAC peak at 3' (dist_3) and at 5' (dist_5). These values are not calculated and they remain to 0 if the atac_score is not 0.

What is suggested for the usage of thsi tool is to first look at the alignment score (which is also used to rank the final output), and only subsequently to look at the ATAC score or even at the TSS_distance if needed.

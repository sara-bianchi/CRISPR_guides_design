#libraries
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(Rbowtie)
library(crisprBwa)
library(dplyr)
library(readxl)

current_time = format(Sys.time(), "%y%m%d")
dir.create(paste0("/home/shared_folder/Output_", current_time))

settings = as.data.frame(read_excel("/home/shared_folder/settings.xlsx", col_names = FALSE))
rownames(settings) = settings$...1
settings = settings %>% dplyr::select(-c("...1", "...2"))
settings = as.data.frame(t(settings))

if(length(unique(na.omit(settings$modality))) != 1 | !as.character(settings$modality[1]) %in% c("CRISPRa", "CRISPRi")){
  stop("invalid argument modality")
} else{modality = as.character(settings$modality[1])}

if(length(unique(na.omit(settings$nuclease))) != 1 | !as.character(settings$nuclease[1]) %in% c("SpCas9", "SaCas9")){
  stop("invalid argument nuclease")
} else{nuclease = as.character(settings$nuclease[1])}

if(length(unique(na.omit(settings$target_window))) != 2 | as.integer(settings$target_window[1]) == "NA" | as.integer(settings$target_window[2]) == "NA"){
  stop("invalid argument target_window")
} else{target_window = sort(c(as.integer(settings$target_window[1]), as.integer(settings$target_window[2])))}

if(length(unique(na.omit(settings$spacer_len))) != 1 | as.integer(settings$spacer_len[1]) == "NA"){
  stop("invalid argument spacer_len")
} else{spacer_len = as.integer(settings$spacer_len[1])}

if(length(unique(na.omit(settings$indexing))) != 1 | as.logical(settings$indexing[1]) == "NA"){
  stop("invalid argument indexing")
} else{alignment = as.logical(settings$indexing[1])}

if(length(unique(na.omit(settings$acessibility))) != 1 | file.exists(settings$acessibility[1]) == FALSE){
  stop("invalid argument for acessibility")
} else{acessibility = as.character(settings$acessibility[1])}

if(alignment == TRUE){
  if(length(unique(na.omit(settings$fasta))) != 1 | file.exists(settings$fasta[1]) == FALSE){
    stop("invalid argument for fasta")
  } else{fasta = as.character(settings$fasta[1])}
} else {
  if(length(unique(na.omit(settings$genome_dir))) != 1 | dir.exists(settings$genome_dir[1]) == FALSE){
    stop("invalid argument for genome_dir")
  } else{genome_dir = as.character(settings$genome_dir[1])}
}

if(length(unique(na.omit(settings$genome))) != 1 | !settings$genome[1] %in% c("hg38", "hg19")){
  stop("invalid argument for genome")
} else{genome = as.character(settings$genome[1])
      if(genome == "hg38"){library(BSgenome.Hsapiens.UCSC.hg38}
      if(genome == "hg19"){BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
                          library(BSgenome.Hsapiens.UCSC.hg19)}

#alignment
if(alignment == TRUE){
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/genome_", genome))
  bowtie_build(fasta,
               outdir = paste0("/home/shared_folder/Output_", current_time, "/genome_", genome, "/"),
               force = TRUE,
               prefix = genome)
  genome_dir = paste0("/home/shared_folder/Output_", current_time, "/genome_", genome)
}

#spacer identification
genes = na.omit(settings$genes)

data("tss_human", package = "crisprDesignData")

ATAC <- read.delim(acessibility, header=FALSE)
colnames(ATAC) = c("chr", "start", "end", "value")
ATAC = ATAC %>% dplyr::mutate(chr = paste0("chr", chr)) %>% dplyr::mutate(peak = 1:length(rownames(ATAC)))


for(gene in genes){
  target_region = queryTss(tss_human,
                          queryColumn = "gene_symbol",
                          queryValue = gene,
                          tss_window = target_window)
  if(genome == "hg38"){
    if(nuclease == "SpCas9"){
      data(SpCas9, package = "crisprBase")
      gs = findSpacers(target_region,
                       crisprNuclease = SpCas9,
                       spacer_len = spacer_len,
                       bsgenome = BSgenome.Hsapiens.UCSC.hg38)
    } else if(nuclease == "SaCas9"){
      data(SaCas9, package = "crisprBase")
      gs = findSpacers(target_region,
                       crisprNuclease = SaCas9,
                       spacer_len = spacer_len,
                       bsgenome = BSgenome.Hsapiens.UCSC.hg38)
    }
    gs = addTssAnnotation(gs,
                          tssObject = tss_human,
                          tss_window = target_window)
    results = as.data.frame(gs@elementMetadata@listData)
    gs2 = addSpacerAlignments(gs,
                              aligner = "bowtie",
                              aligner_index = paste0(genome_dir, "/", genome),
                              bsgenome = BSgenome.Hsapiens.UCSC.hg38,
                              n_mismatches = 3,
                              tssObject = tss_human,
                              tss_window = c(-2000, 500))
  } else if(genome == "hg19"){
    if(nuclease == "SpCas9"){
      data(SpCas9, package = "crisprBase")
      gs = findSpacers(target_region,
                       crisprNuclease = SpCas9,
                       spacer_len = spacer_len,
                       bsgenome = BSgenome.Hsapiens.UCSC.hg19)
    } else if(nuclease == "SaCas9"){
      data(SaCas9, package = "crisprBase")
      gs = findSpacers(target_region,
                       crisprNuclease = SaCas9,
                       spacer_len = spacer_len,
                       bsgenome = BSgenome.Hsapiens.UCSC.hg319)
    }
    gs = addTssAnnotation(gs,
                          tssObject = tss_human,
                          tss_window = target_window)
    results = as.data.frame(gs@elementMetadata@listData)
    gs2 = addSpacerAlignments(gs,
                              aligner = "bowtie",
                              aligner_index = paste0(genome_dir, "/", genome),
                              bsgenome = BSgenome.Hsapiens.UCSC.hg19,
                              n_mismatches = 3,
                              tssObject = tss_human,
                              tss_window = c(-2000, 500))
  }
  
  alignments = as.data.frame(gs2@elementMetadata@listData[["alignments"]])
  spacers = alignments %>% dplyr::select(c("group_name", "spacer")) %>% dplyr::distinct()
  colnames(spacers) = c("tssAnnotation.group_name", "spacer")
  results = results %>% dplyr::left_join(spacers, by = "tssAnnotation.group_name")
  spacers_to_del = c()
  alignments$promoters[is.na(alignments$promoters)] = 0
  alignments_0_1 = alignments %>% dplyr::filter(n_mismatches != 2 & n_mismatches != 3)
  for(i in 1:length(rownames(alignments_0_1))){if(alignments_0_1$promoters[i] != gene){spacers_to_del = append(spacers_to_del, alignments_0_1$group_name[i])}}
  spacers_to_del = unique(spacers_to_del)
  results_filtered = results %>% dplyr::filter(!tssAnnotation.group_name %in% spacers_to_del)
  
  if(length(rownames(results_filtered)) != 0){
    alignments_filtered = alignments %>% dplyr::filter(group_name %in% results_filtered$tssAnnotation.group_name)
    other_alignments = alignments_filtered %>% dplyr::filter(promoters != gene)
    alignments_3 = c()
    alignments_2 = c()
    promoters_3 = c()
    promoters_2 = c()
    scores = c()
    for(i in results_filtered$tssAnnotation.group_name){al_i = other_alignments %>% dplyr::filter(group_name == i);
    al_3 = length(rownames(al_i %>% dplyr::filter(n_mismatches == 3)));
    al_2 = length(rownames(al_i %>% dplyr::filter(n_mismatches == 2)));
    alignments_3 = append(alignments_3, al_3);
    alignments_2 = append(alignments_2, al_2);
    other_prom_3 = al_i %>% dplyr::filter(n_mismatches == 3) %>% dplyr::filter(promoters != 0)
    other_prom_2 = al_i %>% dplyr::filter(n_mismatches == 2) %>% dplyr::filter(promoters != 0)
    prom_3 = length(rownames(other_prom_3));
    prom_2 = length(rownames(other_prom_2));
    proms_3 = ""
    proms_2 = ""
    for(j in unique(other_prom_3$promoters)){proms_3 = paste(proms_3, j)}
    promoters_3 = append(promoters_3, proms_3);
    for(j in unique(other_prom_2$promoters)){proms_2 = paste(proms_2, j)}
    promoters_2 = append(promoters_2, proms_2);
    score_i = prom_3 * 2 + (al_3 - prom_3) + prom_2 * 3 + (al_2 - prom_2) * 2
    scores = append(scores, score_i)}
    results_filtered$alignments_3_mismatches = alignments_3
    results_filtered$alignments_2_mismatches = alignments_2
    results_filtered$other_promoters_3 = promoters_3
    results_filtered$other_promoters_2 = promoters_2
    results_filtered$other_promoters_3[is.na(results_filtered$other_promoters_3)] = 0
    results_filtered$other_promoters_2[is.na(results_filtered$other_promoters_2)] = 0
    results_filtered$scores = scores
    start = c()
    end = c()
    for(i in 1:length(rownames(results_filtered))){
      if(results_filtered$tssAnnotation.strand[i] == "+"){end = append(end, results_filtered$pam_site[i] - 1); start = append(start, results_filtered$pam_site[i] - 1 - spacer_len)}
      else if(results_filtered$tssAnnotation.strand[i] == "-"){start = append(start, results_filtered$pam_site[i] + 1); end = append(end, results_filtered$pam_site[i] + 1 + spacer_len)}
    }
    results_filtered$start = start
    results_filtered$end = end
    results_final = results_filtered %>% select("tssAnnotation.group_name",  "tssAnnotation.chr", "start", "end", "tssAnnotation.gene_symbol", "tssAnnotation.gene_id", "tssAnnotation.tx_id", "tssAnnotation.promoter", "tssAnnotation.tss_id", "tssAnnotation.tss_pos", "tssAnnotation.dist_to_tss", "pam", "spacer", "protospacer", "alignments_2_mismatches", "alignments_3_mismatches", "other_promoters_2", "other_promoters_3", "scores")
    colnames(results_final) = c("spacer_ID", "chr", "start", "end", "gene", "ENSEMBL_ID", "transcript_ID", "promoter_ID", "TSS_ID", "TSS_position", "distance_to_TSS", "PAM", "spacer", "protospacer", "alignments_2_mismatches", "alignments_3_mismatches", "other_promoters_2", "other_promoters_3", "alignment_scores")
    
    atac = ATAC %>% dplyr::filter(chr == results_final$chr[1])
    atac_score = c()
    dist_3 = c()
    dist_5 = c()
    for(i in 1:length(rownames(results_final))){peaks = c(); value = 0
    for(j in results_final$start[i]:results_final$end[i]){
      peak = atac %>% dplyr::filter(start <= j & end > j); if(length(rownames(peak)) == 1){peaks = append(peaks, peak$peak[1]); value = value + peak$value[1]}
    }
    peaks = unique(peaks)
    if(length(peaks != 0)){atac_score = append(atac_score, value); dist_3 = append(dist_3, 0); dist_5 = append(dist_5, 0)} else{atac_score = append(atac_score, 0)
    atac_magg = atac %>% dplyr::filter(start > results_final$end[i]); if(length(rownames(atac_magg)) != 0){atac_magg = head(atac_magg$start, 1); n = atac_magg - results_final$end[i]; dist_3 = append(dist_3, n)} else{dist_3 = append(dist_3, 0)};
    atac_min = atac %>% dplyr::filter(end < results_final$start[i]); if(length(rownames(atac_min)) != 0){atac_min = tail(atac_min$end, 1); n = results_final$start[i] - atac_min; dist_5 = append(dist_5, n)} else{dist_5 = append(dist_5, 0)};
    }
    }
    results_final$atac_score = atac_score
    results_final$dist_5_start = dist_5
    results_final$dist_3_end = dist_3
    results_final = results_final %>% dplyr::arrange(alignment_scores)
    write.csv(results_final, paste0("/home/shared_folder/Output_", current_time, "/", gene, "_", modality, "_results.csv"), quote = F, row.names = T)
  } else {
    print(paste0("no guides found for ", gene))
  }
}

                             


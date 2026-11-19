library(pubmed.mineR)
library(rentrez)
library(xml2)

#setwd("/home/adminsb/disk2/meta-analysis/rnafusion/pubmed_miner")
#out <- "test1"
#dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd("/home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/")
out="/home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/litData_9oct"
dir.create(out)

for (sub in c("fulltext", "tables", "logs", "pubmed_abs")) {
  dir.create(file.path(out, sub), showWarnings = FALSE, recursive = TRUE)
}

# ---------------------------------------------------------- #
# 1. SEARCH
# ---------------------------------------------------------- #
query <- paste0(
  '("acute myeloid leukemia" OR "acute myeloid leukaemia" OR AML) AND ',
  '(fusion OR "gene fusion" OR "RNA fusion" OR "protein fusion" OR "RNA-seq" OR "RNA sequencing" OR ',
  '"fusion transcript" OR "chimeric RNA" OR "chimeric transcript" OR "transcriptome" OR "rearrangement")'
)
x <- entrez_search(db = "pubmed", term = query, retmax = 0)
pmids <- entrez_search(db = "pubmed", term = query, retmax = x$count)$ids
cat("Total PubMed papers:", x$count, "\n")


# ---------------------------------------------------------- #
# 2. DOWNLOAD DATA PUBMED + PMC Full Text
# ---------------------------------------------------------- #

#---------- 2(A) DOWNLOAD PubMed Abstracts ---------- 

# i) XML field extraction helpers
get1 <- function(a,path){z <- xml_find_first(a,path); if(inherits(z,"xml_missing")) return(NA_character_)
v <- xml_text(z,trim=TRUE); if(!nzchar(v)) NA_character_ else v}

getall <- function(a,path){z <- xml_find_all(a,path); if(!length(z)) return(NA_character_)
v <- xml_text(z,trim=TRUE); v <- v[nzchar(v)]; if(!length(v)) NA_character_ else paste(v,collapse=" ")}

# ii) Fetch Data: Download/save XML
batches <- split(pmids,ceiling(seq_along(pmids)/10))

fetch_pubmed_xml <- function(ids,i){
  f <- file.path(out,"pubmed_abs",paste0("batch_",i,".xml")) # path
  if(file.exists(f) && tryCatch({read_xml(f);TRUE},error=function(e)FALSE)) # files exist or not
  {cat("Skip",i,"\n");return(f)} # print progress
  # Attempt: On failing download
  for(r in 1:5){
    cat("Batch",i, "/", length(batches), "attempt",r,"\n")
    z <- tryCatch(entrez_fetch(db="pubmed",id=ids,rettype="xml",retmode="xml"), error=function(e)NULL)
    if(!is.null(z)) tryCatch({writeLines(z,f,useBytes=TRUE);read_xml(f);return(f)}, error=function(e)unlink(f))
    Sys.sleep(2*r)
  }
  stop("Failed batch ",i)
}
xml_files <- mapply(fetch_pubmed_xml,batches,seq_along(batches),SIMPLIFY=TRUE)

# iii) Parse PMID, Title, Abstract, DOI, PMC ID into a data.frame
parse_batch <- function(f){
  doc <- read_xml(f)
  arts <- xml_find_all(doc,".//PubmedArticle")
  do.call(rbind,lapply(arts,function(a)data.frame(
    PMID=get1(a,".//PMID"),
    Title=get1(a,".//ArticleTitle"),
    Abstract=getall(a,".//AbstractText"),
    DOI=get1(a,".//ArticleId[@IdType='doi']"),
    PMC_ID=get1(a,".//ArticleId[@IdType='pmc']"),
    stringsAsFactors=FALSE
  )))
}

pubmed_df <- do.call(rbind,lapply(xml_files,parse_batch))
pubmed_df <- pubmed_df[!duplicated(pubmed_df$PMID),]
cat("Parsed papers:",nrow(pubmed_df),"\n")

# Parse abstracts with pubmed.mineR
myabs <- lapply(xml_files,pubmed.mineR::new_xmlreadabs)

#---------- 2(B) DOWNLOAD PMC  Full Text ---------- 
query <- paste0(
  '("acute myeloid leukemia"[tiab] OR "acute myeloid leukaemia"[tiab] OR AML[tiab]) AND ',
  '(fusion[tiab] OR "gene fusion"[tiab] OR "RNA fusion"[tiab] OR "protein fusion"[tiab] OR ',
  '"RNA-seq"[tiab] OR "RNA sequencing"[tiab] OR "fusion transcript"[tiab] OR ',
  '"chimeric RNA"[tiab] OR "chimeric transcript"[tiab] OR transcriptome[tiab] OR "rearrangement[tiab]")'
)

x1 <- entrez_search(db="pmc",term=query,retmax=0)
pmc_ids <- entrez_search(db="pmc",term=query,retmax=x1$count)$ids
cat("PMC Title/Abstract papers:",x1$count,"\n")

# iii) Fetch PMC full-text: Download/save XML
fetch_pmc_xml <- function(pmc,i){
  f <- file.path(out,"fulltext",paste0("PMC", pmc,".xml")) # path
  if(file.exists(f) && tryCatch({read_xml(f);TRUE},error=function(e)FALSE))
  {cat("Skip",i,"\n");return(f)}
  
  for(r in 1:5){ # retry failed download
    cat("PMC",i,"/",length(pmc_ids),"attempt",r,"\n")
    z <- tryCatch(entrez_fetch(db="pmc",id=pmc,rettype="full",retmode="xml"),
                  error=function(e)NULL)
    if(!is.null(z)) tryCatch({
      writeLines(z,f,useBytes=TRUE);read_xml(f);return(f)
    },error=function(e)unlink(f))
    Sys.sleep(2*r)
  }
  stop("Failed PMC ",pmc)
}

pmc_files <- mapply(fetch_pmc_xml,pmc_ids,seq_along(pmc_ids),SIMPLIFY=TRUE)
# Parse PMC XML: PMID + PMC_ID + article body
parse_pmc <- function(f){
  doc <- read_xml(f)
  data.frame(
    PMID=get1(doc,".//article-id[@pub-id-type='pmid']"),
    PMC_ID=get1(doc,".//article-id[@pub-id-type='pmcid']"),
    FullText=tryCatch(paste(xml_text(xml_find_all(doc,".//body")),collapse=" "),
                      error=function(e) NA_character_),
    stringsAsFactors=FALSE
  )
}

pmc_df <- do.call(rbind,lapply(pmc_files,parse_pmc))
pmc_df <- pmc_df[!is.na(pmc_df$FullText) & pmc_df$FullText != "",]
cat("Parsed papers:",nrow(pmc_df),"\n")

# Common df: PubMed metadata + Abstract + matching PMC full text
df <- merge(pubmed_df,pmc_df,by="PMID",all=TRUE)
df$FullText <- ifelse(is.na(df$FullText) | !nzchar(df$FullText), NA_character_,df$FullText)
# Add SearchText for mining
df$SearchText <- apply(df[, c("Title", "Abstract", "FullText")], 1, paste, collapse = " ")


# ---------------------------------------------------------- #
# 3. FLAGS
# ---------------------------------------------------------- #
txt <- paste(ifelse(is.na(df$Abstract), "", df$Abstract),  ifelse(is.na(df$FullText), "", df$FullText))
# i) AML mention
df$AML <- grepl("acute myeloid leuk(a)?emia|acute myelogenous leuk(a)?emia|\\bAML\\b", txt, ignore.case = TRUE)
# ii) RNA-level evidence (RNA-seq, transcriptome, chimeric/fusion transcript, etc.)
df$RNA_evidence <- grepl(
  "RNA fusion|fusion RNA|fusion transcripts?|chimeric RNA|chimeric transcripts?|RNA-seq|RNA sequencing|transcriptomic?",
  txt, ignore.case = TRUE)
# iii) fusion mention
df$Fusion <- grepl("gene fusions?|fusion", txt, ignore.case = TRUE)
# iv) combined AML + RNA + fusion flag
df$RNA_fusion <- df$AML & df$RNA_evidence & df$Fusion


# ---------------------------------------------------------- #
# 4. EXTRACT GENE1--GENE2 (HGNC-normalized)
# ---------------------------------------------------------- #
# i) load HGNC old/alias -> approved symbol mapping
gene_map <- read.delim("/home/adminsb/disk2/meta-analysis/rnafusion/pubmed_miner/code_pubmedminer/HGNC_symbol_mapping.tsv", stringsAsFactors = FALSE)
gene_map[] <- lapply(gene_map, function(v) toupper(trimws(v)))

# ii) regex-scan text for GENE1-GENE2 / GENE1--GENE2 / GENE1::GENE2 / GENE1/GENE2,
#     then normalize each pair to approved HGNC symbols.
#     text_scan = regex pattern matching step; hgnc_scan = HGNC normalization step.
#     Toggle either off to run just one of the two, or leave both TRUE to run both.
extract_fusion <- function(x, text_scan = TRUE, hgnc_scan = TRUE) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  
  if (text_scan) {
    hits <- regmatches(x, gregexpr(
      "\\b[A-Za-z][A-Za-z0-9.-]{1,29}\\s*(::|--)\\s*[A-Za-z][A-Za-z0-9.-]{1,29}\\b", x, perl = TRUE
    ))[[1]]
    if (!length(hits)) return(NA_character_)
    pairs <- unique(gsub("\\s*(::|--)\\s*", "--", hits, perl = TRUE))
  } else {
    pairs <- x
  }
  
  if (hgnc_scan) {
    normalized <- vapply(strsplit(pairs, "--", fixed = TRUE), function(g) {
      if (length(g) != 2) return(NA_character_)
      approved <- gene_map$approved[match(toupper(g), gene_map$old)]
      if (anyNA(approved)) NA_character_ else paste(approved, collapse = "--")
    }, character(1))
  } else {
    normalized <- pairs
  }
  
  normalized <- unique(normalized[!is.na(normalized)])
  if (length(normalized)) paste(normalized, collapse = ";") else NA_character_
}

# iii) run extraction over every paper (both steps on by default)
df$Fusion_genes <- vapply(df$SearchText, extract_fusion, character(1),
                          text_scan = TRUE, hgnc_scan = TRUE)
# iv) expand to one row per fusion gene pair (duplicating paper metadata)
gene_list <- strsplit(df$Fusion_genes, ";", fixed = TRUE)
keep <- lengths(gene_list) > 0
df <- df[rep(which(keep), lengths(gene_list[keep])), , drop = FALSE]
df$Fusion_genes <- trimws(unlist(gene_list[keep], use.names = FALSE))
rownames(df) <- NULL


# ---------------------------------------------------------- #
# 5. OUTPUTS
# ---------------------------------------------------------- #
# i) restrict to confirmed AML+RNA+fusion papers with an extracted gene pair
rna <- df[trimws(df$Fusion_genes) != "" & !is.na(df$Fusion_genes), , drop = FALSE]
# ii) fusion gene pair frequency table
fusion_counts <- sort(table(rna$Fusion_genes), decreasing = TRUE)
write.table(fusion_counts, file.path(out, "tables", "fusion_counts_with--2semiCols.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("fusion_gene", "paper_count"))

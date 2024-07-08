#### GO Enrichment (KS Test) ####

#---- Setup ----

fun.GO.setup <- function(df,path1) {
  # Select data
  d <- dplyr::select(df,p.col)
  # Remove rows with NA value in fold change column
  d <- d[!is.na(d[["logFC"]]),]
  d[["log2FC"]] <- log2(exp(d[["logFC"]]))
  
  # map gene names to Ensembl
  gene.db <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "hgnc_symbol"),
    mart = useMart(
      "ensembl",
      dataset = "hsapiens_gene_ensembl"
      )
    )
  names(gene.db) <- c(
    "ID","GENE"
    )
  d <- left_join(
    d,
    gene.db,
    by = "GENE"
    )
  d <- d[!is.na(d[["ID"]]),]
  
  # groupwise filtering by CellType
  enrich.filt <- unique(
    d[["CellType"]]
    )
  enrich.list <- lapply(
    enrich.filt, 
    function(x) 
      dplyr::filter(
        d,
        CellType == x
        )
    )
  # Remove duplicate values
  enrich.list <- setNames(lapply(
    enrich.list,
    function (x) {
      x[!duplicated(x[["GENE"]]),]
      }
    ),c(unique(d[["CellType"]])))
  
  ## save object for running GO Enrichment
  saveRDS(
    enrich.list,
    path1
    )
  return(enrich.list)
  }


#---- Stats ----

if(Sys.info()[["sysname"]] == "Windows" & exists("GO.input")) {
  # Setup cluster (if using Windows)
  clus1 <- makeCluster(
    ifelse(
      length(
        GO.input
        ) >
        parallel::detectCores(),
      ceiling(parallel::detectCores()*0.75),
      length(
        GO.input
        )
      )
    )
  ## Cluster libraries/functions
  clusterEvalQ(
    clus1,
    {
      # libs
      library(topGO)
      library(dplyr)
      library(org.Hs.eg.db)
      library(reshape2)
      
      # GO enrichment function
      ## Definitions:
      ### df - data frame to use for enrichment
      ### id1 - column with list of gene names
      ### var1 - variables to use for topGO creation (ensembl ID, fold-change, and p-value)
      ### ont1 - GO ontology to use (BP - biological process, MF - molecular function, or CC - cellular component)
      ### db1 - database to use *In most cases this will be human "org.Hs.eg.db"
      ### map1 - mapping to be used for annotating genes (use "ensembl")
      ### desc1 - character string with description of study
      
      fun.GO.enrich <- function(df,
                            id1,
                            var1,
                            ont1,
                            db1,
                            map1,
                            desc1) {
        # use gene name as row name
        rownames(df) <- df[[id1]]
        # Select gene, log2FC, and p-value columns
        d.go1 <- setNames(df[,var1],c(id1,"FC.log2","raw.p"))
        
        
        # Create the topGO object
        ## genelist
        genelist <- setNames(d.go1[["raw.p"]],c(d.go1[[id1]]))
        
        ## topGO object
        d.go2 <- new(
          "topGOdata",
          ontology = ont1,
          allGenes = genelist,
          geneSelectionFun = function(x) x < 0.10, 
          annot = annFUN.org, 
          mapping = db1, 
          ID = map1,
          nodeSize = 5
          )
        
        # Check sig
        go.sig <- sample(
          usedGO(d.go2),
          10
          )
        go.terms <- termStat(
          d.go2,go.sig
          )
        
        # perform enrichment statistics
        
        ## Update study description
        description(d.go2) <- paste(
          description(d.go2),
          desc1
          )
        
        ## KS test (Kolmogorov-Smirnov Test)
        KS.enrich <- runTest(
          d.go2, 
          algorithm = "weight01", 
          statistic = "KS"
          )
        
        KS.tab <- GenTable(
          d.go2, 
          `P-Value` = KS.enrich,
          topNodes = length(KS.enrich@score),
          numChar = 120
          )
        
        KS.tab2 <- dplyr::select(
          KS.tab,
          GO.ID, 
          Term, 
          Annotated,
          Significant,
          `P-Value`
          )
        
        ## calculate FDR adjusted p-values
        KS.tab.out <- KS.tab2
        KS.tab.out[["FDR"]] <- p.adjust(
          KS.tab.out[["P-Value"]],
          method = "fdr"
          )
        
        ## keep pathways with raw p < 0.05 and at least 3 significant differentially expressed genes
        KS.tab.out <- dplyr::filter(
          KS.tab.out,
          `P-Value` < 0.05 & 
            Significant >= 3
          )
        KS.tab.out[["Description"]] <- desc1
        
        ## Unlist topGO object to obtain genes contained within each pathway
        genes.in.term <- genesInTerm(d.go2)
        gene.terms <- KS.tab.out[["GO.ID"]]
        genes.in.term2 <- genes.in.term[c(gene.terms)]
        genes.in.term3 <- melt(genes.in.term2)
        names(genes.in.term3) <- c(id1, "GO.ID")
        gene.sig <- sigGenes(d.go2)
        
        ## return list of GO terms for each gene
        genes.output <- merge(
          genes.in.term3,
          d.go1,
          by = id1
          )
        
        genes.output.sig <- dplyr::left_join(
          df,
          dplyr::filter(
            genes.output,
            raw.p < 0.05
          ),
          by = id1
          )
        genes.output.sig[["Description"]] <- desc1
        
        genes.output.sig <- dplyr::left_join(
          genes.output.sig,
          KS.tab.out[
            ,
            c("GO.ID","Term")
            ],
          by = "GO.ID"
          )
        
        return(
          list(
            "GO.output" = KS.tab.out,
            "GO.input" = genes.output.sig
            )
          )
        }
      }
    )
  
  ## Export global environment variables to cluster
  clusterExport(clus1,
                varlist = c("p.col"))
  
  
  # Groupwise GO enrichment by Cell Type function
  fun.GO.run.ct <- function(
    c1,
    l1,
    z
    ) {
    
    ggo <- parLapply(
      c1,
      l1,
      function (y) 
        fun.GO.enrich(
          y,
          "ID",
          c("ID","log2FC",p.col[[3]]),
          z,
          "org.Hs.eg.db",
          "ensembl",
          paste("GO",
                z,
                unique(y[[p.col[[4]]]]),
                sep = "_")
          )
      )
    return(ggo)
    }
  
  
  # Run with all ontologies
  GO.output <- list(
    ## Biological Process
    "Biological Process" = fun.GO.run.ct(
      clus1,
      GO.input,
      "BP"
      ),
    ## Molecular Function
    "Molecular Function" = fun.GO.run.ct(
      clus1,
      GO.input,
      "MF"
      ),
    ## Cell Component
    "Cell Component" = fun.GO.run.ct(
      clus1,
      GO.input,
      "CC"
      )
    )
  
  # Save enrichment results
  saveRDS(
    GO.output,
    name1
    )
  }


































## Compute all GO enrichment for each cell type using BP, MF, and CC in parallel

### Define cluster and add libraries/objects necessary for running GO enrichment

core.num <- detectCores()

clus1 <- makeCluster(ifelse(length(enrich.list) >
                              parallel::detectCores(),
                            parallel::detectCores()*0.75,
                            length(enrich.list)
)
)

clusterEvalQ(clus1,
             {
               source("1_Scripts/0_universal_lib_import.R",
                      local = knitr::knit_global())
               
               # GO enrichment function
               ## Definitions:
               ### df - data frame to use for enrichment
               ### id1 - column with list of gene names
               ### var1 - variables to use for topGO creation (ensembl ID, fold-change, and p-value)
               ### ont1 - GO ontology to use (BP - biological process, MF - molecular function, or CC - cellular component)
               ### db1 - database to use *In most cases this will be human "org.Hs.eg.db"
               ### map1 - mapping to be used for annotating genes (use "ensembl")
               ### desc1 - character string with description of study
               
               
               GO.enrich <- function(df,
                                     id1,
                                     var1,
                                     ont1,
                                     db1,
                                     map1,
                                     desc1) {
                 
                 # use gene name as row name
                 
                 rownames(df) <- df[[id1]]
                 
                 # Select gene, log2FC, and p-value columns
                 
                 d.go1 <- df[,var1]
                 
                 names(d.go1) <- c(id1,
                                   "FC.log2",
                                   "raw.p")
                 
                 
                 # Create the topGO object
                 
                 ## genelist
                 
                 genelist <- d.go1$raw.p
                 
                 names(genelist) <- d.go1[[id1]]
                 
                 ## topGO object
                 
                 d.go2 <- new("topGOdata",
                              ontology = ont1,
                              allGenes = genelist,
                              geneSelectionFun = function(x) x < 0.10, 
                              annot = annFUN.org, 
                              mapping = db1, 
                              ID = map1,
                              nodeSize = 5)
                 
                 # Check sig
                 
                 go.sig <- sample(usedGO(d.go2),
                                  10)
                 
                 go.terms <- termStat(d.go2,
                                      go.sig)
                 
                 sigGenes(d.go2)
                 
                 # perform enrichment statistics
                 
                 ## Update study description
                 
                 description(d.go2) <- paste(description(d.go2),
                                             desc1)
                 
                 description(d.go2)
                 
                 ## KS test (Kolmogorov-Smirnov Test)
                 
                 KS.enrich <- runTest(d.go2, 
                                      algorithm = "weight01", 
                                      statistic = "KS")
                 
                 KS.tab <- GenTable(d.go2, 
                                    `P-Value` = KS.enrich,
                                    topNodes = length(KS.enrich@score),
                                    numChar = 120)
                 
                 KS.tab2 <- KS.tab %>% 
                   dplyr::select(GO.ID, 
                                 Term, 
                                 Annotated,
                                 Significant,
                                 `P-Value`)
                 
                 ## calculate FDR adjusted p-values
                 
                 KS.tab.out <- KS.tab2
                 
                 KS.tab.out[["FDR"]] <- p.adjust(KS.tab.out[["P-Value"]],
                                                 method = "fdr")
                 
                 ## keep pathways with raw p < 0.05 and at least 3 significant differentially expressed genes
                 
                 KS.tab.out <- KS.tab.out %>%
                   filter(`P-Value` < 0.05 &
                            Significant >= 3)
                 
                 KS.tab.out[["Description"]] <- desc1
                 
                 ## Unlist topGO object to obtain genes contained within each pathway
                 
                 genes.in.term <- genesInTerm(d.go2)
                 
                 gene.terms <- KS.tab.out$GO.ID
                 
                 genes.in.term2 <- genes.in.term[c(gene.terms)]
                 
                 genes.in.term3 <- melt(genes.in.term2)
                 
                 names(genes.in.term3) <- c(id1, "GO.ID")
                 
                 gene.sig <- sigGenes(d.go2)
                 
                 ## return list of GO terms for each gene
                 
                 genes.output <- merge(genes.in.term3,
                                       d.go1,
                                       by = id1)
                 
                 genes.output.sig <- genes.output %>%
                   filter(raw.p < 0.05)
                 
                 genes.output.sig <- left_join(df,
                                               genes.output.sig,
                                               by = id1)
                 
                 genes.output.sig[["Description"]] <- desc1
                 
                 genes.output.sig <- left_join(genes.output.sig,
                                               KS.tab.out[,c("GO.ID",
                                                             "Term")],
                                               by = "GO.ID")
                 
                 
                 return(list(KS.tab.out,
                             genes.output.sig)
                 )
                 
               }
             })


# Groupwise GO enrichment by celltype function

GO.Group <- function(c1,
                     l1,
                     z) {
  
  ggo <- parLapply(c1,
                   l1,
                   function (y) GO.enrich(y,
                                          "ID",
                                          c("ID","log2FC",p.col[[3]]),
                                          z,
                                          "org.Hs.eg.db",
                                          "ensembl",
                                          paste("GO",
                                                z,
                                                unique(y[[p.col[[4]]]]),
                                                sep = "_")
                   )
  )
  
  return(ggo)
  
}



# Run with all ontologies

GO.comb <- list(
  # Biological Process
  
  GO.Group(clus1,
           enrich.list,
           "BP"),
  
  # Molecular Function
  
  GO.Group(clus1,
           enrich.list,
           "MF"),
  
  # Cell Component
  
  GO.Group(clus1,
           enrich.list,
           "CC")
  
)


# Save enrichment results

saveRDS(GO.comb,
        t.name)












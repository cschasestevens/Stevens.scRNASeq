#### Background DGEA Job ####
#---- Setup ----
# Load libraries
source("Scripts/4.1.processing.libs.themes.R",
       local = knitr::knit_global())

## Global parameters
list.p <- data.frame(
  # Path to integrated Seurat object
  d.path = "Processed/",
  # Number of metadata columns in input data
  md = 17,
  a.path = "Analysis/"
)

# Load data and set parameters
list.analysis <- readRDS(
  paste(
    list.p$a.path,
    "RDS/D28.original.analysis.secretory.rds",
    sep = ""
    )
  )

# Define DGEA subsets
# One variable only
l.subsets <- data.frame(
  "Group2" = c(
    "16.SecretoryCiliated|12.Secretory",
    "Secretory.set1|12.Secretory",
    "Secretory.set2|12.Secretory",
    "SecretoryCiliated|12.Secretory"
    )
  # "Knockout" = c(
  #   "KO|Control","KO|Control"
  #   )
  # "Time" = c(
  #   "D7","D7","D28","D28",
  #   "D7","D28","D7","D28",
  #   "D7|D28","D7|D28","D7|D28","D7|D28"
  # )
  )


l.subsets <- data.frame(
  "Comp.no" = seq(
    1:nrow(
      l.subsets
    )
  ),
  l.subsets,
  "Name" = paste(
    gsub(
      "\\|","vs",
      l.subsets[,1]
    ),
    # gsub(
    #   "\\|","vs",
    #   l.subsets[,2]
    # ),
    # gsub(
    #   "\\|","vs",
    #   l.subsets[,3]
    # ),
    sep = "."
  )
)


## Add MAST contrast and comparison names
l.subsets <- data.frame(
  l.subsets,
  "MAST.cont" = paste(
    as.vector(
      lapply(
        l.subsets[["Comp.no"]],
        function(x)
        {
          
          names(
            l.subsets[grepl(
              l.subsets[x,
                        grepl(
                          "\\|",
                          l.subsets[x,]
                        )],
              l.subsets
            )]
          )[1]
          
        }
      )
    ),
    ifelse(
      as.vector(
        lapply(
          l.subsets[["Comp.no"]],
          function(x)
          {
            
            names(
              l.subsets[grepl(
                l.subsets[x,
                          grepl(
                            "\\|",
                            l.subsets[x,]
                          )],
                l.subsets
              )]
            )[1]
            
          }
        )
      ) == "Group2",
      "12.Secretory",
      "Other"
    ),
    sep = ""
  ),
  "MAST.name" = l.subsets[["Name"]]
)

#---- Call DGEA functions ----
# Create DGE object from Seurat
fun.dge.subset.input <- function(
    so,
    md.list,
    c1,c2,
    md1,md2,ct.col
) {
  
  # Input Seurat
  d <- so
  
  # Subset
  d <- d[,
         rownames(
           d@meta.data
         )[grepl(
           c1,
           d@meta.data[[md1]]
         ) &
           grepl(
             c2,
             d@meta.data[[md2]]
           )]]
  
  # Data matrix
  deg.mat <- as.matrix(
    GetAssayData(
      d,
      "data",
      assay = "RNA"
    )
  )
  
  # Metadata
  deg.cols <- data.frame(
    d@meta.data[,c(
      md.list
    )]
  )
  
  # Output DGEA object (SingleCellExperiment) and Cell Type list
  deg.sc <- FromMatrix(
    deg.mat,
    cData = deg.cols
  )
  
  deg.cel <- unique(
    colData(
      deg.sc
    )[[ct.col]]
  )
  
  return(
    list(
      "SCE" = deg.sc,
      "CellGroup" = deg.cel
    )
  )
  
}


fun.dge.subset.input.group <- function(
    so,
    md.list,
    c1,c2,
    md1
) {
  
  # Input Seurat
  d <- so
  
  # Subset
  d <- d[,
         rownames(
           d@meta.data
         )[grepl(
           c1,
           d@meta.data[[md1]]
         )]]
  
  # Data matrix
  deg.mat <- as.matrix(
    GetAssayData(
      d,
      "data",
      assay = "RNA"
    )
  )
  
  # Metadata
  deg.cols <- data.frame(
    d@meta.data[,c(
      md.list
    )]
  )
  
  # Output DGEA object (SingleCellExperiment) and Cell Type list
  deg.sc <- FromMatrix(
    deg.mat,
    cData = deg.cols
  )
  
  deg.cel <- as.character(unique(
    colData(
      deg.sc
    )[[md1]]
  ))
  
  return(
    list(
      "SCE" = deg.sc,
      "CellGroup" = deg.cel
    )
  )
  
}




# Run DGEA
fun.dge.run <- function(
    ct,
    dge.object,
    form1,
    g.col,
    MAST.comps,
    MAST.comps.name
) {
  
  tryCatch(
    {
      
      # Determine core number (for parallelization)
      num.core <- detectCores()
      
      options(mc.cores = num.core*
                0.9
      )
      
      # List of comparisons for output dataframes
      l.comps <- MAST.comps
      
      # List of comparison names
      l.comp.names <- MAST.comps.name
      
      # Subset data
      s1 <- dge.object[,
                       colData(
                         dge.object
                       )[[g.col]] == ct]
      
      s1.sum <- rowSums(
        assay(
          s1
        ) >
          0
      )
      
      s1 <- s1[s1.sum/
                 ncol(
                   s1
                 ) >= 
                 0.05,]
      
      ### create glm (generalized linear model for each variable)
      
      s1.fit <- zlm(
        formula = form1,
        s1,
        method = "glm",
        ebayes = F,
        parallel = T
      )
      
      ### Output DFs
      
      d.mast.sum.fun <- function(
    comp1,
    ct2,
    comp1.name) {
        
        s1.res <- summary(
          s1.fit,
          doLRT = comp1,
          logFC = T,
          parallel = T
        )
        
        ### make dfs to display summary results by comp
        
        s1.dt <- reshape2::melt(
          dplyr::select(
            dplyr::filter(
              s1.res$datatable,
              contrast == comp1 &
                component != "S"
            ),
            -c(
              "contrast"
            )
          ),
          id.vars = c(
            "primerid","component"
          )
        )
        
        s1.dt[["vars"]] <- paste(
          s1.dt$component,
          s1.dt$variable,
          sep = "."
        )
        
        s1.dt <- dplyr::select(
          dplyr::mutate(
            reshape2::dcast(
              dplyr::select(
                dplyr::filter(
                  s1.dt,
                  vars != "logFC.Pr(>Chisq)" &
                    vars != "H.ci.hi" &
                    vars != "H.ci.lo" &
                    vars != "H.coef" &
                    vars != "H.z"
                ),
                -c(
                  "component",
                  "variable"
                )
              ),
              primerid ~ vars
            ),
            "CellType" = ct2,
            "Comparison" = comp1.name
          ),
          c(
            "CellType","Comparison","primerid",
            "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
            "D.Pr(>Chisq)",everything()
          )
        )
        
        
        names(s1.dt) <- c(
          "CellType","Comparison","GENE",
          "logFC","H.pval","C.pval",
          "D.pval",
          names(
            s1.dt[8:ncol(
              s1.dt
            )]
          )
        )
        
        s1.mis <- s1.dt[is.na(s1.dt$logFC),]
        
        list.s1 <- list(s1.dt,
                        s1.mis)
        
        return(list.s1)
        
      }
      
      dgea.comb <- dplyr::bind_rows(
        lapply(
          seq(
            1:length(
              MAST.comps
            )
          ),
          function(x) {
            
            d.mast.sum.fun(
              MAST.comps[[x]],
              ct,
              MAST.comps.name[[x]]
            )[[1]]
          }
        )
      )
      
      
      
      dgea.sum <- list(
        "DGEA.results" = dplyr::bind_rows(
          lapply(
            seq(
              1:length(
                MAST.comps
              )
            ),
            function(x) {
              
              d.mast.sum.fun(
                MAST.comps[[x]],
                ct,
                MAST.comps.name[[x]]
              )[[1]]
            }
          )
        ),
        "DGEA.missing" = dplyr::bind_rows(
          lapply(
            seq(
              1:length(
                MAST.comps
              )
            ),
            function(x) {
              
              d.mast.sum.fun(
                MAST.comps[[x]],
                ct,
                MAST.comps.name[[x]]
              )[[2]]
            }
          )
        )
      )
      
      return(dgea.sum) 
      
    },
    
    error = function(e) 
    {
      print("DGEA unsuccessful for selected cell type...")
      
    }
  )
  
}




# Run DGEA (single group)
fun.dge.run <- function(
    dge.object,
    form1,
    g.col,
    MAST.comps,
    MAST.comps.name
) {
  
  tryCatch(
    {
      
      # Determine core number (for parallelization)
      num.core <- detectCores()
      
      options(mc.cores = num.core*
                0.9
      )
      
      # List of comparisons for output dataframes
      l.comps <- l.subsets[4,2]
      
      # List of comparison names
      l.comp.names <- l.subsets[4,3]
      
      # Subset data
      s1 <- list.analysis.dgea[[4]]$SCE
      
      s1.sum <- rowSums(
        assay(
          s1
        ) >
          0
      )
      
      s1 <- s1[s1.sum/
                 ncol(
                   s1
                 ) >= 
                 0.05,]
      
      ### create glm (generalized linear model for each variable)
      
      s1.fit <- zlm(
        formula = as.formula(
          paste(
            "~","Group2","+",
            # "Airway","+",
            # "Time","+",
            "nFeature_RNA",
            sep = " "
          )
        ),
        s1,
        method = "glm",
        ebayes = F,
        parallel = T
      )
      
      ### Output DFs
      
      d.mast.sum.fun <- function(
    comp1,
    comp1.name) {
        
        s1.res <- summary(
          s1.fit,
          doLRT = comp1,
          logFC = T,
          parallel = T
        )
        
        ### make dfs to display summary results by comp
        
        s1.dt <- reshape2::melt(
          dplyr::select(
            dplyr::filter(
              s1.res$datatable,
              contrast == comp1 &
                component != "S"
            ),
            -c(
              "contrast"
            )
          ),
          id.vars = c(
            "primerid","component"
          )
        )
        
        s1.dt[["vars"]] <- paste(
          s1.dt$component,
          s1.dt$variable,
          sep = "."
        )
        
        s1.dt <- dplyr::select(
          dplyr::mutate(
            reshape2::dcast(
              dplyr::select(
                dplyr::filter(
                  s1.dt,
                  vars != "logFC.Pr(>Chisq)" &
                    vars != "H.ci.hi" &
                    vars != "H.ci.lo" &
                    vars != "H.coef" &
                    vars != "H.z"
                ),
                -c(
                  "component",
                  "variable"
                )
              ),
              primerid ~ vars
            ),
            "Comparison" = comp1.name
          ),
          c(
            "Comparison","primerid",
            "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
            "D.Pr(>Chisq)",everything()
          )
        )
        
        
        names(s1.dt) <- c(
          "Comparison","GENE",
          "logFC","H.pval","C.pval",
          "D.pval",
          names(
            s1.dt[7:ncol(
              s1.dt
            )]
          )
        )
        
        s1.mis <- s1.dt[is.na(s1.dt$logFC),]
        
        list.s1 <- list(s1.dt,
                        s1.mis)
        
        return(list.s1)
        
      }
      
      dgea.comb4 <- dplyr::bind_rows(
        lapply(
          seq(
            1:length(
              l.subsets[4,"MAST.cont"]
            )
          ),
          function(x) {
            
            d.mast.sum.fun(
              "Group2SecretoryCiliated",
              l.subsets[4,"MAST.name"][[x]]
            )[[1]]
          }
        )
      )
      
      dgea.comb4[["H.qval"]] <- p.adjust(
        dgea.comb4[["H.pval"]],
        method = "BH"
        )
      
      return(dgea.comb)
      
    },
    
    error = function(e) 
    {
      print("DGEA unsuccessful for selected cell type...")
      
    }
  )
  
}



#---- Run DGEA and save ----
## Create input SCE (subsets)
list.analysis.dgea <- setNames(
  lapply(
    l.subsets[["Comp.no"]],
    function(x) {
      
      fun.dge.subset.input.group(
        list.analysis,
        c(
          "Group2",
          # "Knockout",
          # "CellGroup",
          "nFeature_RNA"
        ),
        l.subsets[x,2],
        l.subsets[x,3],
        "Group2"
        # "Knockout",
        # "CellGroup"
      )
    }
  ),
  c(
    l.subsets[["Name"]]
  )
)

## Remove all Seurat objects prior to DGEA
remove(list.analysis)

# Run DGEA (subsets)
list.analysis.dgea.result <- lapply(
  l.subsets[["Comp.no"]],
  function(x) 
    lapply(
      list.analysis.dgea[[x]][["CellGroup"]],
      function(y) 
        fun.dge.run(
          # CellType
          y,
          # SCE input object
          list.analysis.dgea[[x]][["SCE"]],
          # GLM formula
          as.formula(
            paste(
              "~","Group2","+",
              # "Airway","+",
              # "Time","+",
              "nFeature_RNA",
              sep = " "
            )
          ),
          # Grouping column
          "Group2",
          # List of MAST comparisons (uses leading factor for names)
          c(
            l.subsets[x,"MAST.cont"]
          ),
          # List of comparison names for output
          c(
            l.subsets[x,"MAST.name"]
          )
        )
    )
)



dgea.comb5 <- dplyr::bind_rows(dgea.comb,dgea.comb2,dgea.comb3,dgea.comb4)

write.table(
  dgea.comb5,
  file = "D28.original.dgea.result.secretory.groups.txt",
  sep = "\t",
  col.names = T,
  row.names = F
  )


# Save output
saveRDS(
  dgea.comb5,
  "DGEA.result.background.rds"
  )

library(msigdbr)
library(clusterProfiler)
library(RColorBrewer)

# Define function
DGEA <- function(input_data, 
                 spec = "mouse", 
                 category = "H") {
  
  if (spec == "human") {
    oref <- HumanPrimaryCellAtlasData()
  } else if (spec == "mouse") {
    oref <- MouseRNAseqData()
  } 
  
  # Preparing clusterProfiler to perform hypergeometric test on msigdb signatures
  if (spec == "human") {
    if (category == "H") {
      m_t2g_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, human_gene_symbol)
      m_t2n_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
        dplyr::select(gs_id, gs_name) 
    }
    # msigdb signature to use
    msig_gene_set = m_t2g_h
    msig_name = m_t2n_h
    
  } 
  
  else if (spec == "mouse") {
    if (category == "H") {
      
      m_t2g_h <- msigdbr(species = "Mus musculus", category = "H") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_h <- msigdbr(species = "Mus musculus", category = "H") %>% 
        dplyr::select(gs_id, gs_name)
      #m_t2g=rbind(m_t2g_c2,m_t2g_c6)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_h
      msig_name = m_t2n_h
    }
    
    else if (category == "C2") {
      m_t2g_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
        dplyr::select(gs_id, gs_name)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_c2
      msig_name = m_t2n_c2
    }
    
    else if (category == "C7") {
      m_t2g_c7 <- msigdbr(species = "Mus musculus", category = "C7") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_c7 <- msigdbr(species = "Mus musculus", category = "C7") %>% 
        dplyr::select(gs_id, gs_name)
      #m_t2g=rbind(m_t2g_c7,m_t2g_c7)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_c7
      msig_name = m_t2n_c7
    }
  }
  
  # getting log normalized data for specific cluster
  clust_ids = sort(unique(input_data@active.ident))
  new_cluster_ids = rep(NA, length(clust_ids))
  
  # store top 30 pathway enrichment analysis
  em = NULL
  
  for (i in 1:length(clust_ids)){
    clust <- GetAssayData(subset(x = input_data,
                                 idents=clust_ids[i]),
                          slot="input_data")
    label <- rep(clust_ids[i],dim(clust)[2])
    # getting common genes
    common <- intersect(rownames(clust),
                        rownames(oref))
    # use only differential markers
    cluster_markers <- FindMarkers(input_data,
                                   ident.1 = clust_ids[i],
                                   logfc.threshold = 0.15,
                                   only.pos = TRUE)
    common <- intersect(common,
                        rownames(cluster_markers))
    oref_common <- oref[common,]
    # pred_hpca <- SingleR(test = clust, ref = hpca.se.common, labels = hpca.se$label.main,
    #                      method="cluster",clusters=label)
    # new.cluster.ids[i]=pred.hpca$labels
    tmp <- enricher(rownames(cluster_markers),
                    TERM2GENE = msig_gene_set,
                    TERM2NAME = msig_name)
    em[[i]]=tmp@result[,c("ID",
                          "p.adjust")]
  }

  # heatmap of enrichment
  library(pheatmap)
  # get top 10 enrichments
  em_table_top10 = lapply(em,function(x) x[1:10,])
  # create dataframe for heatmap
  em_hm = NULL
  em_hm = data.frame(gene_set=unique(unlist(lapply(em_table_top10,
                                                   function(x) rownames(x))
                                            )
                                     )
                     )
  
  for (i in 1:length(em_hm$gene_set)){
    for (j in 1:length(clust_ids)){
      em_hm[i,j+1] = em[[j]]$p.adjust[match(em_hm$gene_set[i],
                                            em[[j]]$ID)]
    }
  }
  rownames(em_hm)=em_hm[,1]
  em_hm=em_hm[,-1]
  em_hm[is.na(em_hm)]=1
  colnames(em_hm)=new_cluster_ids
  colnames(em_hm)=as.character(clust_ids)

  library(viridisLite)
  col_breaks=seq(-log10(1),
                 min(max(-log10(em_hm))+1,18),
                 by=0.5)
  col=inferno(length(col_breaks)) # library(viridis)
  col=c("white",
        colorRampPalette(brewer.pal(n = 7, name ="Reds"))(50))
  pheatmap(-log10(em_hm[,1:length(clust_ids)]),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           cellwidth = 10,
           cellheight = 12,
           treeheight_row = 0,
           treeheight_col=0,
           color = col,
           scale='none',
           breaks=col_breaks,
           fontsize = 8)
}


#### Two sided Barplot ####

dgea_barplot <- function(input_data, 
                         spec = "mouse", 
                         category = "H",
                         ident_1 = "4 - oHSV+Trabectedin",
                         ident_2 = "2 - oHSV") {
  output_data <- tibble()
  
  if (spec == "human") {
    oref <- HumanPrimaryCellAtlasData()
  } else if (spec == "mouse") {
    oref <- MouseRNAseqData()
  } 
  
  # Preparing clusterProfiler to perform hypergeometric test on msigdb signatures
  if (spec == "human") {
    if (category == "H") {
      m_t2g_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, human_gene_symbol)
      m_t2n_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
        dplyr::select(gs_id, gs_name) 
    }
    # msigdb signature to use
    msig_gene_set = m_t2g_h
    msig_name = m_t2n_h
    
  } 
  
  else if (spec == "mouse") {
    if (category == "H") {
      
      m_t2g_h <- msigdbr(species = "Mus musculus", category = "H") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_h <- msigdbr(species = "Mus musculus", category = "H") %>% 
        dplyr::select(gs_id, gs_name)
      #m_t2g=rbind(m_t2g_c2,m_t2g_c6)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_h
      msig_name = m_t2n_h
    }
    
    else if (category == "C2") {
      m_t2g_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
        dplyr::select(gs_id, gs_name)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_c2
      msig_name = m_t2n_c2
    }
    
    else if (category == "C7") {
      m_t2g_c7 <- msigdbr(species = "Mus musculus", category = "C7") %>%
        dplyr::select(gs_name, gene_symbol)
      m_t2n_c7 <- msigdbr(species = "Mus musculus", category = "C7") %>% 
        dplyr::select(gs_id, gs_name)
      #m_t2g=rbind(m_t2g_c7,m_t2g_c7)
      
      # msigdb signature to use
      msig_gene_set = m_t2g_c7
      msig_name = m_t2n_c7
    }
  }


  for (direction in c("up", "down")) {
    de_genes <- list()
    Idents(input_data) <- input_data$treatment
      de_genes <-
        FindMarkers(input_data,
                    ident.1 = ident_1,
                    ident.2 = ident_2,
                    only.pos = FALSE,
                    min.pct = 0.1)
      # Select genes that change in the direction of interest
      if (direction == "up") {
        de_genes <-
          de_genes[de_genes$p_val < 0.05 &
                     de_genes$avg_log2FC > 0, ]
      } else {
        de_genes <-
          de_genes[de_genes$p_val < 0.05 &
                     de_genes$avg_log2FC < 0, ]
      }
    # Preparing clusterProfiler to perform hypergeometric test on msigdb signatures
    temp <-
      enricher(row.names(de_genes),
               TERM2GENE = msig_gene_set,
               TERM2NAME = msig_name)@result %>%
      as_tibble() %>%
      #select(ID, p.adjust) %>%
      dplyr::rename(pathway = ID) %>%
      arrange(p.adjust) %>%
      slice_head(n = 5) %>%
      mutate(
        #treatment = i,
        label_y = if_else(direction == "up", 1, -1),
        p.adjust = -log10(p.adjust),
        pathway = str_remove(pathway, "HALLMARK_") %>%
          str_replace_all("_", " "),
        order = seq_len(n())
        )
    if (direction == "down") {
      temp$p.adjust <- temp$p.adjust * -1
    }
    output_data <- bind_rows(output_data, temp)
  }
  
  # Two sided Barplot
  # Properly name and order treatment types
  output_data <- output_data %>%
    mutate(
      # treatment = stringr::str_replace(treatment,
      #                                  "2 - oHSV",
      #                                  "oHSV Gene Sets") %>%
      #   stringr::str_replace("4 - oHSV+Trabectedin",
      #                        "oHSV+Trabectedin Gene Sets"),
      pathway = stringr::str_wrap(pathway, width = 25))
  
  # output_data$treatment <- factor(as.factor(output_data$tissue),
  #                                 levels = c("oHSV Gene Sets",
  #                                            "oHSV+Trabectedin Gene Sets"))
  
  # Create tibble that is properly ordered to use in proper plotting of data
  lab4plot <- tibble(y = c(-4, 4),
                     x = c(0.5, 0.5),
                     label = as.factor(
                       c("Downregulated\nin oHSV+Trabectedin compared to oHSV alone",
                         "Upregulated\nin oHSV+Trabectedin compared to oHSV alone")
                     )
                     # tissue = factor(as.factor(c("Tibia Gene Sets",
                     #                             "Tibia Gene Sets",
                     #                             "Lung Gene Sets",
                     #                             "Lung Gene Sets")),
                     #                 levels = c("Tibia Gene Sets",
                     #                            "Lung Gene Sets")
                     #                 )
  )
  
  figure <-
    ggplot() +
    geom_bar(data = output_data,
             aes(x = -1 * order,
                 y = p.adjust,
                 fill = p.adjust > 0),
             stat = "identity",
             alpha = 0.8) +
    coord_flip() +
    # facet_wrap(~ tissue,
    #            ncol = 1,
    #            scales = "free") +
    geom_hline(yintercept = c(-1.5, 1.5),
               color = "gray",
               linetype = 2,
               linewidth = 0.5) +
    geom_hline(yintercept = 0,
               color = "black",
               linewidth = 0.5) +
    geom_text(data = output_data,
              aes(x = -1 * order,
                  y = label_y * 4,
                  label = pathway),
              size = 5) +
    geom_text(data = lab4plot,
              aes(x = x,
                  y = y,
                  label = label),
              fontface = "bold",
              size = 5) +
    scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
    theme_bw() +
    theme(strip.background = element_rect(color = "white", fill = "white"),
          strip.text.x = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(size = 4),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = "",
         y = "",
         x = "") +
    ylim(-7, 7) +
    xlim(-5.5, 1)
  print(figure)
  # Take a moment to ensure reset to null device
  #dev.off()
}
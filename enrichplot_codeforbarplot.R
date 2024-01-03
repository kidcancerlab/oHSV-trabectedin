##' barplot of enrichResult
##'
##'
##' @importFrom graphics barplot
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_col
##  @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 scale_y_discrete
##' @title barplot
##' @param height enrichResult object
##' @param x one of 'Count' and 'GeneRatio'
##' @param color one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of categories to show
##' @param font.size font size
##' @param title plot title
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' by default wraps names longer that 30 characters
##' @param ... other parameter, ignored
##' @method barplot enrichResult
##' @export
##' @return ggplot object
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' barplot(x)
##' # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
##' barplot(x, showCategory = 10)
##' categorys <- c("pre-malignant neoplasm", "intestinal disease",
##'                "breast ductal carcinoma", "non-small cell lung carcinoma")
##' barplot(x, showCategory = categorys)
barplot.enrichResult <- function(height, x="NES", color="NES",
                                 showCategory=10, font.size=12, title="",
                                 label_format=30, ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue", "NES"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  else if (x == "nes" || x == "NES") {
    x <- "NES"
  }
  dframe <- fortify(object, showCategory=showCategory, by=x, ...)
  
  if(colorBy %in% colnames(dframe)) {
    p <- ggplot(dframe, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(low="#357EBDFF", high="#D43F3AFF", name = color,
                            guide=guide_colorbar(reverse=FALSE),
                            limits = c(-3, 3))
  } else {
    p <- ggplot(dframe, aes_string(x = x, y = "Description",
                               fill = "Description")) +
      theme_dose(font.size) +
      theme(legend.position="none")
  }
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  p + geom_col() + # geom_bar(stat = "identity") + coord_flip() +
    scale_y_discrete(labels = label_func) +
    ggtitle(title) + xlab(NULL) + ylab(NULL)
}


barplot.compareClusterResult <- function(height, color="NES",
                                         showCategory=5, by="NES",
                                         includeAll=TRUE, font.size=12,
                                         title="", ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an compareClusterResult object.
  dframe <- fortify(height, showCategory=showCategory, by=by,
                includeAll=includeAll)
  plotting.clusterProfile(dframe, type="bar", colorBy=color, by=by, title=title,
                          font.size=font.size)
}

#### background needed ####

prepare_pie_gene <- function(y) {
  gene_pie <- tibble::as_tibble(y[,c("Cluster", "Description", "geneID")])
  gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
  gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols=geneID))
  gene_pie2 <- unique(gene_pie2)
  prepare_pie_data(gene_pie2, pie =  "equal", type = "gene")
}


##' Prepare pie data for categories in cnetplot/emapplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default)
##' or 'Count'
##' @return a data.frame
##' @noRd
prepare_pie_category <- function(y, pie = "equal") {
  pie <- match.arg(pie, c("equal", "count", "Count"))
  if (pie == "count") pie <- "Count"
  
  pie_data <- y[,c("Cluster", "Description", "Count")]
  pie_data[,"Description"] <- as.character(pie_data[,"Description"])
  prepare_pie_data(pie_data, pie = pie)
}




prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
  if(type == "category"){
    ID_unique <- unique(pie_data[,2])
  } else {
    ID_unique <- unique(pie_data[,3])
  }
  
  Cluster_unique <- unique(pie_data[,1])
  ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
  rownames(ID_Cluster_mat) <- ID_unique
  colnames(ID_Cluster_mat) <- Cluster_unique
  ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
  if(pie == "Count") {
    for(i in seq_len(nrow(pie_data))) {
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
    }
    for(kk in seq_len(ncol(ID_Cluster_mat))) {
      ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
    }
    return(ID_Cluster_mat)
  }
  for(i in seq_len(nrow(pie_data))) {
    if(type == "category"){
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
    } else {
      ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
    }
    
  }
  return(ID_Cluster_mat)
}


##' create color palette for continuous data
##'
##'
##' @title color_palette
##' @param colors colors of length >=2
##' @return color vector
##' @export
##' @examples
##' color_palette(c("red", "yellow", "green"))
##' @author guangchuang yu
color_palette <- function(colors) {
  # has_package("grDevices")
  grDevices::colorRampPalette(colors)(n = 299)
}


sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)
  
  if(x@readable) {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}

# fc_palette <- function(fc) {
# if (all(fc > 0, na.rm=TRUE)) {
# palette <- color_palette(c("blue", "red"))
# } else if (all(fc < 0, na.rm=TRUE)) {
# palette <- color_palette(c("green", "blue"))
# } else {
## palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
# }
# return(palette)
# }

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}

extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title plotting-clusterProfile
##' @param clProf.reshape.df data frame of compareCluster result
##' @param x x variable
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 guide_colorbar
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotting.clusterProfile <- function(clProf.reshape.df,
                                    x = ~Cluster,
                                    type = "dot",
                                    colorBy = "NES",
                                    by = "rowPercentage",
                                    title="",
                                    font.size=12) {
  Description <- Percentage <- Count <- Cluster <- GeneRatio <- NES <- p.adjust <- pvalue <- NULL # to satisfy codetools
  if (type == "bar") {
    if (by == "percentage") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Percentage, fill=Cluster))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Count, fill=Cluster))
    } else if (by == "NES") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = NES, fill=Cluster))
    } else if (by == "geneRatio") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = geneRatio, fill=Cluster))
    } else {
      print("define what to plot using `by=`")
    }
    p <- p +
      geom_bar() +
      coord_flip()
  }
  if (type == "dot") {
    if (by == "rowPercentage") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Percentage))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Count))
    } else if (by == "geneRatio") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~GeneRatio))
    } else if (by == "NES") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~NES))
    } else {
      ## nothing here
    }
    if (any(colnames(clProf.reshape.df) == colorBy)) {
      p <- p +
        geom_point() +
        aes_string(color=colorBy) +
        scale_color_continuous(low="blue", high="red",
                               guide=guide_colorbar(reverse=FALSE))
      ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = sig_palette)
    } else {
      p <- p + geom_point(colour="steelblue")
    }
  }
  p <- p + xlab(paste(by)) + ylab("") + ggtitle(title) +
    theme_dose(font.size)
  ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
  ##     theme(axis.text.y = element_text(colour="black",
  ##           size=font.size, hjust = 1)) +
  ##               ggtitle(title)+theme_bw()
  ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
  ##                    hjust=hjust.axis.x,
  ##                    vjust=vjust.axis.x))
  return(p)
}




##' Get the distance of the label
##'
##' @param dimension one of 1 and 2
##' @param label_location label_location
##' @noRd
get_label_diss <- function(dimension, label_location) {
  nn <- nrow(label_location)
  label_dis <- matrix(NA, nrow = nn, ncol = nn)
  colnames(label_dis) <- rownames(label_dis) <- label_location$label
  for (i in seq_len(nn - 1)) {
    for (j in (i + 1):nn) {
      label_dis[i ,j] <- label_location[i, dimension] -  label_location[j, dimension]
    }
  }
  label_diss <- reshape2::melt(label_dis)
  label_diss <- label_diss[label_diss[,1] != label_diss[,2], ]
  label_diss <- label_diss[!is.na(label_diss[,3]), ]
  label_diss[, 1] <- as.character(label_diss[, 1])
  label_diss[, 2] <- as.character(label_diss[, 2])
  return(label_diss)
}



# adjust_location <- function(label_location, x_adjust, y_adjust) {
# label_diss_x <- get_label_diss(1, label_location)
# label_diss_y <- get_label_diss(2, label_location)

# label_diss_large <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
# intersect(which(label_diss_y[, 3] > 0)) %>%
# intersect(which(abs(label_diss_x[, 3]) < x_adjust))

# label_diss_small <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
# intersect(which(label_diss_y[, 3] < 0)) %>%
# intersect(which(abs(label_diss_x[, 3]) < x_adjust))

# label_location[label_diss_y[label_diss_large, 1], 2] <- label_location[label_diss_y[label_diss_large, 2], 2] + y_adjust
# label_location[label_diss_y[label_diss_small, 1], 2] <- label_location[label_diss_y[label_diss_small, 2], 2] - y_adjust
# return(label_location)
# }


#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @noRd
ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `ep_str_wrap`
#' @noRd
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}


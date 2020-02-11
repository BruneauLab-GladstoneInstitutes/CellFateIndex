library(Seurat)

#' Create start and target centroids
#' @param sobj
#' @param start
#' @param s_degs
#' @param target
#' @param t_degs
#' @return centroid_list
#' @export 
create.centroids <- function(sobj, start, s_degs, target, t_degs){
    centroid_list <- list(
        'start_centroid_s_degs' = AverageExpression(object=subset(sobj, idents=start), features= row.names(s_degs), 
                                                     assays = 'RNA'),
        'start_centroid_t_degs' = AverageExpression(object=subset(sobj, idents=start), features= row.names(t_degs), 
                                                     assays = 'RNA'),
        'target_centroid_t_degs' = AverageExpression(object=subset(sobj, idents=target), features= row.names(t_degs), 
                                                     assays = 'RNA'),
        'target_centroid_s_degs' = AverageExpression(object=subset(sobj, idents=target), features= row.names(s_degs), 
                                                     assays = 'RNA'))
    return(centroid_list)
}

#' Create weighted gene vector
#' @param centroid_list
#' @return weighted_gene_vector_list
#' @export
create.weighted.gene.vector <- function(centroid_list){
    weighted_gene_vector_list <- lapply(seq_along(centroid_list),
                                     function(y,i){apply(y[[i]]$RNA, 2,
                                                         function(t) t/colSums(y[[i]]$RNA))}, 
                                                         y=centroid_list)
    return(weighted_gene_vector_list)
}

#' Calculate cell fate index j
#' @param r_ij
#' @param gene_weight
#' @return r_ij
#' @export                                                   
calc.cfi.j <- function(r_ij, gene_weight){
    if (r_ij < 0){
        r_ij = 0
    } else if (r_ij > 0 & r_ij < 1){
        r_ij = r_ij * gene_weight
    } else {
        r_ij = 1 * gene_weight
    }
    return (r_ij)
}

#' Run cell fate index j
#' @param exp_matrix
#' @param centroid_1
#' @param centroid_2
#' @param weighted_gene_vector
#' @return cfi_j
#' @export                  
run.cfi.j <- function(exp_matrix, centroid_1, centroid_2, weighted_gene_vector){
    cfi_j <- lapply(seq_along(row.names(centroid_1$'RNA')),
                   function(y,i){
                       cfi_j <- setNames(as.data.frame((exp_matrix[y[i],] - centroid_1$'RNA'[y[i],])/
                                                     (centroid_2$'RNA'[y[i],] - centroid_1$'RNA'[y[i],])), c('cfi'))
                       cfi_j$cfi <- sapply(cfi_j$cfi, function(x) calc.cfi.j(x, weighted_gene_vector[y[i],]))
                   },
                   y=row.names(centroid_1$'RNA'))
    # Bind list of dataframes by columns and add the gene names (columns) and cell barcodes (rows)
    cfi_j <- bind_cols(cfi_j)
    names(cfi_j) <- row.names(centroid_1$'RNA')
    row.names(cfi_j) <- colnames(exp_matrix)
    return(cfi_j)
}

#' Run cell fate index i
#' @param cfi_j
#' @return cfi_is
#' @export
run.cfi.i <- function(cfi_j){
    cfi_i <- setNames(rowSums(cfi_j) %>% as.data.frame(), c('cfi'))
    return(cfi_i)
}

#' Calculate CFI
#' @param sobj
#' @param start
#' @param s_degs
#' @param target
#' @param t_degs
#' @param cfi_col_name
#' @export
cfi <- function(sobj, start, s_degs, target, t_degs, cfi_col_name){
    centroid_list <- create.centroids(sobj, start,  s_degs, target, t_degs)
    weighted_gene_vector_list <- create.weighted.gene.vector(centroid_list)
    cfi_j <- run.cfi.j(GetAssayData(sobj), centroid_list$start_centroid_t_degs, centroid_list$target_centroid_t_degs, 
                       weighted_gene_vector_list[[3]])
    cfi_i <- run.cfi.i(cfi_j)
    sobj <- AddMetaData(sobj, cfi_i, cfi_col_name)
    return(sobj)
}
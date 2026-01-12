#' Plot overall heatmaps and export rankings and top lists (custom nshow)
#'
#' @description
#' Given three metric lists (`all_result_acc`, `all_result_recall`, `all_result_FS`),
#' this function builds matrices, moves the training column to the end, renames
#' columns to "(train set)/(val set)/(test set)", ranks models by the mean over
#' test-set columns, draws heatmaps, and writes:
#' (1) full ranking CSV; (2) filtered ranking CSV; (3) Top-N models CSV
#' (emitted inside the first heatmap). Two PDFs per metric are written:
#' the full heatmap and the Top-N heatmap. All files are saved under `out_dir`.
#'
#' @param all_result_acc list
#'   Accuracy results per method. Each element is a named numeric vector whose
#'   names are dataset columns and values are scores.
#' @param all_result_recall list
#'   Recall results per method, same structure as `all_result_acc`.
#' @param all_result_FS list
#'   F-score results per method, same structure as `all_result_acc`.
#' @param exp_files character
#'   Expression filename used to derive display names with
#'   `gsub(".txt", "", exp_files)`.
#' @param heatmapcolor character
#'   Color vector for the column annotation. Length should be >= number of cols.
#' @param nshow integer
#'   Number of models in the Top-N heatmap and the Top-N CSV. Default is 40.
#' @param namesS character
#'   Names of the three metrics, in order. Default:
#'   `c("Accuracy","Recall","F-score")`.
#' @param out_dir character
#'   Output directory for all files. Will be created if it does not exist.
#'   Default is `"heatmap"`.
#'
#' @details
#' Files written under `out_dir`:
#' * `all_<Metric>_ordered_models.csv`        full ranking (by test mean)
#' * `filtered_<Metric>_ordered_models.csv`   filtered ranking (same order)
#' * `<nshow>_top_<Metric>_models.csv`        Top-N model names
#' * `all_<Metric>_statical.pdf`              full heatmap PDF
#' * `<nshow>.show_<Metric>_statical.pdf`     Top-N heatmap PDF
#'
#' Column renaming convention:
#' the first column is tagged with " (train set)",
#' the second with " (val set)",
#' columns 3..k with " (test set)".
#'
#' @return
#' Invisibly returns a list of length 3, the matrices used to plot each
#' metric heatmap, in the same order as `namesS`.
#'
#' @seealso
#' \code{\link{fh_roc}}, \code{\link{fh_plot_model}}
#'
#' @family FeatureHunter plotting
#'
#' @import ComplexHeatmap
#' @importFrom grid unit gpar grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' heatmap_cols <- c("#8ecae6","#219ebc","#ffc8dd","#ffafcc",
#'                   "#023047","#ffb703","#fb8500","#bde0fe")
#'
#' fh_all_heatmap(
#'   all_result_acc    = acc_list,
#'   all_result_recall = rec_list,
#'   all_result_FS     = fs_list,
#'   exp_files         = "DatasetA.txt",
#'   heatmapcolor      = heatmap_cols,
#'   nshow             = 40,
#'   namesS            = c("Accuracy","Recall","F-score"),
#'   out_dir           = "heatmap"
#' )
#' }
#'
#' @export
fh_all_heatmap <- function(
    all_result_acc,
    all_result_recall,
    all_result_FS,
    exp_files,
    heatmapcolor,
    nshow = 40,
    namesS = c("Accuracy", "Recall", "F-score"),
    out_dir           = "heatmap"
){
  if (!dir.exists("heatmap")) dir.create("heatmap")
  # ---- 与原脚本一致：派生 name_exp ----
  name_exp <- gsub(".txt", "", exp_files)
  
  # ---- 与原脚本一致：先把三个列表打包 ----
  all_result_list <- list(all_result_acc, all_result_recall, all_result_FS)
  
  # ---- 完全照你的 lapply 块构造 all_result_tt ----
  all_result_tt <- lapply(
    all_result_list,
    function(listdata, methodname, name_exp){
      all_k <- t(as.data.frame(listdata))
      rownames(all_k) <- methodname
      colnames(all_k) <- name_exp
      return(all_k)
    },
    methodname = names(all_result_acc),
    name_exp   = name_exp
  )
  
  # ---- 三大评分绘图循环（保持你的原逻辑）----
  for (si in 1:3) {
    statical_mat <- all_result_tt[[si]]
    colnames(statical_mat)[1] <- paste0(colnames(statical_mat)[1], " (train set)")
    colnames(statical_mat)[2] <- paste0(colnames(statical_mat)[2], " (val set)")
    if (ncol(statical_mat) >= 3) {
      colnames(statical_mat)[3:ncol(statical_mat)] <- paste0(
        colnames(statical_mat)[3:ncol(statical_mat)], " (test set)"
      )
    }
    
    ordered_models <- rownames(statical_mat)[order(
      apply(statical_mat[, 3:ncol(statical_mat), drop = FALSE], 1, mean),
      decreasing = TRUE
    )]
    utils::write.csv(
      data.frame(Model = ordered_models),
      file = file.path("heatmap",paste0("all_", namesS[si], "_ordered_models.csv")),
      row.names = FALSE
    )
    
    DTCol <- heatmapcolor[1:ncol(statical_mat)]
    names(DTCol) <- colnames(statical_mat)
    
    filtered_models <- rownames(statical_mat)[order(
      apply(statical_mat[, 3:ncol(statical_mat), drop = FALSE], 1, mean),
      decreasing = TRUE
    )]
    utils::write.csv(
      data.frame(Model = filtered_models),
      file = file.path("heatmap",paste0("filtered_", namesS[si], "_ordered_models.csv")),
      row.names = FALSE
    )
    
    avg_statical <- apply(statical_mat[, 3:ncol(statical_mat), drop = FALSE], 1, mean)
    avg_statical <- sort(avg_statical, decreasing = TRUE)
    statical_mat <- statical_mat[names(avg_statical), , drop = FALSE]
    avg_statical <- as.numeric(format(avg_statical, digits = 3, nsmall = 3))
    
    # 与原脚本一致：将 nshow 与行数取最小
    nshow_use <- min(nshow, nrow(statical_mat))
    
    hm <- ComplexHeatmap::Heatmap(
      as.matrix(statical_mat),
      name = namesS[si],
      col = c("#20ACBD", "#FFFFFF", "#F69896"),
      row_gap = grid::unit(1, "mm"), column_gap = grid::unit(3, "mm"),
      rect_gp = grid::gpar(col = "grey", lwd = 1),
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = grid::gpar(col = "black", border = "grey"),
      row_names_max_width = grid::unit(15, "cm"),
      width  = grid::unit(ncol(statical_mat) + 3, "cm"),
      height = grid::unit(nrow(statical_mat)/1.5, "cm"),
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 15),
        legend_height = grid::unit(15, "cm"),
        grid_width = grid::unit(0.8, "cm")
      ),
      column_split = factor(colnames(statical_mat), levels = colnames(statical_mat)),
      row_split    = factor(rownames(statical_mat), levels = rownames(statical_mat)),
      cluster_columns = FALSE,
      cluster_rows    = FALSE,
      row_title = NULL,
      column_title = NULL,
      right_annotation = ComplexHeatmap::rowAnnotation(
        'Average of test set' = ComplexHeatmap::anno_barplot(
          avg_statical, bar_width = 0.8, border = TRUE,
          gp = grid::gpar(fill = "#B84D64", col = "grey", border = "grey"),
          add_numbers = TRUE, numbers_offset = grid::unit(-10, "mm"),
          axis_param = list("labels_rot" = 0),
          numbers_gp = grid::gpar(fontsize = 10, col = "white"),
          width = grid::unit(4, "cm"), height = grid::unit(1.1, "cm")
        ),
        show_annotation_name = TRUE
      ),
      top_annotation = ComplexHeatmap::columnAnnotation(
        "Data set" = colnames(statical_mat),
        col = list("Data set" = DTCol),
        show_annotation_name = FALSE,
        annotation_legend_param = list(
          title_gp = grid::gpar(fontsize = 15),
          grid_height = grid::unit(1, "cm"),
          grid_width  = grid::unit(1, "cm")
        )
      ),
      layer_fun = function(j, i, x, y, w, h, fill) {
        labs <- format(statical_mat[cbind(i, j)], digits = 3, nsmall = 3)
        grid::grid.text(label = labs, x = x, y = y, gp = grid::gpar(fontsize = 10))
        if (exists("nshow")) {
          top_models <- rownames(statical_mat)[1:nshow_use]
          utils::write.csv(
            data.frame(Model = top_models),
            file = file.path("heatmap", paste0(nshow_use, "_top_", namesS[si], "_models.csv")),
            row.names = FALSE
          )
        }
      }
    )
    
    grDevices::pdf(
      file = file.path("heatmap",paste0("all_", namesS[si], "_statical.pdf")),
      width = ncol(statical_mat) * 1.2 + 12,
      height = nrow(statical_mat) / 3.5
    )
    ComplexHeatmap::draw(hm)
    invisible(grDevices::dev.off())
    
    # ---- 第二轮：Top-N 热图（与你脚本的第二块一致）----
    statical_mat <- all_result_tt[[si]]
    colnames(statical_mat)[1] <- paste0(colnames(statical_mat)[1], " (train set)")
    colnames(statical_mat)[2] <- paste0(colnames(statical_mat)[2], " (val set)")
    if (ncol(statical_mat) >= 3) {
      colnames(statical_mat)[3:ncol(statical_mat)] <- paste0(
        colnames(statical_mat)[3:ncol(statical_mat)], " (test set)"
      )
    }
    DTCol <- heatmapcolor[1:ncol(statical_mat)]
    names(DTCol) <- colnames(statical_mat)
    
    avg_statical <- apply(statical_mat[, 3:ncol(statical_mat), drop = FALSE], 1, mean)
    avg_statical <- sort(avg_statical, decreasing = TRUE)
    statical_mat <- statical_mat[names(avg_statical), , drop = FALSE]
    avg_statical <- as.numeric(format(avg_statical, digits = 3, nsmall = 3))
    
    statical_mat <- statical_mat[1:nshow_use, , drop = FALSE]
    avg_statical <- avg_statical[1:nshow_use]
    
    nshow_use2 <- min(nshow_use, nrow(statical_mat))
    
    hm <- ComplexHeatmap::Heatmap(
      as.matrix(statical_mat),
      name = namesS[si],
      col = c("#20ACBD", "#FFFFFF", "#F69896"),
      row_gap = grid::unit(1, "mm"), column_gap = grid::unit(3, "mm"),
      rect_gp = grid::gpar(col = "grey", lwd = 1),
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = grid::gpar(col = "black", border = "grey"),
      row_names_max_width = grid::unit(15, "cm"),
      width  = grid::unit(ncol(statical_mat) * 2.5, "cm"),
      height = grid::unit(nrow(statical_mat)/1.5, "cm"),
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 15),
        legend_height = grid::unit(15, "cm"),
        grid_width = grid::unit(0.8, "cm")
      ),
      column_split = factor(colnames(statical_mat), levels = colnames(statical_mat)),
      row_split    = factor(rownames(statical_mat), levels = rownames(statical_mat)),
      cluster_columns = FALSE,
      cluster_rows    = FALSE,
      row_title = NULL,
      column_title = NULL,
      right_annotation = ComplexHeatmap::rowAnnotation(
        'Average of test set' = ComplexHeatmap::anno_barplot(
          avg_statical, bar_width = 0.8, border = TRUE,
          gp = grid::gpar(fill = "#B84D64", col = "grey", border = "grey"),
          add_numbers = TRUE, numbers_offset = grid::unit(-10, "mm"),
          axis_param = list("labels_rot" = 0),
          numbers_gp = grid::gpar(fontsize = 10, col = "white"),
          width = grid::unit(4, "cm"), height = grid::unit(1.1, "cm")
        ),
        show_annotation_name = TRUE
      ),
      top_annotation = ComplexHeatmap::columnAnnotation(
        "Data set" = colnames(statical_mat),
        col = list("Data set" = DTCol),
        show_annotation_name = FALSE,
        annotation_legend_param = list(
          title_gp = grid::gpar(fontsize = 15),
          grid_height = grid::unit(1, "cm"),
          grid_width  = grid::unit(1, "cm")
        )
      ),
      layer_fun = function(j, i, x, y, w, h, fill) {
        labs <- format(statical_mat[cbind(i, j)], digits = 3, nsmall = 3)
        grid::grid.text(label = labs, x = x, y = y, gp = grid::gpar(fontsize = 10))
      }
    )
    
    grDevices::pdf(
      file = file.path("heatmap",paste0(nshow_use2, ".show_", namesS[si], "_statical.pdf")),
      width = ncol(statical_mat) * 1.2 + 12,
      height = nrow(statical_mat) / 3.5
    )
    ComplexHeatmap::draw(hm)
    invisible(grDevices::dev.off())
  }
  
  invisible(all_result_tt)
}
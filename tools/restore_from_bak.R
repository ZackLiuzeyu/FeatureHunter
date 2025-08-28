# restore_from_bak.R
# 说明：
# - 在指定目录(默认 R/) 里找到所有 *.bak
# - 把它们还原成去掉 .bak 的原名（覆盖当前文件）
# - 可 dry-run 预览；也可在覆盖前把当前文件再备份成 .now

restore_from_bak <- function(
    dir = "R",
    recursive = TRUE,
    dry_run = TRUE,
    backup_current_as_now = TRUE
) {
  bak_files <- list.files(dir, pattern = "\\.bak$", recursive = recursive, full.names = TRUE)
  if (!length(bak_files)) {
    message("No .bak files found under: ", normalizePath(dir))
    return(invisible(character()))
  }
  
  msg <- function(...) cat(sprintf(...), "\n")
  
  to_restore <- data.frame(
    bak      = bak_files,
    original = sub("\\.bak$", "", bak_files),
    stringsAsFactors = FALSE
  )
  
  if (dry_run) {
    msg("[Dry-run] %d files would be restored:", nrow(to_restore))
    print(to_restore)
    msg("Nothing changed. Set dry_run = FALSE to actually restore.")
    return(invisible(to_restore))
  }
  
  restored <- character(0)
  for (i in seq_len(nrow(to_restore))) {
    bak <- to_restore$bak[i]
    ori <- to_restore$original[i]
    
    # 如目标已存在，可先备份一份 .now
    if (file.exists(ori) && backup_current_as_now) {
      now <- paste0(ori, ".now")
      ok_now <- tryCatch(file.copy(ori, now, overwrite = TRUE), error = function(e) FALSE)
      if (!ok_now) warning("Could not backup current file to: ", now)
    }
    
    ok <- tryCatch(file.copy(bak, ori, overwrite = TRUE), error = function(e) FALSE)
    if (ok) {
      restored <- c(restored, ori)
      message("[restored] ", ori)
    } else {
      warning("Failed to restore: ", bak, " -> ", ori)
    }
  }
  
  message("Done. Restored: ", length(restored), " files.")
  invisible(restored)
}

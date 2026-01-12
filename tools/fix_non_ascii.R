# 在包根目录执行（含 DESCRIPTION、R/、NAMESPACE）
escape_non_ascii <- function(path) {
  x <- readBin(path, "raw", file.info(path)$size)  # raw 读，避免编码干扰
  s <- rawToChar(x, multiple = FALSE)
  Encoding(s) <- "UTF-8"  # 明确标注
  
  # 把所有非 ASCII 字符替换成 \uXXXX
  esc <- function(txt) {
    ints <- utf8ToInt(txt)
    paste0(vapply(
      ints,
      function(i) if (i > 127) sprintf("\\u%04x", i) else intToUtf8(i),
      character(1)
    ), collapse = "")
  }
  
  s2 <- esc(s)
  if (!identical(s, s2)) {
    writeLines(s2, path, useBytes = TRUE)
    message("[fixed] ", path)
  }
}

# 遍历 R/ 下所有 .R 与 NAMESPACE
files <- c(
  list.files("R", pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE),
  "NAMESPACE"
)
files <- files[file.exists(files)]

lapply(files, escape_non_ascii)

# 可视化哪些文件仍含非 ASCII（应为空）
to_check <- files[file.exists(files)]
non_ascii_left <- Filter(function(f) length(tools::showNonASCIIfile(f)) > 0, to_check)
if (length(non_ascii_left)) {
  message("Still has non-ASCII: ", paste(non_ascii_left, collapse = ", "))
} else {
  message("All code/NAMESPACE are now ASCII-only (via \\uXXXX).")
}



fix_rd_unicode <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # 常见的几个
  x <- gsub("\\\\u201c|\\\\u201d", "\"", x)   # 左/右弯双引号 → ASCII "
  x <- gsub("\\\\u2018|\\\\u2019", "'", x)   # 左/右弯单引号 → ASCII '
  x <- gsub("\\\\u2013|\\\\u2014", "--", x)  # 短/长破折号 → --
  x <- gsub("\\\\u2026", "...", x)           # 省略号 → ...
  writeLines(x, path, useBytes = TRUE)
  message("[fixed Rd] ", path)
}

rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
invisible(lapply(rd_files, fix_rd_unicode))



#（撤销）
# 在包根目录执行（包含 DESCRIPTION、R/、NAMESPACE、man/ 等）
unescape_unicode_in_files <- function(
    paths = c("R", "NAMESPACE", "man"),         # 需要处理的路径/文件
    backup = TRUE,                               # 是否备份为 .bak
    verbose = TRUE
) {
  # 列出待处理文件
  expand_targets <- function(p) {
    if (file.exists(p) && file.info(p)$isdir) {
      list.files(p, pattern = "\\.(r|R|Rd|rd|Rmd|qmd|Rnw|tex)$", full.names = TRUE, recursive = TRUE)
    } else if (file.exists(p)) {
      p
    } else {
      character(0)
    }
  }
  files <- unique(unlist(lapply(paths, expand_targets)))
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    message("No files to process. Check your 'paths' argument.")
    return(invisible(character(0)))
  }
  
  # 将 \uXXXX / \UXXXXXXXX 还原为实际字符
  unescape_unicode <- function(txt) {
    # \UXXXXXXXX（8位）要先处理，避免与 \uXXXX 冲突
    txt <- gsubfn::gsubfn(
      "\\\\U([0-9a-fA-F]{8})",
      function(h) intToUtf8(as.integer(strtoi(h, 16L))),
      txt
    )
    # 再处理 \uXXXX（4位）
    txt <- gsubfn::gsubfn(
      "\\\\u([0-9a-fA-F]{4})",
      function(h) intToUtf8(as.integer(strtoi(h, 16L))),
      txt
    )
    txt
  }
  
  changed <- character(0)
  for (f in files) {
    # 以原始方式读，避免编码干扰
    raw <- readBin(f, what = "raw", n = file.info(f)$size)
    s   <- rawToChar(raw, multiple = FALSE)
    Encoding(s) <- "UTF-8"  # 统一视为 UTF-8 文本
    
    s_new <- unescape_unicode(s)
    
    if (!identical(s, s_new)) {
      if (backup && !file.exists(paste0(f, ".bak"))) {
        file.copy(f, paste0(f, ".bak"), overwrite = FALSE)
      }
      # 覆写为 UTF-8
      con <- file(f, open = "wb", encoding = "UTF-8")
      writeLines(s_new, con, useBytes = TRUE)
      close(con)
      if (verbose) message("[unescaped] ", f)
      changed <- c(changed, f)
    }
  }
  
  if (verbose) {
    if (length(changed)) {
      message("Done. Changed files: ", length(changed))
    } else {
      message("Nothing changed (no \\uXXXX / \\UXXXXXXXX sequences found).")
    }
  }
  invisible(changed)
}

# 依赖 gsubfn（只在脚本运行时用，不影响包 DESCRIPTION）
if (!requireNamespace("gsubfn", quietly = TRUE)) {
  install.packages("gsubfn")
}

# 执行：R/、NAMESPACE、man/ 都还原
unescape_unicode_in_files(paths = c("R", "NAMESPACE", "man"))

### 0) 重启 R 之后执行本段（很重要）

### 1) 列出所有库路径，彻底删除 FeatureHunter 及锁目录
suppressWarnings(try(remove.packages("FeatureHunter"), silent = TRUE))

all_libs <- .libPaths()
for (lib in all_libs) {
  unlink(file.path(lib, "FeatureHunter"), recursive = TRUE, force = TRUE)
  unlink(file.path(lib, "00LOCK-FeatureHunter"), recursive = TRUE, force = TRUE)
}

# 有些系统会留下一堆 00LOCK*
locks <- unique(unlist(lapply(all_libs, list.files, pattern="^00LOCK", full.names=TRUE)))
if (length(locks)) unlink(locks, recursive = TRUE, force = TRUE)

# 额外：把所有库里可能残留的 .rdb/.rdx 单独扫一遍（保险起见）
rdbs <- unique(unlist(lapply(all_libs, function(lib) {
  dir(file.path(lib, "FeatureHunter"), pattern="\\.rd(b|x)$", full.names=TRUE)
})))
if (length(rdbs)) unlink(rdbs, force = TRUE)

cat("清理完成。库中是否还存在包目录：\n")
for (lib in all_libs) if (dir.exists(file.path(lib, "FeatureHunter"))) cat(" - ", lib, "\n")

### 2) 把“用户库”放第一位（以后都往这里装）
userlib <- "~/Library/R/arm64/4.4/library"
dir.create(path.expand(userlib), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(path.expand(userlib), setdiff(.libPaths(), path.expand(userlib))))
print(.libPaths())

### 3) 干净重装（不往系统库写）
if (!"remotes" %in% rownames(installed.packages())) install.packages("remotes")
Sys.setenv(R_INSTALL_STAGED = "false")           # 避免 staged install 偶发抽风
remotes::install_local(".", force = TRUE, upgrade = "never", build_vignettes = FALSE)

### 4) 验证加载的确来自“用户库”
library(FeatureHunter)
cat("包实际位置：", system.file(package="FeatureHunter"), "\n")


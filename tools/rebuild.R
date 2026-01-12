## 0) 看看现在 R 会去哪些库装包
.libPaths()

## 1) 确保没有加载
if ("FeatureHunter" %in% loadedNamespaces()) {
  try(detach("package:FeatureHunter", unload=TRUE, character.only=TRUE), silent=TRUE)
  try(unloadNamespace("FeatureHunter"), silent=TRUE)
}

## 2) 删除所有库路径里的 FeatureHunter 目录 + 遗留的 00LOCK*
for (lib in .libPaths()) {
  cat("Cleaning lib:", lib, "\n")
  try(unlink(file.path(lib, "FeatureHunter"), recursive = TRUE, force = TRUE), silent = TRUE)
  try(unlink(Sys.glob(file.path(lib, "00LOCK*")), recursive = TRUE, force = TRUE), silent = TRUE)
}

## 3)（可选）你的系统库路径是：
## /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
## 如果上面 unlink 没权限，确认你当前用户对该目录有写权限，
## 或者把 .libPaths() 指到你有权限的用户库:
## .libPaths("~/Library/R/arm64/4.4/library")

## 4) 清理开发目录里的旧检查/构建残留（在你的包源代码目录里执行）
## 如果你是在包源码目录里运行：
try(unlink("FeatureHunter.Rcheck", recursive = TRUE, force = TRUE), silent = TRUE)
try(unlink("00LOCK-FeatureHunter", recursive = TRUE, force = TRUE), silent = TRUE)

## 5) 强制一次“干净安装”
## 建议先生成文档，再安装；避免 vignette 拖慢或失败

rm(list = ls()[sapply(ls(), function(x) is.function(get(x)))])

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::document(quiet = TRUE)
  ## 不建 vignette；不升级依赖；避免多架构标志
  devtools::install(upgrade = "never",
                    dependencies = NA,
                    build_vignettes = FALSE,
                    args = c("--no-multiarch"))
} else {
  install.packages("devtools")
}







































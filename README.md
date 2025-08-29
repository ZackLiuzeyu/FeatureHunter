# FeatureHunter


**FeatureHunter** 是一个“机器学习 + 深度学习”一体化的特征发现与可解释性分析工具箱。  
它提供 **28 个模型接口**（**1 个 MLP + 27 个经典ML**），并通过核心函数 **`fh_hunter()`** 从排行榜里获取最优模型，统一计算三类特征重要性（模型内置、Permutation、SHAP），用稳健 z 分数（median/MAD）做 **复合打分**，再结合多次重复训练（不同种子）和自举置信区间完成 **稳健特征筛选**、可视化与下游验证（UMAP、Logistic 回归公式等）。

## 安装

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)
devtools::install_github("ZackLiuzeyu/FeatureHunter")
library(FeatureHunter)
```

**数据要求**，
## 数据准备 / Input Data
1. 输入至少 **4 个数据集**：  
   - 1 个作为训练集，命名为 **DatasetA**  
   - 1 个作为验证集，命名为 **DatasetB**  
   - 2 个作为测试集，命名为 **DatasetC** 和 **DatasetD**  
   - 文件分别保存为：`DatasetA.txt`、`DatasetB.txt`、`DatasetC.txt`、`DatasetD.txt`  
2. 输入至少 **4 个样本信息表**，格式类似于示例文档：  
   - 命名为：`SampleA.txt`、`SampleB.txt`、`SampleC.txt`、`SampleD.txt`
3. 输入一个 **公共基因列表文件**：  
   - 命名为：`Common_genes.txt`
  
最后把以上所有文档单独放到一个干净的文件夹内，并将文件夹命名为raw_data,在使用时在相应位置填写该文件夹路径，即可开始分析
  
     
1. Input at least **4 datasets**:  
      - 1 as training set called **DatasetA**  
      - 1 as validation set called **DatasetB**  
      - 2 as test sets called **DatasetC** and **DatasetD**  
      - Save as: `DatasetA.txt`, `DatasetB.txt`, `DatasetC.txt`, `DatasetD.txt`  
2. Input at least **4 sample information tables**, similar to the example document:  
      - Named: `SampleA.txt`, `SampleB.txt`, `SampleC.txt`, `SampleD.txt`  
3. Input a **common gene list file**:  
      - Named: `Common_genes.txt`
        
Finally, put all the above documents into a clean folder, and name the folder raw_data, fill in the folder path in the corresponding location when using, you can start the analysis

---

**补充说明 / Notes**：  
- 如果提供超过 4 个数据集，从第 **3 个数据集开始**及其之后的所有数据集都将作为测试集。  
- *If more than 4 datasets are provided, the 3rd and subsequent datasets will be treated as test sets.*  

## 模型接口一览（28）
- **深度学习**：`fh_mlp()`（SMLP 权重 / 输入梯度 / fastshap SHAP）  
- **经典 ML（节选）**：rf / xgboost / svm / glmnet / gbm / C5.0 / nnet / earth / glmStepAIC / sda / kknn / bayesglm / catboost 等，均带有 `*_best` 或 `lasso_*` 变体以支持 LASSO 预筛选与 CV 调参。  
> 完整函数清单见本仓库 `R/` 目录与下文“API 参考”。

## 工作流总览（`fh_hunter()`）
1. **自动选型**：从排行榜（CSV）中按指标（默认 F1）选择第 `pick_index` 行对应的最优模型；支持任意模型类型（MLP / RF / GLMNET / XGBoost / SVM / LDA / QDA / NaiveBayes…）。  
2. **重复训练与集成**：使用 `num_runs` 个随机种子重复训练；先对预测做 **集成（ensemble-first）**，再在集成预测器上计算 Permutation / SHAP，降低方差。  
3. **三路重要性**：
   - **模型内置**（系数/权重/gain 等）  
   - **Permutation**（打乱单特征造成的性能下降）  
   - **SHAP**（模型无关解释，`fastshap`）  
4. **复合打分与置信区间**：对每次 run 的三路重要性做 **稳健 z 分数**；按 `method_weights` 加权融合为 **Composite**；对每个基因的复合分数跨 run 进行 **bootstrap**，得到均值与 `ci_level` 置信区间。  
5. **筛选控制（可选）**：用一组 “**严格模式**” 开关控制最终 Top 集合：
   - `strict_min_freq`：跨 run 命中频率阈值  
   - `strict_min_effect`：复合均值下限  
   - `strict_ci_gate`：要求 CI 下界 > 0  
   - `strict_knee` + `strict_nmax`：在排序曲线自动找折点（knee）并设最大数量上限  
   - `apply_selection_to_plots`：将筛选后的 Top 集用于柱图/UMAP（默认仅展示原始 Top）  
6. **可视化 & 下游**：

### 图表产出总览

运行 `fh_mlp()` / `fh_run_ml_models()` / `fh_hunter()` 过程中，会自动生成下列核心图表（保存至 `out_dir`）：

- **LASSO 预筛选回归路径图**：展示逐步引入特征的系数路径变化  
- **模型性能热图**：所有模型的 Accuracy、F-score、Recall 三项指标排序热图  
- **ROC 曲线图**：各模型的 ROC 曲线与 AUC 值对比  
- **PR 曲线图**：各模型的 Precision–Recall 曲线与 AUC-PR 值对比  
- **混淆矩阵**：逐模型的分类混淆矩阵，附带 ACC / F1 / REC 指标  
- **训练曲线**：MLP 的训练 Accuracy / Loss 随 epoch 变化曲线  
- **特征重要性柱图**：按复合分数排序的 Top 特征及置信区间柱状图  
- **三方法密度分布图**：模型内置 / Permutation / SHAP 三种重要性的分布对比  
- **UMAP 降维图**：基于 Top 基因特征的二维可视化，区分样本类别  
- **稳定性热图**：跨 run 的 Top 基因出现频率热图  

> 这些图表既能用于方法学报告（补充材料/附录），也能直接嵌入论文结果部分。


### 前置筛选（两条路）
- **路线 A · MLP 专项筛选（`fh_mlp()`）**  
- **重要性三件套**：SMLP 权重、输入梯度、fastshap SHAP；先对每 run 做稳健 z-score，再融合为复合分数。
- **阈值**：支持 F1 / Youden 自动或手动；适配类别不平衡。
- **复现性**：多 seed 重复、可选 bootstrap 置信区间；固定 `seed` 可复现。

- **路线 B · 27 个经典 ML 批量筛选（`fh_run_ml_models()`）**  
- **覆盖面**：一次性并行/批量跑 **27 个经典 ML**（含 `*_best` 调参款与 `lasso_*` 预筛选款）。
- **调参策略**：内置交叉验证或栅格搜索，选择 `best` 变体的最优超参。
- **统一评测**：以 Accuracy / Recall / F1（可含 ROC-AUC / PR-AUC）打分，生成 **leaderboard CSV**。
- **可控随机性**：固定随机种子，多次重复以降低方差。
- **落盘**：结果与关键图表/CSV 自动保存，供 `fh_hunter()` 后续读取复用。

> 两条路都会产出一个**排行榜（CSV）**与若干可视化。`fh_hunter()` 会 **自动读取该排行榜**，根据你设定的指标与行号（`score_index`/`pick_index`）锁定候选最优模型，然后进入它的统一重要性计算与下游流程。

   - Top20 箱线图、Composite 柱图（均值+95%CI）、三法分布对比密度、UMAP（Top 基因）、稳定性热图（跨 run）  
   - 选定 Top 基因后进行 **Logistic 回归**，输出 **系数表与可直接引用的公式**（便于撰写论文方法/公式部分）。

## 生成的文件（默认输出到 `out_dir`）
- `FI_Boxplot_Top20.pdf` — Top-20 特征重要性分布箱线图  
- `FI_Bar_TopComposite.pdf` — 复合重要性 Top-N 柱图（均值 + 95%CI）  
- `FI_Density_AllMethods.pdf` — 三种重要性方法的分布对比  
- `UMAP_TopSignatureGenes.pdf` — 依据 Top 基因的 UMAP 投影  
- `Stability_TopGenes_Heatmap.pdf` — 跨 run 的 Top 基因稳定性热图  

## `fh_hunter()` 参数（节选）
- `train_exp` — Numeric matrix of predictors (samples x features).
- `train_labels` — Response vector (0/1, factor, or convertible).
- `nshow` — Number of models shown in leaderboard.
- `namesS` — Character vector of metric names shown in leaderboard
- `score_index` — Integer index of metric to rank by (default: 3 = F-score).
- `pick_index` — Row index of leaderboard to pick (default: 1). Works for any model type.
- `top_models_csv` — Path to leaderboard CSV (auto-inferred if \code{NULL}).
- `num_runs` — Number of repeat runs (seeds) (default: 10).
- `num_coregene` — Number of core genes to keep in stability/importance plots.
- `n_likes` — Number of top genes to use in downstream analysis (UMAP, logistic).
- `n_interest` — Number of genes to show in the stability heatmap.
- `seed` — Random seed (default: 424).
- `out_dir` — Output directory for plots (default: working dir).
- `perm_metric` — Performance metric for permutation importance (\code{"f1"} or \code{"prauc"}).
- `perm_nrep` — Number of repetitions for permutation (default: 3).
- `shap_nsim` — Number of Monte Carlo simulations for SHAP (default: 10).
- `shap_subsample` — Number of samples for SHAP subsampling (default: 100).
- `svm_cost` — Cost parameter for SVM (default: 1).
- `svm_gamma` — Gamma for SVM-RBF (default: \code{1/p} if \code{NULL}).
- `shuffle` — Logical; whether to shuffle data each epoch in MLP training (default: \code{TRUE}).
- `standardize` — logical; if TRUE, standardize using train-set mean/sd and apply to val/test.
- `hidden_units` — integer vector; number of units per hidden layer,
- `activation` — character; activation function for hidden layers,
- `use_batchnorm` — logical; whether to insert BatchNorm after each hidden Dense.
- `l2` — numeric; L2 regularization strength (e.g. 1e-4). Default: NULL (auto).
- `gaussian_noise_sd` — numeric; stddev for input GaussianNoise layer.
- `min_lr` — numeric; minimum learning rate for ReduceLROnPlateau.
- `plateau_factor` — numeric; factor for ReduceLROnPlateau (e.g. 0.5).
- `plateau_patience` — integer; patience (epochs) for ReduceLROnPlateau.
- `early_patience` — integer; patience (epochs) for EarlyStopping.
- `imbalance_thresh` — numeric in \eqn{[0,1]}; threshold for enabling imbalance handling.
- `auto_th_method` — character; "youden", "f1", or "auto" for thresholding.
- `method_weights` — length-3 named numeric vector for composite weighting,
- `ci_level` — numeric; CI level for composite bootstrap over runs (default 0.95).
- `strict_min_freq` — numeric; minimum cross-run hit frequency in per-run Top-K to keep a gene
- `strict_min_effect` — numeric; minimum composite mean to keep a gene on the composite scale.
- `strict_ci_gate` — logical; if TRUE, require composite CI lower bound > 0.
- `strict_knee` — logical; if TRUE, apply a knee cutoff on sorted composite means.
- `strict_nmax` — integer; optional hard cap after knee (e.g. 10). NULL disables.
- `apply_selection_to_plots` — logical; if TRUE, bar/UMAP use selected genes;

> **Composite 权重**：`method_weights = c(model = 1, perm = 1, shap = 1)`；可按任务调节三者权重。  
> **随机性**：通过 `num_runs` + `seed` 控制；推荐固定 `seed` 并报告 `num_runs`。  

## 返回值
- `params`：解析后的参数与元数据  
- `importance_df`：每个基因的三路重要性与复合分数  
- `top_list`：每次 run 的 Top 基因列表  
- `final_top`：应用严格筛选后的最终 Top 集合  
- `composite_mat`：逐 run 的复合分数矩阵  
- `glm_summary`：基于最终 Top 基因的 Logistic 回归结果（含公式）  

## 快速上手示例

```r
# X: 样本 x 基因 矩阵；Y: 二分类标签（0/1 或 factor）
res <- fh_hunter(
  train_exp = X, train_labels = Y,
  nshow = 40, score_index = 3, pick_index = 1,
  num_runs = 5, num_coregene = 50, n_likes = 30, n_interest = 30,
  strict_min_freq = 0.65, strict_ci_gate = TRUE, strict_knee = TRUE, strict_nmax = 10
)

head(res$importance_df)
res$glm_summary$formula   # 论文可直接引用的公式
```

> 如果你只想直接用 MLP：请看 `fh_mlp()`，支持 SMLP/Gradient/SHAP 三路重要性与固定风格绘图；超参既可从排行榜自动解析，也可手动覆盖。

## 实用提示
- **阈值自动化**：`auto_th_method` 根据类别不平衡与可用指标在 F1 / Youden 之间自动选择，亦可手动指定。  
- **PR/ROC**：内部已切换为专业包（如 `precrec` / `pROC`）计算，确保与业界指标对齐。  
- **并行**：若配置了 `doParallel` / `foreach`，部分耗时阶段可并行；请在调用前注册集群。  
- **可复现性**：固定 `seed`、记录 `method_weights` 与严格筛选门限，建议在论文附录公开 `importance_df` 与公式。  

---

最后更新：2025-08-28 16:24  
欢迎提 Issue / PR：**ZackLiuzeyu/FeatureHunter**

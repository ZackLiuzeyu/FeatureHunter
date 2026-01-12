str(besttune)
n_trees <- besttune$n.trees
c(mode = mode(n_trees), class = class(n_trees), value = as.character(n_trees))
as.integer(n_trees)  # 看看它到底是多少


mean(as_binary01(train_labels))

setdiff(model$var.names, colnames(val_x))        # 训练用到但验证缺的列
setdiff(colnames(val_x), model$var.names)        # 验证多出来的列（一般无所谓）

summary(sapply(val_x, function(col) c(var = stats::var(col, na.rm=TRUE), na = sum(is.na(col)))))


length(comd_val); head(comd_val); length(unique(comd_val))
stopifnot(all(rownames(val_x) %in% comd_val))






# 截距对应的概率（gbm 的 initF 是 logit 标度）
plogis(model$initF)

# 训练集上的预测是否也“全常数”？
prob_tr_check <- as.numeric(gbm::predict.gbm(model,
                                             newdata = as.data.frame(train_exp)[, lassogene, drop = FALSE],
                                             type = "response", n.trees = as.integer(besttune$n.trees)
))
length(unique(round(prob_tr_check, 6)))
head(unique(round(prob_tr_check, 6)), 10)

str(besttune)                      # 你贴过是 numeric=100，这里再确认一次
n_trees <- as.integer(besttune$n.trees)

# 逐步试一下 1、10、n_trees 棵树，看看唯一值个数
u1  <- length(unique(round(as.numeric(gbm::predict.gbm(model,
                                                       newdata = val_x, type = "response", n.trees = 1L)), 6)))
u10 <- length(unique(round(as.numeric(gbm::predict.gbm(model,
                                                       newdata = val_x, type = "response", n.trees = 10L)), 6)))
uB  <- length(unique(round(as.numeric(gbm::predict.gbm(model,
                                                       newdata = val_x, type = "response", n.trees = n_trees)), 6)))
c(u1=u1, u10=u10, uBest=uB)


table(train_expd_final$labels, useNA = "ifany")
str(train_expd_final$labels)

imp <- summary(model, plotit = FALSE)
head(imp, 10)




# 看验证集每个 LASSO gene 的方差
sort(sapply(val_x, var))

# 换一个测试集 C / D 来预测，看看是不是还全常数
prob_c <- as.numeric(predict(model,
                             newdata = as.data.frame(test_exp[[2]])[, lassogene, drop=FALSE],
                             type="response", n.trees=n_trees))
length(unique(round(prob_c, 6)))

qq_train <- t(sapply(as.data.frame(train_expd_final)[, lassogene, drop=FALSE],
                     function(v) quantile(v, c(.01,.5,.99), na.rm=TRUE)))
qq_val   <- t(sapply(val_x, function(v) quantile(v, c(.01,.5,.99), na.rm=TRUE)))

head(round(cbind(
  train_p01 = qq_train[,1], train_p50 = qq_train[,2], train_p99 = qq_train[,3],
  val_p01   = qq_val[,1],   val_p50   = qq_val[,2],   val_p99   = qq_val[,3]
), 6))

rng_train <- sapply(as.data.frame(train_expd_final)[, lassogene, drop=FALSE], range, na.rm=TRUE)
rng_val   <- sapply(val_x, range, na.rm=TRUE)
round(rbind(train_min = rng_train[1,], val_min = rng_val[1,],
            train_max = rng_train[2,], val_max = rng_val[2,]), 6)[, 1:6]  # 只看前6个基因

















































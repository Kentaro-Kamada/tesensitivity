library(tidyverse)
library(haven)
library(quantreg)
library(WeightIt)
library(kamaken)


data <- 
  read_dta('tesensitivityStataPackage/lalonde1986.dta') |> 
  filter(sample1 == 1) |> 
  mutate(id = row_number())


nodes <- 20

# 傾向スコアの推定
logit <- glm(treat ~ married + age + black + hispanic + education + re74 + re75 + re74pos + re75pos, data = data, family = binomial)

p1 <- predict(logit, type = "response")
p0 <- 1 - p1

# 分位点回帰
# 分位点回帰を行い、その後にそれらを積分してATEを求める
# 数値積分する際の区分点＝補完ノード：今回はチェビシェフノードを用いる（Wikipedia参照）

# 補完ノードの設定
k <- 0:(nodes-1)
xnodes <- 1/2 + 1/2 * cos((2*k + 1)*pi / (2 * nodes))

# 各ノードに対して分位点回帰を行う
res <- 
  map(xnodes, \(tau) rq(
    formula = re78 ~ treat + married + age + black + hispanic + education + re74 + re75 + re74pos + re75pos,
    data = data,
    tau = tau
  )) |> 
  map(broom::augment) |> 
  bind_rows() |> 
  arrange(.tau) |> 
  distinct(.tau, .fitted)
  


pracma::pchipfun(res$.tau, res$.fitted) |> 
  integrate(0.05, 0.95)

# c-valueとも組み合わせる
  cross_join(tibble(c = seq(0, 1, 0.1))) 


# c-valueとtauの組み合わせを作成
cross_join(
  tibble(c = seq(0, 1, 0.1)),
  tibble(tau = xnodes) 
) |> 
  cross_join(
    tibble(
      id = 1:length(p1),
      p1 = p1,
      p0 = p0
    ) 
  ) |> 
  # t(tau,p1,p0)の値を計算
  mutate(
    t1_upper = pmap_dbl(list(c, tau, p1) , \(c, tau, p1) min(tau + c*min(tau, 1-tau)/p1, tau/p1, 1)),
    t1_lower = pmap_dbl(list(c, tau, p1) , \(c, tau, p1) max(tau - c*min(tau, 1-tau)/p1, (tau - 1)/p1 + 1, 0)),
    t0_upper = pmap_dbl(list(c, tau, p0) , \(c, tau, p0) min(tau + c*min(tau, 1-tau)/p0, tau/p0, 1)),
    t0_lower = pmap_dbl(list(c, tau, p0) , \(c, tau, p0) max(tau - c*min(tau, 1-tau)/p0, (tau - 1)/p0 + 1, 0))
  ) 
  mutate(
    
  )



# パッケージの読み込み
library(pracma)

# データ点
x <- c(0, 1, 1, 2, 3, 4, 5)
y <- c(0, 1, 1, 4, 2, 3, 5)

# 補間点
x_interp <- seq(0, 5, length.out = 100)

# PCHIP補間
y_interp <- pracma::pchip(x, y, x_interp)


tibble(x_interp, y_interp) |> 
  ggplot(aes(x = x_interp, y = y_interp)) +
  geom_line() +
  geom_point(data = tibble(x, y), aes(x = x, y = y), color = "red")



integrate(pracma::pchipfun(x, y), 0, 5)


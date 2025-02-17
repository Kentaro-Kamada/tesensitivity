


# 必要なパッケージの読み込み
library(MASS)
library(quantreg)
library(stats)

# 関数の定義
tesen_cpi <- function(data, treatment, outcome, covariates, qcovariates = NULL, stat = 'ate',
                      median = FALSE, quantile = 0.5, cgrid = 40, creference = FALSE, 
                      breakdown = 0, nobreakdown = FALSE, nodes = 100, tol = 0.001, 
                      verbose = FALSE, debug = FALSE) {
  
  # =========================================================================
  # 1. Define temporary variables
  # =========================================================================
  
  # Temporary variables in R
  xw0 <- xw1 <- pindex <- p0 <- p1 <- pmax <- p11 <- p10 <- xnodes <- ynodes <- coef <- NULL
  tau_n_0 <- tau_n_1 <- exp_y0 <- exp_y1 <- up0 <- up1 <- c_table <- c_search_table <- NULL
  lb <- ub <- b <- V <- wsupp <- call <- cref <- bscoef0 <- bscoef1 <- NULL
  
  # =========================================================================
  # 2. Parse Input
  # =========================================================================
  
  if (debug) {
    cat("command line input is:", deparse(substitute(data)), "\n")
  }
  
  # Parse model specification
  Y <- outcome
  X <- treatment
  W <- covariates
  
  # Check if covariates are the same for both models
  if (!all(W == qcovariates)) {
    tmvarlist <- union(W, qcovariates)
    cat("Currently, covariates in the outcome and treatment models must be the same, will use union of both varlists.\n")
    cat("Covariates used are:", tmvarlist, "\n")
  }
  
  if (is.null(qcovariates)) {
    qmodel <- "interaction"
  } else {
    qmodel <- ""
  }
  
  # calculate dimensions
  nobs <- nrow(data)
  K <- length(W)
  nvars <- K + 2
  
  # Process verbose + debug option
  verbose <- as.logical(verbose)
  debug <- as.logical(debug)
  
  # Process the specification of q(x,w) for the linear quantile regression
  if (qmodel == "interaction") {
    qmvarlist <- c(W, paste0("c.", X, "#c.(", W, ")"))
    nvars <- 2 * length(W) + 2
    if (debug) {
      cat("Covariates are", W, "\n")
      cat("Full varlist is", qmvarlist, "\n")
    }
  } else {
    qmvarlist <- W
    nvars <- length(W) + 2
    if (debug) {
      cat("Full varlist is", qmvarlist, "\n")
    }
  }
  
  # Process nobreakdown option
  if (nobreakdown) {
    breakdown <- NA
  }
  
  # Check if treatment variable is binary
  if (!all(data[[X]] %in% c(0, 1))) {
    stop("Treatment variable must be binary (0 or 1).")
  }
  
  # Check if outcome variable is binary
  binary_outcome <- all(data[[Y]] %in% c(0, 1))
  if (binary_outcome && !all(data[[Y]] %in% c(0, 1))) {
    stop("Outcome variable must be binary (0 or 1).")
  }
  
  # Check for multiple statistics and store statistic type
  if (length(stat) > 1) {
    stop("Only one statistic can be calculated in a single call to tesensitivity cpi.")
  }
  cond_stat <- any(c("cate", "cqte") %in% stat)
  qtl_stat <- any(c("qte", "cqte") %in% stat)
  
  # Store the name of the c_function to use
  c_function <- paste0(stat, "_bounds_c")
  if (binary_outcome) {
    c_function <- paste0("binary_", c_function)
  }
  
  # Check if binary outcome, check only ate specified
  if (binary_outcome && stat != "ate") {
    stop("For binary outcome, only ate is currently supported.")
  }
  
  # Form covariate matrix for conditional statistics
  if (cond_stat) {
    if (median || covariates == "median" || qcovariates == "median") {
      wsupp <- apply(data[W], 2, median)
    } else {
      wsupp <- colMeans(data[W])
    }
    
    # Overwrite covariate matrix with custom direct values if specified
    if (length(covariates) == 1 && !covariates %in% c("mean", "median")) {
      wsupp <- matrix(covariates, nrow = 1)
      colnames(wsupp) <- W
    } else if (length(covariates) > 1) {
      wsupp <- matrix(covariates, nrow = 1)
      colnames(wsupp) <- W
    }
    
    # Overwrite covariate matrix with custom quantile values if specified
    if (length(qcovariates) == 1 && qcovariates != "median") {
      for (i in seq_along(W)) {
        pct <- qcovariates[i] * 100
        wsupp[i] <- quantile(data[[W[i]]], probs = pct / 100)
      }
    } else if (length(qcovariates) > 1) {
      for (i in seq_along(W)) {
        pct <- qcovariates[i] * 100
        wsupp[i] <- quantile(data[[W[i]]], probs = pct / 100)
      }
    }
    
    # Add a constant to covariate support
    wsupp <- cbind(wsupp, 1)
    colnames(wsupp) <- c(W, "_const")
    
    if (debug) {
      cat("Covariate support to be used is:\n")
      print(wsupp)
    }
  } else {
    wsupp <- as.matrix(data[W])
    wsupp <- cbind(wsupp, 1)
    colnames(wsupp) <- c(W, "_const")
  }
  
  # =========================================================================
  # 3. Calculate propensity score(s)
  # =========================================================================
  logit_model <- glm(as.formula(paste(X, "~", paste(W, collapse = "+"))), data = data, family = binomial)
  # pindex <- wsupp %*% coef(logit_model)
  # p1 <- plogis(pindex)
  p1 <- predict(logit_model, type = "response")
  p0 <- 1 - p1
  pmax <- max(c(p1, p0))
  
  if (cond_stat && debug) {
    cat("Propensity scores are p0:", p0, "p1:", p1, "\n")
  }

  # =========================================================================
  # 4. Steps required for particular statistics
  # =========================================================================
  
  # =========================================================================
  # 4.1 (Averaged stats: calculate interpolation nodes)
  # =========================================================================
  if (!qtl_stat && !binary_outcome) {
    if (verbose) cat("Calculating interpolation nodes...\n")
    
    xnodes <- ynodes <- NULL
    for (j in seq_len(nodes)) {
      x_j <- -cos((2 * j - 1) * pi / (2 * nodes))
      x_j <- (x_j + 1) / 2
      qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = x_j)
      xnodes <- c(xnodes, x_j)
      ynodes <- rbind(ynodes, coef(qreg_model))
      if (verbose) cat(".")
    }
  }
  
  # pchip係数を計算します
  pchipcoefs <- pchip_all_coef(xnodes, ynodes, nvars)

  # =========================================================================
  # 4.2 (Quantile stats: approximate quantile function limits)
  # =========================================================================
  if (stat == "cqte") {
    alpha_accel <- 0.5
    iterate <- 30
    
    # Approximate quantile coefficients for tau -> 1
    tau <- 0.01
    qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = 1 - tau)
    bscoef1 <- coef(qreg_model)
    if (verbose) cat("Calculating approximate upper limit quantile function...\n")
    for (k in seq_len(iterate)) {
      if (verbose) cat(".")
      tau <- alpha_accel * tau
      qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = 1 - tau)
      bscoef1 <- coef(qreg_model)
    }
    
    # Approximate quantile coefficients for tau -> 0
    tau <- 0.01
    qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = tau)
    bscoef0 <- coef(qreg_model)
    if (verbose) cat
    
    # Approximate quantile coefficients for tau -> 0
    tau <- 0.01
    qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = tau)
    bscoef0 <- coef(qreg_model)
    if (verbose) cat("Calculating approximate lower limit quantile function...\n")
    for (k in seq_len(iterate)) {
      if (verbose) cat(".")
      tau <- alpha_accel * tau
      qreg_model <- rq(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, tau = tau)
      bscoef0 <- coef(qreg_model)
    }
  }
  
  # =========================================================================
  # 4.3 (ATET: unconditional expectations)
  # =========================================================================
  if (stat == "atet") {
    reg_model <- lm(as.formula(paste(Y, "~", X)), data = data)
    exp_y1 <- predict(reg_model, newdata = data.frame(X = 1))
    exp_y0 <- predict(reg_model, newdata = data.frame(X = 0))
    
    up1 <- mean(data[[X]])
    up0 <- 1 - up1
  }
  
  # =========================================================================
  # 4.4 (Binary ATE: outcome probabilities)
  # =========================================================================
  if (stat == "ate" && binary_outcome) {
    logit_model <- glm(as.formula(paste(Y, "~", paste(c(X, qmvarlist), collapse = "+"))), data = data, family = binomial)
    
    nobs <- nrow(wsupp)
    nwsupp <- ncol(wsupp)
    
    if (qmodel == "interaction") {
      xw0 <- cbind(0, wsupp[, 1:(nwsupp - 1)], 0 * wsupp[, 1:(nwsupp - 1)], 1)
      xw1 <- cbind(1, wsupp[, 1:(nwsupp - 1)], 1 * wsupp[, 1:(nwsupp - 1)], 1)
    } else {
      xw0 <- cbind(0, wsupp)
      xw1 <- cbind(1, wsupp)
    }
    
    p11 <- plogis(xw1 %*% coef(logit_model))
    p10 <- plogis(xw0 %*% coef(logit_model))
  }
  
  # =========================================================================
  # 5. Form matrix of c-dependence values
  # =========================================================================
  
  cgrid_values <- seq(0, 1, length.out = cgrid)
  call <- matrix(cgrid_values, ncol = 1)
  
  if (creference) {
    cref <- max(c(p1, p0))
    call <- rbind(call, cref)
  }
  
  call <- rbind(call, pmax)
  cnum <- nrow(call)
  
  # =========================================================================
  # 6. Calculate c-dependence bounds on statistics
  # =========================================================================

  c_function_args <- list()
  if (stat == "cqte") {
    c_function_args <- list(quantile, p0, p1, wsupp, bscoef0, bscoef1, Y, X, qmvarlist, data)
  } else if (stat == "atet") {
    c_function_args <- list(p0, p1, wsupp, xnodes, ynodes, nvars, nobs, exp_y0, exp_y1, up0, up1)
  } else if (stat == "qte") {
    c_function_args <- list(quantile, p0, p1, wsupp, tol, Y, X, qmvarlist, data)
  } else if (stat == "ate" && binary_outcome) {
    c_function_args <- list(p0, p1, p10, p11)
  } else {
    c_function_args <- list(p0, p1, wsupp, pchipcoefs, nvars, nobs)
  }
  
  c_table <- matrix(NA, nrow = cnum, ncol = 2)

  for (i in seq_len(cnum)) {
    c_value <- call[i, 1]
    
    if (c_value <= pmax) {
      bounds <- do.call(c_function, c(c_value, c_function_args))
      c_table[i, 1] <- bounds$lower
      c_table[i, 2] <- bounds$upper
    } else {
      c_table[i, 1] <- c_table[1, 1]
      c_table[i, 2] <- c_table[1, 2]
    }
    
    if (verbose) cat(".")
  }
  
  # =========================================================================
  # 6.1 Monotonize the output
  # =========================================================================
  c_table[, 1] <- cummin(c_table[, 1])
  c_table[, 2] <- cummax(c_table[, 2])
  
  # =========================================================================
  # 7. Calculate breakdown point
  # =========================================================================
  if (!is.na(breakdown)) {
    if (verbose) cat("Calculating breakdown point...\n")
    breakdown_point <- calculate_breakdown(c_table, tol, breakdown, c_function, c_function_args)
  } else {
    breakdown_point <- NA
  }
  
  # =========================================================================
  # 8. Return results
  # =========================================================================
  
  result <- list(
    c_table = c_table,
    breakdown_point = breakdown_point,
    p0 = p0,
    p1 = p1,
    wsupp = wsupp,
    call = call
  )
  
  return(result)
}


# 関数の定義
cqte_bounds_c <- function(c, qtl, p0, p1, wsupp, bscoef0, bscoef1, varlist, data) {
  
  # =========================================================================
  # 1. Calculate the arguments to the quantile function given c, p1, p0
  # =========================================================================
  
  qtls <- cqtl_arg_c(c, p0, p1, qtl)
  
  # =========================================================================
  # 2. Calculate the upper and lower bound on the quantile function for each potential outcome
  # =========================================================================
  
  qtl_bounds <- matrix(NA, nrow = 4, ncol = 1)
  
  for (j in 1:4) {
    # Select the right covariates for each set of coefficients
    nobs <- nrow(wsupp)
    nwsupp <- ncol(wsupp) # Dimension of W, including constant
    
    # bscoef0 contains constant and X
    if (ncol(bscoef0) > (nwsupp + 1)) {
      if (j %% 2 == 0) {
        xw <- cbind(rep(0, nobs), wsupp[, 1:(nwsupp - 1)], rep(0, nobs) * wsupp[, 1:(nwsupp - 1)], rep(1, nobs))
      } else {
        xw <- cbind(rep(1, nobs), wsupp[, 1:(nwsupp - 1)], rep(1, nobs) * wsupp[, 1:(nwsupp - 1)], rep(1, nobs))
      }
    } else {
      if (j %% 2 == 0) {
        xw <- cbind(rep(0, nobs), wsupp)
      } else {
        xw <- cbind(rep(1, nobs), wsupp)
      }
    }
    
    # Calculate the quantile regression if not on the boundary
    # Use the approximate boundary coefficients if not
    qtl_value <- qtls[j, 1]
    if (qtl_value > 0 && qtl_value < 1) {
      qreg_model <- rq(as.formula(paste(varlist, collapse = " + ")), data = data, tau = qtl_value)
      qtl_coef <- coef(qreg_model)
    } else if (qtl_value == 0) {
      qtl_coef <- bscoef0
    } else if (qtl_value == 1) {
      qtl_coef <- bscoef1
    }
    
    # Calculate the quantile function value
    qtl_bounds[j, 1] <- xw %*% qtl_coef
  }
  
  # Take differences to calculate CQTE
  qtl_bounds <- matrix(qtl_bounds, nrow = 2)
  cqtl_bounds <- qtl_bounds[, 1] - qtl_bounds[, 2]
  
  lower <- cqtl_bounds[1, 1]
  upper <- cqtl_bounds[2, 1]
  
  return(list(lower = lower, upper = upper))
}


# 必要なパッケージの読み込み
library(mgcv)

# 関数の定義
cate_bounds_c <- function(c, p0, p1, wsupp, pchipcoefs, nvars, nobs) {
  
  # =========================================================================
  # 1. Convert inputs to matrices
  # =========================================================================
  
  p0 <- as.matrix(p0)
  p1 <- as.matrix(p1)
  wsupp <- as.matrix(wsupp)
  
  # =========================================================================
  # 2. Calculate bounds using cate_integral_c function
  # =========================================================================
  
  bds <- cate_integral_c(c, p0, p1, wsupp, pchipcoefs, nvars, 1, 0)
  
  lower <- bds[1, 1]
  upper <- bds[2, 1]
  
  return(list(lower = lower, upper = upper))
}


# 関数の定義
ate_bounds_c <- function(c, p0, p1, wsupp, pchipcoefs, nvars, nobs) {
  
  # =========================================================================
  # 1. Initialize bounds matrix
  # =========================================================================
  
  bds <- matrix(NA, nrow = 2, ncol = nobs)
  
  p0 <- as.matrix(p0)
  p1 <- as.matrix(p1)
  wsupp <- as.matrix(wsupp)
  
  # =========================================================================
  # 2. Calculate bounds for each observation
  # =========================================================================
  
  for (n in 1:nobs) {
    bds[, n] <- cate_integral_c(c, p0[n, 1], p1[n, 1], wsupp[n, ], pchipcoefs, nvars, 1, 0)
  }
  
  bounds <- colMeans(bds)
  
  lower <- bounds[1]
  upper <- bounds[2]
  
  return(list(lower = lower, upper = upper))
}



# 関数の定義
binary_ate_bounds_c <- function(c, p0, p1, p10, p11) {
  
  # =========================================================================
  # 1. Calculate bounds using binary_ate_c function
  # =========================================================================
  
  bds <- binary_ate_c(c, p0, p1, p10, p11)
  
  lower <- bds[1, 1]
  upper <- bds[1, 2]
  
  return(list(lower = lower, upper = upper))
}


# 関数の定義
atet_bounds_c <- function(c, p0, p1, wsupp, pchipcoefs, nvars, nobs, exp_y0, exp_y1, up0, up1) {
  
  # =========================================================================
  # 1. Initialize bounds matrix
  # =========================================================================
  
  bds <- matrix(NA, nrow = 2, ncol = nobs)
  
  p0 <- as.matrix(p0)
  p1 <- as.matrix(p1)
  wsupp <- as.matrix(wsupp)
  
  # =========================================================================
  # 2. Calculate bounds for each observation
  # =========================================================================
  
  for (n in 1:nobs) {
    bds[, n] <- cate_integral_c(c, p0[n, 1], p1[n, 1], wsupp[n, ], pchipcoefs, nvars, 0, 1)
  }
  
  bds <- colMeans(bds)
  
  lower <- exp_y1 - ((bds[2] - up0 * exp_y0) / up1)
  upper <- exp_y1 - ((bds[1] - up0 * exp_y0) / up1)
  
  return(list(lower = lower, upper = upper))
}


# 関数の定義
qte_bounds_c <- function(c, qtl, p0, p1, wsupp, tol, Y, X, qmvarlist, data) {
  
  # =========================================================================
  # 1. Initialize variables
  # =========================================================================
  
  nobs <- nrow(wsupp)
  nwsupp <- ncol(wsupp) # Dimension of W, including constant
  ncov <- length(qmvarlist) # Dimension of q(X,W), excluding X
  
  # Compare the dimensions to determine whether XW is included
  if (ncov > (nwsupp - 1)) {
    xw0 <- cbind(rep(0, nobs), wsupp[, 1:(nwsupp - 1)], rep(0, nobs) * wsupp[, 1:(nwsupp - 1)], rep(1, nobs))
    xw1 <- cbind(rep(1, nobs), wsupp[, 1:(nwsupp - 1)], rep(1, nobs) * wsupp[, 1:(nwsupp - 1)], rep(1, nobs))
  } else {
    xw0 <- cbind(rep(0, nobs), wsupp)
    xw1 <- cbind(rep(1, nobs), wsupp)
  }
  
  p0 <- as.matrix(p0)
  p1 <- as.matrix(p1)
  
  # =========================================================================
  # 2. Calculate LQ and UQ for QTE
  # =========================================================================
  
  LQ0 <- LQ(c, Y, X, qmvarlist, data, tol, qtl, 0)
  UQ1 <- UQ(c, Y, X, qmvarlist, data, tol, qtl, 1)
  UQTE <- UQ1 - LQ0
  
  LQ1 <- LQ(c, Y, X, qmvarlist, data, tol, qtl, 1)
  UQ0 <- UQ(c, Y, X, qmvarlist, data, tol, qtl, 0)
  LQTE <- LQ1 - UQ0
  
  # =========================================================================
  # 3. Clean up temporary variables
  # =========================================================================
  
  rm(xw0, xw1, p0, p1)
  
  # =========================================================================
  # 4. Return results
  # =========================================================================
  
  return(list(upper = UQTE, lower = LQTE))
}



# 関数の定義
calculate_breakdown <- function(c_table, tol, breakdown_y, function_name, arguments) {
  
  # =========================================================================
  # 1. Determine the side of the breakdown point
  # =========================================================================
  
  if (c_table[1, 2] > breakdown_y) {
    side <- "lower"
  } else {
    side <- "upper"
  }
  
  breakdown_start_points <- function(c_table, breakdown_y, side) {
    # Calculate the start points for the breakdown
    # ここに具体的なアルゴリズムを実装します
    return(list(clow = 0.5, chi = 0.5)) # 仮の結果
  }
  
  start_points <- breakdown_start_points(c_table, breakdown_y, side)
  clow <- start_points$clow
  chi <- start_points$chi
  
  # =========================================================================
  # 2. Check if breakdown point = 1 (i.e., no breakdown), otherwise run bisection algorithm
  # =========================================================================
  
  if (clow == 1) {
    output <- 1
  } else {
    bisection <- function(clow, chi, side, tol, breakdown_y, function_name, arguments) {
      # Bisection algorithm to find the breakdown point
      # ここに具体的なアルゴリズムを実装します
      return(list(output = 0.5)) # 仮の結果
    }
    
    result <- bisection(clow, chi, side, tol, breakdown_y, function_name, arguments)
    output <- result$output
  }
  
  return(output)
}



# 関数の定義
bisection <- function(clow, chi, side, tol, y, function_name, arguments) {
  cmid <- (chi + clow) / 2
  
  # 呼び出し関数の実行
  fcmid <- do.call(function_name, c(list(cmid), arguments))[[side]]
  
  if (abs(fcmid - y) < tol || abs((chi - clow) / 2) < tol) {
    output <- cmid
  } else {
    if (fcmid < y) {
      clow <- cmid
    } else {
      chi <- cmid
    }
    output <- bisection(clow, chi, side, tol, y, function_name, arguments)
  }
  
  return(output)
}


# 関数の定義
UF1 <- function(gridy, gridc, Y, X, qmvarlist, data) {
  # Estimate F_{Y|X,W}(y | x,w)
  data$indicatorY <- as.numeric(data[[Y]] <= gridy)
  logit_model <- glm(indicatorY ~ ., data = data[c(X, qmvarlist)], family = binomial, control = list(maxit = 2000))
  Fy1w <- plogis(predict(logit_model, newdata = data))
  UF1W <- pmin(((p1 * Fy1w) / (p1 - gridc)) * (p1 > gridc) + (p1 <= gridc), (p1 * Fy1w + gridc) / (p1 + gridc), p1 * Fy1w + (1 - p1))
  UF1 <- mean(UF1W)
  return(UF1)
}

UF0 <- function(gridy, gridc, Y, X, qmvarlist, data) {
  # Estimate F_{Y|X,W}(y | x,w)
  data$indicatorY <- as.numeric(data[[Y]] <= gridy)
  logit_model <- glm(indicatorY ~ ., data = data[c(X, qmvarlist)], family = binomial, control = list(maxit = 2000))
  Fy0w <- plogis(predict(logit_model, newdata = data))
  UF0W <- pmin(((p0 * Fy0w) / (p0 - gridc)) * (p0 > gridc) + (p0 <= gridc), (p0 * Fy0w + gridc) / (p0 + gridc), p0 * Fy0w + (1 - p0))
  UF0 <- mean(UF0W)
  return(UF0)
}

LF1 <- function(gridy, gridc, Y, X, qmvarlist, data) {
  # Estimate F_{Y|X,W}(y | x,w)
  data$indicatorY <- as.numeric(data[[Y]] <= gridy)
  logit_model <- glm(indicatorY ~ ., data = data[c(X, qmvarlist)], family = binomial, control = list(maxit = 2000))
  Fy1w <- plogis(predict(logit_model, newdata = data))
  LF1W <- pmax((p1 * Fy1w) / (p1 + gridc), ((p1 * Fy1w - gridc) / (p1 - gridc)) * (p1 > gridc), p1 * Fy1w)
  LF1 <- mean(LF1W)
  return(LF1)
}

LF0 <- function(gridy, gridc, Y, X, qmvarlist, data) {
  # Estimate F_{Y|X,W}(y | x,w)
  data$indicatorY <- as.numeric(data[[Y]] <= gridy)
  logit_model <- glm(indicatorY ~ ., data = data[c(X, qmvarlist)], family = binomial, control = list(maxit = 2000))
  Fy0w <- plogis(predict(logit_model, newdata = data))
  LF0W <- pmax((p0 * Fy0w) / (p0 + gridc), ((p0 * Fy0w - gridc) / (p0 - gridc)) * (p0 > gridc), p0 * Fy0w)
  LF0 <- mean(LF0W)
  return(LF0)
}


# 関数の定義
UQ <- function(gridc, Y, X, qmvarlist, data, tol, tau, treat) {
  # Estimate F_{Y|X,W}(y | x,w)
  ymax <- max(data[[Y]], na.rm = TRUE)
  ymin <- min(data[[Y]], na.rm = TRUE)
  
  LF_treat <- numeric(10)
  for (i in 2:9) {
    y_i <- (ymax - ymin) * ((i - 1) / 9) + ymin
    LF_treat[i] <- LF(treat, y_i, gridc, Y, X, qmvarlist, data)
  }
  
  y1 <- ymin
  y10 <- ymax
  LF_treat[1] <- 0 # h(y_min)=0
  LF_treat[10] <- 1 # h(y_max)=1
  
  vecLF <- LF_treat
  index <- which(vecLF < tau)
  ylow <- max(index)
  yhi <- ylow + 1
  
  yhi_val <- (ymax - ymin) * ((yhi - 1) / 9) + ymin
  ylow_val <- (ymax - ymin) * ((ylow - 1) / 9) + ymin
  
  output <- bisection(ylow_val, yhi_val, "bound", tol, tau, LF, list(treat, gridc, Y, X, qmvarlist, data))
  return(output)
}

LQ <- function(gridc, Y, X, qmvarlist, data, tol, tau, treat) {
  # Estimate F_{Y|X,W}(y | x,w)
  ymax <- max(data[[Y]], na.rm = TRUE)
  ymin <- min(data[[Y]], na.rm = TRUE)
  
  UF_treat <- numeric(10)
  for (i in 2:9) {
    y_i <- (ymax - ymin) * ((i - 1) / 9) + ymin
    UF_treat[i] <- UF(treat, y_i, gridc, Y, X, qmvarlist, data)
  }
  
  y1 <- ymin
  y10 <- ymax
  UF_treat[1] <- 0 # h(y_min)=0
  UF_treat[10] <- 1 # h(y_max)=1
  
  vecUF <- UF_treat
  index <- which(vecUF < tau)
  ylow <- max(index)
  yhi <- ylow + 1
  
  yhi_val <- (ymax - ymin) * ((yhi - 1) / 9) + ymin
  ylow_val <- (ymax - ymin) * ((ylow - 1) / 9) + ymin
  
  output <- bisection(ylow_val, yhi_val, "bound", tol, tau, UF, list(treat, gridc, Y, X, qmvarlist, data))
  return(output)
}


# 関数の定義


# 必要なパッケージの読み込み
library(matrixStats)

# FUNCTION: Binary ATE Bounds
binary_ate_c <- function(c, p0, p1, p10, p11) {
  UP1 <- rowMins(cbind(
    ((p11 * p1) / (p1 - c)) * (p1 > c) + (p1 <= c),
    ((p11 * p1 + c)) / (p1 + c),
    (p11 * p1) + (1 - p1)
  ))
  LP1 <- rowMaxs(cbind(
    (p11 * p1) / (p1 + c),
    (((p11 * p1) - c) / (p1 - c)) * (p1 > c),
    p11 * p1
  ))
  UP0 <- rowMins(cbind(
    ((p10 * p0) / (p0 - c)) * (p0 > c) + (p0 <= c),
    ((p10 * p0 + c)) / (p0 + c),
    (p10 * p0) + (1 - p0)
  ))
  LP0 <- rowMaxs(cbind(
    (p10 * p0) / (p0 + c),
    (((p10 * p0) - c) / (p0 - c)) * (p0 > c),
    p10 * p0
  ))
  return(colMeans(cbind(LP1 - UP0, UP1 - LP0)))
}

# FUNCTION: Conditional Quantile Argument
cqtl_arg_c <- function(c, p0, p1, qtl) {
  limits <- qtl_bound(c, p0, p1)
  coef <- qtl_coef(c, p0, p1)
  
  qtl_select <- (limits[, 1] <= qtl) & (qtl < limits[, 2])
  coef <- coef[qtl_select, ]
  
  return(coef * c(1, qtl))
}

# FUNCTION: CATE integral, single c dependence
cate_integral_c <- function(c, p0, p1, cov, pchipcoefs, nvars, cate, catt) {
  int_segs <- matrix(NA, nrow = 8, ncol = 1)
  qtl_coef <- matrix(NA, nrow = 4, ncol = nvars)
  n <- length(cov)
  
  if (nvars > (n + 2)) {
    control_cov <- c(0, cov[1:(n - 1)], 0 * cov[1:(n - 1)], 1)
    treat_cov <- c(1, cov[1:(n - 1)], 1 * cov[1:(n - 1)], 1)
  } else {
    control_cov <- c(0, cov)
    treat_cov <- c(1, cov)
  }
  
  limits <- qtl_bound(c, p0, p1)
  coef <- qtl_coef(c, p0, p1)
  
  # calculate the integrals
  # outer loop: covariates
  # inner loop: intervals of integrals
  # TODO: is there a more efficient way to do this to avoid inner loop?
  for (k in 1:nvars) {
    for (i in 1:8) {
      int_segs[i, 1] <- analyint(limits[i, 1], limits[i, 2], coef[i, 2], coef[i, 1], pchipcoefs[[k]])
    }
    qtl_coef[, k] <- rowSums(matrix(int_segs, nrow = 4, byrow = TRUE))
  }
  
  if (cate) {
    cate_bounds <- numeric(2)
    cate_bounds[1] <- sum(qtl_coef[1, ] * treat_cov - qtl_coef[2, ] * control_cov)
    cate_bounds[2] <- sum(qtl_coef[3, ] * treat_cov - qtl_coef[4, ] * control_cov)
    return(cate_bounds)
  }
  
  if (catt) {
    exp_y0_bounds <- numeric(2)
    exp_y0_bounds[1] <- sum(qtl_coef[4, ] * control_cov)
    exp_y0_bounds[2] <- sum(qtl_coef[2, ] * control_cov)
    return(exp_y0_bounds)
  }
}

# FUNCTION: Quantile Coefficients
qtl_coef <- function(c, p0, p1) {
  if ((c < p0) & (c < p1)) {
    coef <- matrix(c(
      0, 1 - c / p1, -c / p1, 1 + c / p1,
      0, 1 + c / p0, c / p0, 1 - c / p0,
      0, 1 + c / p1, c / p1, 1 - c / p1,
      0, 1 - c / p0, -c / p0, 1 + c / p0
    ), ncol = 2, byrow = TRUE)
  } else if ((p0 <= c) & (c < p1)) {
    coef <- matrix(c(
      0, 1 - c / p1, 1 - 1 / p1, 1 / p1,
      0, 1 + c / p0, 1, 0,
      0, 1 / p1, c / p1, 1 - c / p1,
      0, 0, -c / p0, 1 + c / p0
    ), ncol = 2, byrow = TRUE)
  } else if ((p1 <= c) & (c < p0)) {
    coef <- matrix(c(
      0, 0, -c / p1, 1 + c / p1,
      0, 1 / p0, c / p0, 1 - c / p0,
      0, 1 + c / p1, 1, 0,
      0, 1 - c / p0, 1 - 1 / p0, 1 / p0
    ), ncol = 2, byrow = TRUE)
  } else if ((p0 <= c) & (p1 <= c)) {
    coef <- matrix(c(
      0, 0, 1 - 1 / p1, 1 / p1,
      0, 1 / p0, 1, 0,
      0, 1 / p1, 1, 0,
      0, 0, 1 - 1 / p0, 1 / p0
    ), ncol = 2, byrow = TRUE)
  }
  return(coef)
}

# FUNCTION: Quantile Bounds
qtl_bound <- function(c, p0, p1) {
  if ((c < p0) & (c < p1)) {
    bound1 <- 0.5
    bound2 <- 0.5
    bound3 <- 0.5
    bound4 <- 0.5
  } else if ((p0 <= c) & (c < p1)) {
    bound1 <- (1 - p1) / (1 - (p1 - c))
    bound2 <- 1 / (1 + (c / p0))
    bound3 <- c / (1 - (p1 - c))
    bound4 <- (c / p0) / (1 + c / p0)
  } else if ((p1 <= c) & (c < p0)) {
    bound1 <- (c / p1) / (1 + c / p1)
    bound2 <- c / (1 - (p0 - c))
    bound3 <- 1 / (1 + (c / p1))
    bound4 <- (1 - p0) / (1 - (p0 - c))
  } else if ((p0 <= c) & (p1 <= c)) {
    bound1 <- p0
    bound2 <- p0
    bound3 <- p1
    bound4 <- p1
  }
  
  bounds <- matrix(c(0, bound1,
                     bound1, 1,
                     0, bound2,
                     bound2, 1,
                     0, bound3,
                     bound3, 1,
                     0, bound4,
                     bound4, 1), 
                   ncol = 2, byrow = TRUE)
  
  return(bounds)
}

# FUNCTION: PCHIP All Coefficients
pchip_all_coef <- function(xnodes, ynodes, nvars) {
  pchipcoefs <- vector("list", nvars)
  
  for (k in 1:nvars) {
    pchipcoefs[[k]] <- pchip(xnodes, ynodes[, k])
  }
  
  coef1 <- sapply(pchipcoefs, function(coef) pchipval(coef, 1))
  coef0 <- sapply(pchipcoefs, function(coef) pchipval(coef, 0))
  
  xnodes <- c(0, xnodes, 1)
  ynodes <- rbind(coef0, ynodes, coef1)
  
  for (k in 1:nvars) {
    pchipcoefs[[k]] <- pchip(xnodes, ynodes[, k])
  }
  
  return(pchipcoefs)
}

# FUNCTION: PCHIP Coefficients
pchip <- function(x, y) {
  n <- length(x)
  h <- diff(x)
  delta <- diff(y) / h
  d <- pchipslopes(h, delta)
  
  c <- (3 * delta - 2 * d[1:(n - 1)] - d[2:n]) / h
  b <- (d[1:(n - 1)] - 2 * delta + d[2:n]) / (h^2)
  
  pchipcoef <- cbind(b, c, d[1:(n - 1)], y[1:(n - 1)], x[1:(n - 1)])
  return(pchipcoef)
}

# FUNCTION: PCHIP slopes
pchipslopes <- function(h, delta) {
  n <- length(h) + 1
  d <- numeric(n)
  k <- 1 + which(sign(delta[1:(n - 2)]) * sign(delta[2:(n - 1)]) > 0)
  w1 <- 2 * h[k] + h[k - 1]
  w2 <- h[k] + 2 * h[k - 1]
  d[k] <- (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])
  d[1] <- pchipend(h[1], h[2], delta[1], delta[2])
  d[n] <- pchipend(h[n - 1], h[n - 2], delta[n - 1], delta[n - 2])
  
  return(d)
}

# FUNCTION: PCHIP end
pchipend <- function(h1, h2, del1, del2) {
  d <- ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2)
  if (sign(d) != sign(del1)) d <- 0
  else if (sign(del1) != sign(del2) && abs(d) > abs(3 * del1)) d <- 3 * del1
  
  return(d)
}

# FUNCTION: PCHIP values
pchipval <- function(coef, u) {
  nu <- length(u)
  n <- nrow(coef)
  k <- rep(1, nu)
  b <- coef[1:(n - 1), 1]
  c <- coef[1:(n - 1), 2]
  d <- coef[1:(n - 1), 3]
  y <- coef[1:(n - 1), 4]
  x <- coef[, 5]
  
  if (nu == 1) {
    for (j in 2:(n - 1)) {
      if (x[j] <= u) k <- j
    }
  } else {
    for (j in 2:(n - 1)) {
      which <- which(x[j] <= u)
      k[which] <- j
    }
  }
  s <- u - x[k]
  return(y[k] + s * (d[k] + s * (c[k] + s * b[k])))
}

# FUNCTION: Analytical Integral
analyint <- function(a, b, c, d, pchipcoef) {
  if (c == 0) {
    # TODO: avoid recalculating each time...
    return((b - a) * pchipval(pchipcoef, d))
  }
  
  # この部分で引っかかっている
  
  lower_bound <- c * a + d
  upper_bound <- c * b + d
  difflower <- matrix(lower_bound, nrow = nrow(pchipcoef), ncol = 1) - pchipcoef[, 5]
  diffupper <- matrix(upper_bound, nrow = nrow(pchipcoef), ncol = 1) - pchipcoef[, 5]
  indexlower <- which(difflower <= 0)[1]
  indexupper <- tail(which(diffupper >= 0), 1)
  
  if (indexlower == indexupper) {
    midint <- 0
  } else {
    yk1 <- pchipcoef[(indexlower + 1):indexupper, 4]
    yk <- pchipcoef[indexlower:(indexupper - 1), 4]
    hk <- pchipcoef[(indexlower + 1):indexupper, 5] - pchipcoef[indexlower:(indexupper - 1), 5]
    dk1 <- pchipcoef[(indexlower + 1):indexupper, 3]
    dk <- pchipcoef[indexlower:(indexupper - 1), 3]
    integral <- (1 / 2) * (yk1 + yk) * hk - (1 / 12) * (dk1 - dk) * (hk^2)
    midint <- sum(integral)
  }
  
  if (indexlower == 1) {
    leftint <- 0
  } else {
    y1 <- pchipcoef[indexlower, 4]
    y0 <- pchipcoef[indexlower - 1, 4]
    d1 <- pchipcoef[indexlower, 3]
    d0 <- pchipcoef[indexlower - 1, 3]
    h0 <- pchipcoef[indexlower, 5] - pchipcoef[indexlower - 1, 5]
    x0 <- pchipcoef[indexlower - 1, 5]
    u0 <- lower_bound - x0
    piece1 <- (y1 / (h0^2)) * u0^3 - (y1 / (2 * h0^3)) * u0^4 + y0 * u0 - (y0 / h0^2) * u0^3 + (y0 / (2 * h0^3)) * u0^4
    piece2 <- (d1 / h0^2) * ((1 / 4) * u0^4 - (h0 / 3) * u0^3) + (d0 / h0^2) * ((1 / 4) * u0^4 - (2 * h0 / 3) * u0^3 + ((h0^2) / 2) * u0^2)
    leftint <- (1 / c) * (piece1 + piece2)
  }
  
  if (indexupper == nrow(pchipcoef)) {
    rightint <- 0
  } else {
    y1 <- pchipcoef[indexupper + 1, 4]
    y0 <- pchipcoef[indexupper, 4]
    d1 <- pchipcoef[indexupper + 1, 3]
    d0 <- pchipcoef[indexupper, 3]
    h0 <- pchipcoef[indexupper + 1, 5] - pchipcoef[indexupper, 5]
    x0 <- pchipcoef[indexupper, 5]
    u0 <- upper_bound - x0
    piece1 <- (y1 / (h0^2)) * u0^3 - (y1 / (2 * h0^3)) * u0^4 + y0 * u0 - (y0 / h0^2) * u0^3 + (y0 / (2 * h0^3)) * u0^4
    piece2 <- (d1 / h0^2) * ((1 / 4) * u0^4 - (h0 / 3) * u0^3) + (d0 / h0^2) * ((1 / 4) * u0^4 - (2 * h0 / 3) * u0^3 + ((h0^2) / 2) * u0^2)
    rightint <- (1 / c) * (piece1 + piece2)
  }
  
  analyint <- (1 / c) * (leftint + midint + rightint)
  return(analyint)
}

# FUNCTION: breakdown_start_points
breakdown_start_points <- function(c_table, bd_y, side) {
  if (side == "lower") {
    index <- which(c_table[, 2] > bd_y)
    c1 <- "chi"
    c2 <- "clow"
  } else if (side == "upper") {
    index <- which(c_table[, 3] < bd_y)
    c1 <- "clow"
    c2 <- "chi"
  }
  
  if (length(index) == nrow(c_table)) {
    clow <- 1
  } else {
    clim <- max(index)
    clow <- c_table[clim, 1]
    chi <- c_table[clim + 1, 1]
  }
  
  return(list(clow = clow, chi = chi))
}

# FUNCTION: minimum 2
min2 <- function(A, B) {
  D <- pmin(A, B)
  return(D)
}

# FUNCTION: minimum 3
min3 <- function(A, B, C) {
  D <- pmin(A, B)
  E <- pmin(B, C)
  F <- pmin(D, E)
  return(F)
}

# FUNCTION: maximum 2
max2 <- function(A, B) {
  D <- pmax(A, B)
  return(D)
}

# FUNCTION: maximum 3
max3 <- function(A, B, C) {
  D <- pmax(A, B)
  E <- pmax(B, C)
  F <- pmax(D, E)
  return(F)
}

apply_cumfunc <- function(matname, col, f) {
  mat <- get(matname)
  for (row in 1:nrow(mat)) {
    mat[row, col] <- f(mat[1:row, col])
  }
  assign(matname, mat, envir = .GlobalEnv)
}

# FUNCTION: Sort Stata matrix
sort_st_matrix <- function(matname, sortcol) {
  mat <- get(matname)
  perm <- order(mat[, sortcol])
  sorted_mat <- mat[perm, ]
  assign(matname, sorted_mat, envir = .GlobalEnv)
}

# FUNCTION: Matrix Equal
mat_eq <- function(a, b, out) {
  result <- all.equal(get(a), get(b))
  assign(out, as.numeric(result == TRUE), envir = .GlobalEnv)
}



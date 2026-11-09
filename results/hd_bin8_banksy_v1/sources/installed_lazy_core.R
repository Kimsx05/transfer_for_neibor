function (data_own, knn_list = NULL, group_knn = NULL, group_idx = NULL, 
    lambda, npcs = 50L, split_scale = FALSE, scale_max = 10, 
    seed = NULL, pca_backend = c("cpp", "r"), verbose = TRUE) 
{
    n_genes <- nrow(data_own)
    n_cells <- ncol(data_own)
    gene_names <- rownames(data_own)
    cell_names <- colnames(data_own)
    lambdas <- getLambdas(lambda, n_harmonics = 1)
    has_groups <- !is.null(group_knn)
    max_npcs <- min(2L * n_genes, n_cells) - 1L
    if (npcs >= max_npcs) {
        stop("npcs (", npcs, ") must be < min(2*n_genes, n_cells) = ", 
            max_npcs + 1L)
    }
    if (min(2L * n_genes, n_cells) < 6L) {
        stop("lazy PCA requires at least 3 genes and 6 cells")
    }
    if (has_groups) 
        .checkGroupIdx(group_idx)
    if (has_groups) {
        if (verbose) 
            message("Building per-group weight matrices")
        W_list <- lapply(seq_along(group_idx), function(gr) .buildWeightMatrix(group_knn[[gr]][[1]], 
            length(group_idx[[gr]])))
        rm(group_knn)
        gc(verbose = FALSE)
        if (verbose) 
            message("Building per-group expression matrices")
        gcm_list <- lapply(group_idx, function(cid) as(data_own[, 
            cid], "dgCMatrix"))
        own_result <- .lazy_own_scaling(NULL, split_scale, group_idx, 
            groups = NULL, ugroups = NULL, n_genes, n_cells, 
            scale_max, verbose, gcm_list = gcm_list)
        h0_result <- .lazy_h0_scaling_grouped(NULL, W_list, split_scale, 
            group_idx, n_genes, n_cells, scale_max, verbose, 
            gcm_list = gcm_list)
        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) 
                length(own_result$excess@x)
            else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) 
                length(h0_result$excess@x)
            else 0
            message("Clipping corrections: own=", n_own_clip, 
                " H0=", n_h0_clip, " entries")
        }
        w_csc <- lapply(W_list, function(w) if (is(w, "dgCMatrix")) 
            w
        else as(w, "dgCMatrix"))
        w_slots <- list(i = lapply(w_csc, slot, "i"), p = lapply(w_csc, 
            slot, "p"), x = lapply(w_csc, slot, "x"), ncol = vapply(w_csc, 
            ncol, integer(1)))
        rm(w_csc)
        gcm_slots <- list(i = lapply(gcm_list, slot, "i"), p = lapply(gcm_list, 
            slot, "p"), x = lapply(gcm_list, slot, "x"))
        rm(gcm_list, data_own)
        gc(verbose = FALSE)
        banksy_op <- structure(list(gcm_slots = gcm_slots, W_list = W_list, 
            w_slots = w_slots, mu = list(own_result$mu, h0_result$mu), 
            sd = list(own_result$sd, h0_result$sd), lam = lambdas, 
            excess = list(own_result$excess, h0_result$excess), 
            valid = list(own_result$valid, h0_result$valid), 
            has_groups = TRUE, split_scale = split_scale, group_idx = group_idx, 
            n_genes = n_genes, n_cells = n_cells), class = "BanksyLazy")
        rm(own_result, h0_result)
    }
    else {
        if (verbose) 
            message("Building sparse weight matrix")
        W <- .buildWeightMatrix(knn_list[[1]], n_cells)
        rm(knn_list)
        gc(verbose = FALSE)
        own_result <- .lazy_own_scaling(data_own, split_scale, 
            group_idx, groups = NULL, ugroups = NULL, n_genes, 
            n_cells, scale_max, verbose)
        h0_result <- .lazy_h0_scaling(data_own, W, split_scale, 
            group_idx, n_genes, n_cells, scale_max, verbose)
        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) 
                length(own_result$excess@x)
            else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) 
                length(h0_result$excess@x)
            else 0
            message("Clipping corrections: own=", n_own_clip, 
                " H0=", n_h0_clip, " entries")
        }
        banksy_op <- structure(list(gcm = data_own, W = W, mu = list(own_result$mu, 
            h0_result$mu), sd = list(own_result$sd, h0_result$sd), 
            lam = lambdas, excess = list(own_result$excess, h0_result$excess), 
            valid = list(own_result$valid, h0_result$valid), 
            has_groups = FALSE, split_scale = split_scale, group_idx = group_idx, 
            n_genes = n_genes, n_cells = n_cells), class = "BanksyLazy")
        rm(own_result, h0_result)
    }
    gc(verbose = FALSE)
    verbose.seed(seed)
    pca_backend <- match.arg(pca_backend, c("cpp", "r"))
    if (pca_backend == "cpp") {
        work <- as.integer(min(npcs + 7L, n_cells - 1L))
        maxit <- 1000L
        tol <- 1e-05
        svtol <- tol
        if (verbose) 
            message("Computing BANKSY PCA (", npcs, " PCs) via ", 
                "C++ irlba (work=", work, ")")
        m <- 2L * n_genes
        n <- n_cells
        k <- npcs
        eps <- .Machine$double.eps
        eps23 <- eps^(2/3)
        sqrteps <- sqrt(eps)
        V_buf <- numeric(as.double(n) * work)
        W_mat <- matrix(0, m, work)
        F_vec <- numeric(n)
        B <- NULL
        q1 <- .lazy_start_vector(n, seed = seed)
        lanczos_set_col(V_buf, 0L, q1, n)
        rm(q1)
        t0 <- proc.time()["elapsed"]
        mprod <- 0L
        iter <- 0L
        Smax <- 1
        Smin <- NULL
        lastsv <- numeric(0)
        converged <- FALSE
        restart <- FALSE
        while (iter < maxit) {
            iter <- iter + 1L
            j <- if (iter == 1L && !restart) 
                1L
            else k + 1L
            VJ <- lanczos_get_col(V_buf, j - 1L, n)
            avj <- as.numeric(.banksy_lazy_mult(banksy_op, matrix(VJ, 
                ncol = 1)))
            rm(VJ)
            W_mat[, j] <- avj
            mprod <- mprod + 1L
            if (iter > 1L && j > 1L) 
                W_mat[, j] <- W_mat[, j] - as.numeric(W_mat[, 
                  1:(j - 1), drop = FALSE] %*% crossprod(W_mat[, 
                  1:(j - 1), drop = FALSE], W_mat[, j]))
            S <- sqrt(sum(W_mat[, j]^2))
            if (S < eps23) {
                W_mat[, j] <- .lazy_start_vector(m, j + 1L, seed)
                if (j > 1L) 
                  W_mat[, j] <- W_mat[, j] - as.numeric(W_mat[, 
                    1:(j - 1), drop = FALSE] %*% crossprod(W_mat[, 
                    1:(j - 1), drop = FALSE], W_mat[, j]))
                W_mat[, j] <- W_mat[, j]/sqrt(sum(W_mat[, j]^2))
                S <- 0
            }
            else {
                W_mat[, j] <- W_mat[, j]/S
            }
            while (j <= work) {
                F_vec <- as.numeric(.banksy_lazy_mult(banksy_op, 
                  matrix(W_mat[, j], ncol = 1), transpose = TRUE))
                mprod <- mprod + 1L
                F_vec <- F_vec - S * lanczos_get_col(V_buf, j - 
                  1L, n)
                F_vec <- lanczos_reorth(V_buf, j, F_vec, n)
                if (j + 1L <= work) {
                  R <- lanczos_norm(F_vec)
                  if (R < eps23) {
                    F_vec <- .lazy_start_vector(n, j + 1L, seed)
                    F_vec <- lanczos_reorth(V_buf, j, F_vec, 
                      n)
                    lanczos_set_col(V_buf, j, F_vec/lanczos_norm(F_vec), 
                      n)
                    R <- 0
                  }
                  else {
                    lanczos_set_col(V_buf, j, F_vec/R, n)
                  }
                  if (is.null(B)) {
                    B <- matrix(c(S, R), nrow = 1)
                  }
                  else {
                    B <- rbind(cbind(B, 0), c(rep(0, ncol(B) - 
                      1), S, R))
                  }
                  VJP1 <- lanczos_get_col(V_buf, j, n)
                  W_mat[, j + 1L] <- as.numeric(.banksy_lazy_mult(banksy_op, 
                    matrix(VJP1, ncol = 1)))
                  rm(VJP1)
                  mprod <- mprod + 1L
                  W_mat[, j + 1L] <- W_mat[, j + 1L] - R * W_mat[, 
                    j]
                  if (j > 1L) 
                    W_mat[, j + 1L] <- W_mat[, j + 1L] - as.numeric(W_mat[, 
                      1:j, drop = FALSE] %*% crossprod(W_mat[, 
                      1:j, drop = FALSE], W_mat[, j + 1L]))
                  S <- sqrt(sum(W_mat[, j + 1L]^2))
                  if (S < eps23) {
                    W_mat[, j + 1L] <- .lazy_start_vector(m, 
                      j + 2L, seed)
                    if (j > 0L) 
                      W_mat[, j + 1L] <- W_mat[, j + 1L] - as.numeric(W_mat[, 
                        1:j, drop = FALSE] %*% crossprod(W_mat[, 
                        1:j, drop = FALSE], W_mat[, j + 1L]))
                    W_mat[, j + 1L] <- W_mat[, j + 1L]/sqrt(sum(W_mat[, 
                      j + 1L]^2))
                    S <- 0
                  }
                  else {
                    W_mat[, j + 1L] <- W_mat[, j + 1L]/S
                  }
                }
                else {
                  B <- rbind(B, c(rep(0, j - 1), S))
                }
                j <- j + 1L
            }
            Bsz <- nrow(B)
            R_F <- lanczos_norm(F_vec)
            F_vec <- F_vec/R_F
            Bsvd <- svd(B)
            Smax <- max(Smax, Bsvd$d[1])
            R_vals <- R_F * Bsvd$u[Bsz, , drop = FALSE]
            ct_conv <- all(abs(R_vals[1:k]) < tol * Smax)
            sv_conv <- length(lastsv) >= k && all(abs(Bsvd$d[1:k] - 
                lastsv[1:k]) < svtol * Bsvd$d[1:k])
            lastsv <- Bsvd$d
            if (verbose && (iter <= 2L || iter%%5 == 0 || ct_conv || 
                sv_conv)) {
                elapsed <- proc.time()["elapsed"] - t0
                message(sprintf("  iter=%d  mprod=%d  sv[%d]=%.4e  t=%.0fs", 
                  iter, mprod, k, Bsvd$d[k], elapsed))
            }
            if (ct_conv || sv_conv) {
                converged <- TRUE
                if (verbose) 
                  message(sprintf("  Converged: iter=%d, mprod=%d", 
                    iter, mprod))
                break
            }
            V_rot <- lanczos_extract(V_buf, Bsvd$v[, 1:k, drop = FALSE], 
                n, Bsz, k)
            for (i in seq_len(k)) lanczos_set_col(V_buf, i - 
                1L, V_rot[, i], n)
            rm(V_rot)
            lanczos_set_col(V_buf, k, F_vec, n)
            B <- cbind(diag(Bsvd$d[1:k], nrow = k), R_F * Bsvd$u[Bsz, 
                1:k])
            W_mat[, 1:k] <- W_mat[, 1:Bsz, drop = FALSE] %*% 
                Bsvd$u[, 1:k, drop = FALSE]
            restart <- TRUE
        }
        if (!converged) 
            warning("C++ irlba did not converge; try increasing maxit or work")
        Bsvd <- svd(B)
        V <- lanczos_extract(V_buf, Bsvd$v[, 1:npcs, drop = FALSE], 
            n, nrow(B), npcs)
        rm(V_buf)
        gc(verbose = FALSE)
        U <- W_mat[, 1:nrow(B), drop = FALSE] %*% Bsvd$u[, 1:npcs, 
            drop = FALSE]
        rm(W_mat)
        d <- Bsvd$d[1:npcs]
        rm(Bsvd, B, F_vec)
        pca <- list(u = U, v = V, d = d)
    }
    else {
        if (verbose) 
            message("Computing BANKSY PCA (", npcs, " PCs) via irlba")
        .irlba_iter <- 0L
        .irlba_t0 <- proc.time()["elapsed"]
        .instrumented_mult <- function(A, x, transpose = FALSE) {
            .irlba_iter <<- .irlba_iter + 1L
            if (.irlba_iter%%20 == 0) 
                gc(verbose = FALSE)
            if (.irlba_iter%%20 == 1) {
                rss <- tryCatch({
                  l <- readLines("/proc/self/status", warn = FALSE)
                  v <- grep("^VmRSS:", l, value = TRUE)
                  as.numeric(gsub("[^0-9]", "", v))/1024^2
                }, error = function(e) NA)
                elapsed <- proc.time()["elapsed"] - .irlba_t0
                message(sprintf("  irlba iter=%d  t=%.0fs  RSS=%.1fGB", 
                  .irlba_iter, elapsed, rss))
            }
            .banksy_lazy_mult(A, x, transpose)
        }
        pca <- irlba::irlba(banksy_op, nv = npcs, mult = .instrumented_mult, 
            v = .lazy_start_vector(n_cells, seed = seed))
    }
    embeddings <- sweep(pca$v, 2, pca$d, `*`)
    rownames(embeddings) <- cell_names
    colnames(embeddings) <- paste0("PC", seq_len(npcs))
    loadings <- pca$u
    feat_names <- c(gene_names, paste0(gene_names, ".m0"))
    rownames(loadings) <- feat_names
    colnames(loadings) <- paste0("PC", seq_len(npcs))
    total_var <- sum(pca$d^2)
    stdev <- pca$d/sqrt(max(1, n_cells - 1))
    list(embeddings = embeddings, loadings = loadings, stdev = stdev, 
        total_var = total_var)
}

prepare_observations <- function(data.mtx)
{
  data.vec <- as.vector(data.mtx)
  data.vec[ data.vec == 0 ] <- NA
  data.obs.mask <- !is.na( data.vec )
  data.idx <- ifelse( data.obs.mask, cumsum( data.obs.mask ), -cumsum( !data.obs.mask ) )
  data.idx.arr <- matrix( data.idx, ncol = ncol( data.mtx ) )
  list( Nobserved = sum(data.obs.mask),
        oData = array( data.vec[ data.obs.mask ] ),
        dataIndex = as.array(data.idx.arr)
  )
}

obj2effect <- function(expXeff, obj_class, objs, data2obj, data2exp) {
  effXexp <- t(expXeff)
  expXeff_mask <- apply(effXexp, c(2, 1), function(x) x != 0.0)
  objXexp_mask <- matrix(FALSE, nrow = length(objs), ncol = ncol(effXexp))
  if (length(data2obj)>0) { for (i in 1:length(data2obj)) {
    objXexp_mask[data2obj[[i]], data2exp[[i]]] <- TRUE
  } }
  objXeff_mask <- objXexp_mask %*% expXeff_mask
  dnams <- list(objs, dimnames(effXexp)[[1]])
  names(dnams) <- c(obj_class, names(dimnames(effXexp))[[1]])
  dimnames(objXeff_mask) <- dnams
  if (length(objXeff_mask) > 0L) {
    res <- as.data.frame(as.table(t(objXeff_mask))) %>% dplyr::filter(Freq != 0) %>% dplyr::select(-Freq)
  } else {
    res <- data.frame() %>% dplyr::mutate(x = character(), y = character())
    names(res) <- c(obj_class, names(dnams)[[2]])
  }
  res[,paste(colnames(res), collapse='_')] <- seq_len(nrow(res))
  res
}

effects_cumsum <- function(df, group_col) {
  n_effects <- df %>% dplyr::arrange_(group_col) %>% dplyr::group_by_(group_col) %>%
    dplyr::summarize(n_effects = n()) %>% .$n_effects
  if (length(n_effects)>0) {
    res <- cumsum(as.integer(n_effects))
    c(1L, res+1L)
  } else {
    return(0L)
  }
}

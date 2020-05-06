

# functions ---------------------------------------------------------------


maximize_ess <- function (x) {
  if (length(x)<1) return (length(x))
  sample_every <- 1
  if (length(x)>1500) sample_every <- floor(length(x)/1500)
  x <- x[seq(1, length(x), by=sample_every)]
  num_tests <- min(1000, length(x)-2)
  burnin_num <- round(seq(2, (length(x)-num_tests), length.out=num_tests))
  ess <- mcmapply(seq, list(-1), as.list(-burnin_num), SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
    mcmapply(`[`, list(x), ., SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
    lapply(effectiveSize) %>% unlist()
  c(burnin=burnin_num[which.max(ess)]*sample_every, ess=max(ess))
}

get_num_causing <- function (R, k, causing=.8) {
  prob_seq <- seq(0, 0.999, length.out=1000)
  if (length(R)==1) {
    qvalues <- qnbinom(prob_seq, size=k, mu=R)
    1-rev(prob_seq)[min(which(cumsum(rev(qvalues))>=(causing*sum(qvalues))))]
  }
  threadx <- rep(1:detectCores(), each=ceiling(length(R)/detectCores()))[1:length(R)]
  mcmapply(function (Rx, kx) {
    lapply(1:length(Rx), function (i) qnbinom(prob_seq, size=kx[i], mu=Rx[i])) %>% 
      lapply(function (x) 1-rev(prob_seq)[min(which(cumsum(rev(x))>=(causing*sum(x))))]) %>%
      unlist()
  }, split(R, threadx), split(k, threadx), SIMPLIFY=FALSE, mc.cores=detectCores()) %>%
    unlist()
}


get_change_point <- function(input_df, change_points, dt) {
  input_df <- filter(input_df, !is.na(value))
  input_locations <- change_points[names(loc_dict)[match(input_df$location, loc_dict)]] %>% 
    lapply(`*`, dt) %>%
    mapply(`+`, input_df$start_date, .) %>%
    mapply(c, as.list(input_df$start_date), .)
  input_change_dates <- mcmapply(`[`, input_locations, gsub("[^\\d]+", "", input_df$name, perl=TRUE) %>% 
                                   as.numeric() %>% `+`(1) %>% as.list(), mc.cores=detectCores())
  input_df %>% 
    mutate(change_point=input_change_dates) %>%
    mutate(change_point_date=as.Date(date_decimal(change_point)))
}

save_plot_to_file <- function (plot_obj, file_prefix, ext, figwidth, figheight, fig_dir="figures") {
  lapply(ext, function (ext) {
    plot_obj %>%
      ggsave(file.path(fig_dir, paste0(file_prefix, ".", ext)), .,
             width=figwidth, height=figheight)
  })
}



align_traj_to_data <- function (infected, reported_data, sum_every_num) {
  multinom_prob <- reported_data
  if (sum(reported_data)==0) {
    multinom_prob <- reported_data+1
  } else if (max(reported_data)<5) {
    multinom_prob <- reported_data+max(reported_data)
  }
  if (sum(infected>=sum(reported_data))>20) infected <- infected[infected>=sum(reported_data)]
  outdf <- lapply(infected, function (x) {
    output <- reported_data
    if (x==0) return (output)
    if (x<sum(reported_data)) rmultinom(1, x, prob=multinom_prob)
    else c(output + rmultinom(1, x-sum(reported_data), prob=multinom_prob))
  })
  sum_every_num <- min(length(outdf[[1]]), sum_every_num)
  outdf %>%
    lapply(sum_every, sum_every_num) %>%
    bind_cols() %>%
    apply(1, quantile, c(0.5, 0.25, 0.75, 0.025, 0.975)) %>%
    t() %>%
    data.frame() %>%
    setNames(quantile_cols) %>%
    mutate(data=sum_every(reported_data, sum_every_num))
}



library(ggplot2)
library(data.table)
mean_vec <- c(10,20,5,25)
sim_fun_list <- list(
  constant_changes=function(N){
    rep(mean_vec, each=N/4)
  },
  linear_changes=function(N){
    rep(rep(mean_vec, each=10), l=N)
  })
lapply(sim_fun_list, function(f)f(40))

N_data_vec <- c(40, 400)
sim_data_list <- list()
sim_changes_list <- list()
sim_segs_list <- list()
for(N_data in N_data_vec){
  for(simulation in names(sim_fun_list)){
    sim_fun <- sim_fun_list[[simulation]]
    data_mean_vec <- sim_fun(N_data)
    end <- which(diff(data_mean_vec) != 0)
    set.seed(1)
    data_value <- rnorm(N_data, data_mean_vec, 2)
    cum.vec <- c(0, cumsum(data_value))
    n.segs <- sum(diff(data_mean_vec)!=0)+1
    wfit <- fpopw::Fpsn(data_value, n.segs)
    end <- wfit$t.est[n.segs, 1:n.segs]
    start <- c(1, end[-length(end)]+1)
    sim_segs_list[[paste(N_data, simulation)]] <- data.table(
      N_data, simulation, start.pos=start-0.5, end.pos=end+0.5,
      mean=(cum.vec[end+1]-cum.vec[start])/(end-start+1))
    sim_data_list[[paste(N_data, simulation)]] <- data.table(
      N_data, simulation, data_i=seq_along(data_value), data_value)
    sim_changes_list[[paste(N_data, simulation)]] <- data.table(
      N_data, simulation, end)
  }
}
addSim <- function(DT)DT[, Simulation := paste0("\n", simulation)][]
(sim_changes <- addSim(rbindlist(sim_changes_list)))
one_sim <- sim_data_list[["40 constant_changes"]]
gg <- ggplot()+
  theme_bw()+
  theme(text=element_text(size=14))+
  geom_point(aes(
    data_i, data_value),
    color="grey50",
    data=one_sim)+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=10))
gg
wfit <- fpopw::Fpop(one_sim$data_value, 50)
cum.vec <- c(0, cumsum(one_sim$data_value))
end <- wfit$t.est
start <- c(1, end[-length(end)]+1)
one_sim_means <- data.table(
  start.pos=start-0.5, end.pos=end+0.5,
  mean=(cum.vec[end+1]-cum.vec[start])/(end-start+1))
gg+
  geom_vline(aes(
    xintercept=start.pos),
    color="green",
    linetype="dashed",
    data=one_sim_means[-1])+
  geom_segment(aes(
    start.pos, mean,
    xend=end.pos, yend=mean),
    color="green",
    data=one_sim_means)
(sim_data <- addSim(rbindlist(sim_data_list)))
gg <- ggplot()+
  theme_bw()+
  theme(text=element_text(size=14))+
  geom_point(aes(
    data_i, data_value),
    color="grey50",
    data=sim_data)+
  facet_grid(Simulation ~ N_data, labeller=label_both, scales="free_x", space="free")+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=20))
gg
gg+
  geom_vline(aes(
    xintercept=end+0.5),
    data=sim_changes)
sim_segs <- addSim(rbindlist(sim_segs_list))
gg+
  geom_vline(aes(
    xintercept=start.pos),
    color="green",
    linetype="dashed",
    size=1,
    data=sim_segs[1<start.pos])+
  geom_segment(aes(
    start.pos, mean,
    xend=end.pos, yend=mean),
    data=sim_segs,
    color="green",
    size=2)
PELT <- function(sim.mat, penalty, prune=TRUE, verbose=FALSE){
  if(!is.matrix(sim.mat))sim.mat <- cbind(sim.mat)
  cum.data <- rbind(0, apply(sim.mat, 2, cumsum))
  sum_trick <- function(m, start, end)m[end+1,,drop=FALSE]-m[start,,drop=FALSE]
  cost_trick <- function(start, end){
    sum_vec <- sum_trick(cum.data, start, end)
    N <- end+1-start
    -sum_vec^2/N
  }
  N.data <- nrow(sim.mat)
  pelt.change.vec <- rep(NA_integer_, N.data)
  pelt.cost.vec <- rep(NA_real_, N.data+1)
  pelt.cost.vec[1] <- -penalty
  pelt.candidates.vec <- rep(NA_integer_, N.data)
  candidate.vec <- 1L
  if(verbose)index_dt_list <- list()
  for(up.to in 1:N.data){
    N.cand <- length(candidate.vec)
    pelt.candidates.vec[up.to] <- N.cand
    last.seg.cost <- rowSums(cost_trick(candidate.vec, rep(up.to, N.cand)))
    prev.cost <- pelt.cost.vec[candidate.vec]
    cost.no.penalty <- prev.cost+last.seg.cost
    total.cost <- cost.no.penalty+penalty
    if(verbose)index_dt_list[[up.to]] <- data.table(
      data_i=up.to, change_i=candidate.vec, cost=total.cost/up.to)
    best.i <- which.min(total.cost)
    pelt.change.vec[up.to] <- candidate.vec[best.i]
    total.cost.best <- total.cost[best.i]
    pelt.cost.vec[up.to+1] <- total.cost.best
    keep <- if(isTRUE(prune))cost.no.penalty < total.cost.best else TRUE
    candidate.vec <- c(candidate.vec[keep], up.to+1L)
  }
  list(
    change=pelt.change.vec,
    cost=pelt.cost.vec,
    candidates=pelt.candidates.vec,
    index_dt=if(verbose)rbindlist(index_dt_list))
}
decode <- function(best.change){
  seg.dt.list <- list()
  last.i <- length(best.change)
  while(last.i>0){
    first.i <- best.change[last.i]
    seg.dt.list[[paste(last.i)]] <- data.table(
      first.i, last.i)
    last.i <- first.i-1L
  }
  rbindlist(seg.dt.list)[seq(.N,1)]
}
pelt_info_list <- list()
pelt_segs_list <- list()
both_index_list <- list()
penalty <- 100
algo.prune <- c(
  OPART=FALSE,
  PELT=TRUE)
for(N_data in N_data_vec){
  for(simulation in names(sim_fun_list)){
    N_sim <- paste(N_data, simulation)
    data_value <- sim_data_list[[N_sim]]$data_value
    for(algo in names(algo.prune)){
      prune <- algo.prune[[algo]]
      fit <- PELT(data_value, penalty=penalty, prune=prune, verbose = TRUE)
      both_index_list[[paste(N_data, simulation, algo)]] <- data.table(
        N_data, simulation, algo, fit$index_dt)
      fit_seg_dt <- decode(fit$change)
      pelt_info_list[[paste(N_data,simulation,algo)]] <- data.table(
        N_data,simulation,algo,
        candidates=fit$candidates,
        cost=fit$cost[-1],
        data_i=seq_along(data_value))
      pelt_segs_list[[paste(N_data,simulation,algo)]] <- data.table(
        N_data,simulation,algo,
        fit_seg_dt)
    }
  }
}
(pelt_info <- addSim(rbindlist(pelt_info_list)))
(pelt_segs <- addSim(rbindlist(pelt_segs_list)))

algo.colors <- c(
  OPART="grey50",
  PELT="red",
  FPOP="blue",
  DUST="deepskyblue")
cat(paste0("\\definecolor{", names(algo.colors), "}{HTML}{", sub("#", "", animint2::toRGB(algo.colors)), "}\n"))
ggplot()+
  theme_bw()+
  theme(
    legend.position=c(0.8, 0.2),
    panel.spacing=grid::unit(1,"lines"),
    text=element_text(size=15))+
  geom_vline(aes(
    xintercept=end+0.5),
    data=sim_changes)+
  scale_color_manual(
    breaks=names(algo.colors),
    values=algo.colors)+
  geom_point(aes(
    data_i, candidates, color=algo),
    data=pelt_info)+
  facet_grid(Simulation ~ N_data, labeller=label_both, scales="free_x", space="free")+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=20))+
  scale_y_continuous(
    breaks=c(10, 100, 400))+
  theme(panel.grid.minor=element_blank())+
  coord_cartesian(expand=FALSE)

remotes::install_github("tdhock/fpopw/fpopw")
fpop_info_list <- list()
fpop_segs_list <- list()
for(N_data in N_data_vec){
  for(simulation in names(sim_fun_list)){
    N_sim <- paste(N_data, simulation)
    data_value <- sim_data_list[[N_sim]]$data_value
    pfit <- fpopw::Fpop(
      data_value, penalty, verbose_file=tempfile())
    both_index_list[[paste(N_data, simulation, "FPOP")]] <- data.table(
      N_data, simulation, algo="FPOP", pfit$model
    )[, let(
      data_i=data_i+1L,
      change_i=ifelse(change_i == -10, 0, change_i)+1
    )][]
    end <- pfit$t.est
    start <- c(1, end[-length(end)]+1)
    count_dt <- pfit$model[, .(
      intervals=.N,
      cost=min(cost)
    ), by=data_i]
    fpop_info_list[[paste(N_data,simulation)]] <- data.table(
      N_data,simulation,
      algo="FPOP",
      candidates=count_dt$intervals,
      cost=count_dt$cost,
      data_i=seq_along(data_value))
    fpop_segs_list[[paste(N_data,simulation)]] <- data.table(
      N_data,simulation,start,end)
  }
}
(fpop_info <- addSim(rbindlist(fpop_info_list)))
(fpop_segs <- rbindlist(fpop_segs_list))
both_info <- rbind(pelt_info, fpop_info)
ggplot()+
  theme_bw()+
  theme(
    legend.position=c(0.8, 0.2),
    panel.spacing=grid::unit(1,"lines"),
    text=element_text(size=15))+
  geom_vline(aes(
    xintercept=end+0.5),
    data=sim_changes)+
  scale_color_manual(
    breaks=names(algo.colors),
    values=algo.colors)+
  geom_point(aes(
    data_i, candidates, color=algo),
    data=both_info)+
  facet_grid(Simulation ~ N_data, labeller=label_both, scales="free_x", space="free")+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=20))+
  scale_y_continuous(
    breaks=c(10, 100, 400))+
  theme(panel.grid.minor=element_blank())+
  coord_cartesian(expand=FALSE)

both_info[data_i==max(data_i)][order(simulation, algo)]

remotes::install_github("vrunge/dust@910f8c67f99354fdb5ff7740e6436eb487d9efa6")
set.seed(1)
ex_data <- rnorm(7)
ex_penalty <- 1
dfit <- dust::dust.1D(ex_data, ex_penalty)
wfit <- fpopw::Fpop(ex_data, ex_penalty)
pfit <- PELT(ex_data, ex_penalty)
list(PELT=pfit$change, FPOP=wfit$path, DUST=dfit$changepoints)
rbind(
  PELT=pfit$cost[-1]/seq_along(ex_data),
  FPOP=(wfit$cost-ex_penalty)/seq_along(ex_data),
  DUST=2*dfit$costQ/seq_along(ex_data))

dust_info_list <- list()
dust_segs_list <- list()
for(N_data in N_data_vec){
  for(simulation in names(sim_fun_list)){
    N_sim <- paste(N_data, simulation)
    data_value <- sim_data_list[[N_sim]]$data_value
    dfit <- dust::dust.object.1D("gauss")
    candidates <- rep(NA_integer_, N_data)
    for(data_i in seq_along(data_value)){
      dfit$append_c(data_value[[data_i]], penalty)
      dfit$update_partition()
      pinfo <- dfit$get_partition()
      change_i <- pinfo$lastIndexSet[-1]+1L
      candidates[[data_i]] <- length(change_i)
      mean_cost <- pinfo$costQ/seq_along(pinfo$costQ)
      both_index_list[[paste(
        N_data, simulation, "DUST", data_i
      )]] <- data.table(
        N_data, simulation, algo="DUST", data_i,
        change_i, cost=mean_cost[change_i])
    }
    end <- pinfo$changepoints
    start <- c(1, end[-length(end)]+1)
    dust_info_list[[paste(N_data,simulation)]] <- data.table(
      N_data,simulation,
      algo="DUST",
      candidates,
      cost=pinfo$costQ,
      data_i=seq_along(data_value))
    dust_segs_list[[paste(N_data,simulation)]] <- data.table(
      N_data,simulation,start,end)
  }
}
(dust_info <- addSim(rbindlist(dust_info_list)))
(dust_segs <- rbindlist(dust_segs_list))
three_info <- rbind(pelt_info, fpop_info, dust_info)
ggplot()+
  theme_bw()+
  theme(
    legend.position=c(0.8, 0.2),
    panel.spacing=grid::unit(1,"lines"),
    text=element_text(size=15))+
  geom_vline(aes(
    xintercept=end+0.5),
    data=sim_changes)+
  scale_color_manual(
    breaks=names(algo.colors),
    values=algo.colors)+
  geom_point(aes(
    data_i, candidates, color=algo),
    data=three_info)+
  facet_grid(Simulation ~ N_data, labeller=label_both, scales="free_x", space="free")+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=20))+
  scale_y_continuous(
    breaks=c(10, 100, 400))+
  coord_cartesian(expand=FALSE)

algo.levs <- c("OPART","PELT","FPOP","DUST")
(both_index <- addSim(rbindlist(both_index_list, use.names=TRUE))[, let(
  Algorithm = factor(algo, algo.levs)
)][])
for(a in algo.levs){
  gg <- ggplot()+
    ggtitle(a)+
    theme_bw()+
    theme(
      text=element_text(size=20))+
    coord_equal()+
    geom_tile(aes(
      data_i, change_i, fill=cost),
      data=both_index[Algorithm==a])+
    scale_fill_gradient(low="black", high="red")+
    facet_grid(Simulation ~ N_data, label=label_both)
  print(gg)
}

ggplot()+
  theme_bw()+
  theme(
    legend.position=c(0.8, 0.2),
    text=element_text(size=15))+
  geom_tile(aes(
    data_i, change_i, fill=Algorithm),
    data=both_index)+
  scale_fill_manual(values=algo.colors)+
  scale_x_continuous(
    breaks=seq(0,max(N_data_vec),by=20))+
  facet_grid(Simulation ~ N_data, label=label_both, scales="free", space="free")

N_data <- 1e5
sim_fun <- sim_fun_list$constant_changes
data_mean_vec <- sim_fun(N_data)
set.seed(1)
data_value <- rnorm(N_data, data_mean_vec, 2)
wfit <- fpopw::Fpop(data_value, penalty)
dfit <- dust::dust.1D(data_value, penalty)
rbind(
  FPOP=wfit$t.est,
  DUST=dfit$changepoints)

base_N <- c(100,200,400,800)
(all_N <- unlist(lapply(10^seq(0,5), function(x)x*base_N)))

grid_args <- list(
  list(simulation=names(sim_fun_list)),
  DUST=quote({
    dfit <- dust::dust.1D(data_list[[simulation]], penalty)
    with(dfit, data.frame(
      mean_candidates=mean(nb),
      segments=length(changepoints),
      max_seg_size=max(diff(c(0,changepoints)))))
  }),
  FPOP=quote({
    pfit <- fpopw::Fpop(
      data_list[[simulation]], penalty, verbose_file=tempfile())
    count_dt <- pfit$model[, .(intervals=.N), by=data_i]
    data.frame(
      mean_candidates=mean(count_dt$intervals),
      segments=length(pfit$t.est),
      max_seg_size=max(diff(c(0,pfit$t.est))))
  }))
for(algo in names(algo.prune)){
  prune <- algo.prune[[algo]]
  grid_args[[algo]] <- substitute({
    data_value <- data_list[[simulation]]
    fit <- PELT(data_value, penalty, prune=prune)
    dec <- decode(fit$change)
    data.frame(
      mean_candidates=mean(fit$candidates),
      segments=nrow(dec),
      max_seg_size=max(diff(c(0,dec$last.i))))
  }, list(prune=prune))
}
(expr.list <- do.call(atime::atime_grid, grid_args))

cache.rds <- "2025-04-15-PELT-vs-fpopw.rds"
if(file.exists(cache.rds)){
  atime_list <- readRDS(cache.rds)
}else{
  atime_list <- atime::atime(
    N=all_N,
    setup={
      data_list <- list()
      for(simulation in names(sim_fun_list)){
        sim_fun <- sim_fun_list[[simulation]]
        data_mean_vec <- sim_fun(N)
        set.seed(1)
        data_list[[simulation]] <- rnorm(N, data_mean_vec, 2)
      }
    },
    expr.list=expr.list,
    seconds.limit=1,
    result=TRUE)
  saveRDS(atime_list, cache.rds)
}

refs_list <- atime::references_best(atime_list)
ggplot()+
  theme(text=element_text(size=12))+
  directlabels::geom_dl(aes(
    N, empirical, color=expr.grid, label=expr.grid),
    data=refs_list$measurements,
    method="right.polygons")+
  geom_line(aes(
    N, empirical, color=expr.grid),
    data=refs_list$measurements)+
  scale_color_manual(
    "algorithm",
    guide="none",
    values=algo.colors)+
  facet_grid(unit ~ simulation, scales="free")+
  scale_x_log10(
    "N = number of data in sequence",
    breaks=10^seq(2,6),
    limits=c(NA, 1e8))+
  scale_y_log10("")

refs_list$meas[, let(
  Simulation = sub("_","\n",simulation),
  Algorithm = factor(expr.grid, algo.levs)
)][]
ggplot()+
  theme(
    legend.position="none",
    text=element_text(size=12))+
  directlabels::geom_dl(aes(
    N, empirical, color=Simulation, label=Simulation),
    data=refs_list$measurements,
    method="right.polygons")+
  geom_line(aes(
    N, empirical, color=Simulation),
    data=refs_list$measurements)+
  facet_grid(unit ~ Algorithm, scales="free")+
  scale_x_log10(
    "N = number of data in sequence",
    breaks=10^seq(2,6),
    limits=c(NA, 1e8))+
  scale_y_log10("")

expr.levs <- CJ(
  algo=algo.levs,
  simulation=names(sim_fun_list)
)[algo.levs, paste0(algo,"\n",simulation), on="algo"]
edit_expr <- function(DT)DT[
, expr.name := factor(sub(" simulation=", "\n", expr.name), expr.levs)]
edit_expr(refs_list$meas)
edit_expr(refs_list$plot.ref)
plot(refs_list)+
  theme(panel.spacing=grid::unit(1.5,"lines"))

pred_list <- predict(refs_list)
ggplot()+
  theme_bw()+
  theme(text=element_text(size=20))+
  geom_line(aes(
    N, empirical, color=expr.grid),
    data=pred_list$measurements[unit=="seconds"])+
  directlabels::geom_dl(aes(
    N, unit.value, color=expr.grid, label=sprintf(
      "%s\n%s\nN=%s",
      expr.grid,
      ifelse(expr.grid %in% c("FPOP","DUST"), "C++", "R"),
      format(round(N), big.mark=",", scientific=FALSE, trim=TRUE))),
    data=pred_list$prediction,
    method=list(cex=1.5,directlabels::polygon.method("top",offset.cm = 0.5)))+
  scale_color_manual(
    "algorithm",
    guide="none",
    values=algo.colors)+
  facet_grid(. ~ simulation, scales="free", labeller=label_both)+
  scale_x_log10(
    "N = number of data in sequence",
    breaks=10^seq(2,7),
    limits=c(NA, 5e7))+
  scale_y_log10(
    "Computation time (seconds)",
    breaks=10^seq(-3,0),
    limits=10^c(-3,1))

cand_dt <- pred_list$measurements[unit=="mean_candidates"]
ggplot()+
  theme_bw()+
  theme(text=element_text(size=20))+
  geom_line(aes(
    N, empirical, color=expr.grid),
    data=cand_dt)+
  directlabels::geom_dl(aes(
    N, empirical, color=expr.grid, label=expr.grid),
    data=cand_dt,
    method=list(cex=2, "right.polygons"))+
  scale_color_manual(
    "algorithm",
    guide="none",
    values=algo.colors)+
  facet_grid(. ~ simulation, scales="free", labeller=label_both)+
  scale_x_log10(
    "N = number of data in sequence",
    breaks=10^seq(2,7),
    limits=c(NA, 1e8))+
  scale_y_log10(
    "Mean number of candidate change-points",
    limits=10^c(0, 4))

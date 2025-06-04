```{r}
mean_vec <- c(10,15,8)
data_mean_vec <- rep(mean_vec, each=20)
library(data.table)
N_data_vec <- c(40, 400)
N_data <- 60
end <- which(diff(data_mean_vec) != 0)
set.seed(3)
data_value <- rnorm(N_data, data_mean_vec, 2)
one_sim <- data.table(data_i=seq_along(data_value), data_value)
```
library(animint2)
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
cum.vec <- c(0, cumsum(one_sim$data_value))
(fpop_means <- data.table(penalty=10^seq(0, 3, by=0.2))[, {
  wfit <- fpopw::Fpop(one_sim$data_value, penalty)
  end <- wfit$t.est
  start <- c(1, end[-length(end)]+1)
  data.table(
    start.pos=start-0.5, end.pos=end+0.5,
    mean=(cum.vec[end+1]-cum.vec[start])/(end-start+1))
}, by=.(log10.penalty=log10(penalty))][])
fpop_changes <- fpop_means[, {
  diff.mean <- diff(mean)
  rbind(
    data.table(variable="changes", value=sum(diff.mean!=0)),
    data.table(variable="L1norm", value=sum(abs(diff.mean))))
}, by=log10.penalty][variable=="changes"]
animint(
  data=ggplot()+
    geom_point(aes(
      data_i, data_value),
      data=one_sim)+
    geom_vline(aes(
      xintercept=start.pos,
      key=paste(log10.penalty, start.pos)),
      color='green',
      size=1,
      linetype="dashed",
      showSelected='log10.penalty',
      data=fpop_means[1<start.pos])+
    geom_text(aes(
      x=0, y=4, label=sprintf("penalty=%.2f", 10^log10.penalty)),
      data=fpop_changes,
      hjust=0,
      showSelected='log10.penalty',
      size=15)+
    geom_segment(aes(
      start.pos, mean,
      key=paste(log10.penalty, start.pos),
      xend=end.pos, yend=mean),
      color='green',
      size=2,
      showSelected='log10.penalty',
      data=fpop_means),
  overview=ggplot()+
    scale_y_continuous("changes")+
    geom_text(aes(
      log10.penalty, value+0.5, label=value),
      data=fpop_changes)+
    geom_point(aes(
      log10.penalty, value),
      data=fpop_changes)+
    make_tallrect(fpop_changes, "log10.penalty"),
  time=list(
    variable='log10.penalty',
    ms=400),
  source="https://github.com/tdhock/functional-pruning-theory/"
)
```

```{r sim-data-model}

Kmax <- 15
wfit <- fpopw::Fpsn(one_sim$data_value, Kmax)
(fpsn_means <- data.table(segments=1:Kmax)[,{
  end <- wfit$t.est[segments, 1:segments]
  start <- c(1, end[-length(end)]+1)
  data.table(
    start.pos=start-0.5, end.pos=end+0.5,
    mean=(cum.vec[end+1]-cum.vec[start])/(end-start+1))
}, by=segments][])
fpsn_changes <- fpsn_means[, {
  diff.mean <- diff(mean)
  rbind(
    data.table(variable="changes", value=sum(diff.mean!=0)),
    data.table(variable="L1norm", value=sum(abs(diff.mean))))
}, by=segments][variable=="changes"]
animint(
  data=ggplot()+
    geom_point(aes(
      data_i, data_value),
      data=one_sim)+
    geom_vline(aes(
      xintercept=start.pos,
      key=paste(segments, start.pos)),
      color='green',
      size=1,
      linetype="dashed",
      showSelected='segments',
      data=fpsn_means[1<start.pos])+
    geom_segment(aes(
      start.pos, mean,
      key=paste(segments, start.pos),
      xend=end.pos, yend=mean),
      color='green',
      size=2,
      showSelected='segments',
      data=fpsn_means),
  overview=ggplot()+
    scale_y_continuous("changes")+
    geom_text(aes(
      segments, value+0.5, label=value),
      data=fpsn_changes)+
    geom_point(aes(
      segments, value),
      data=fpsn_changes)+
    make_tallrect(fpsn_changes, "segments"),
  time=list(
    variable='segments',
    ms=400)
)

flsa_mean_dt <- data.table(penalty=10^seq(1, 2, by=0.025))[, data.table(
  mean=as.numeric(flsa::flsa(one_sim$data_value, lambda2=penalty)),
  data_i=1:nrow(one_sim)
), by=.(segments=log10(penalty))]
flsa_mean_changes <- flsa_mean_dt[, {
  diff.mean <- diff(mean)
  rbind(
    data.table(variable="changes", value=sum(diff.mean!=0)),
    data.table(variable="L1norm", value=sum(abs(diff.mean))))
}, by=segments]
(flsa_segs_dt <- flsa_mean_dt[, {
  is.diff <- diff(mean)!=0
  start <- 1+c(0, which(is.diff))
  end <- c(start[-1]-1, .N)
  data.table(start.pos=start-0.5,end.pos=end+0.5,mean=mean[start])
}, by=segments][])
animint(
  data=ggplot()+
    geom_point(aes(
      data_i, data_value),
      data=one_sim)+
    geom_vline(aes(
      xintercept=start.pos,
      key=paste(segments, start.pos)),
      color='green',
      size=1,
      linetype="dashed",
      showSelected='segments',
      data=flsa_segs_dt[1<start.pos])+
    geom_segment(aes(
      start.pos, mean,
      key=paste(segments, start.pos),
      xend=end.pos, yend=mean),
      color='green',
      size=2,
      showSelected='segments',
      data=flsa_segs_dt),
  overview=ggplot()+
    theme_animint(width=600)+
    scale_y_continuous("")+
    geom_text(aes(
      segments, value+0.3, label=value),
      data=flsa_mean_changes[variable=="changes"])+
    facet_grid(variable ~ ., scales="free")+
    geom_point(aes(
      segments, value),
      data=flsa_mean_changes)+
    scale_x_continuous(breaks=seq(1, 2, by=0.1))+
    make_tallrect(flsa_mean_changes, "segments"),
  time=list(
    variable='segments',
    ms=400)
)
    


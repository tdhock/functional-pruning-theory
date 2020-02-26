data(Mono27ac, package="PeakSegDisk", envir=environment())
library(data.table)
loss.list <- list()
N.data.vec <- 10^seq(1, 3)
for(penalty in c(0, 1e2, 1e4, 1e6)){
  for(N.data in N.data.vec){
    Mono27ac$coverage[, real := count]
    Mono27ac$coverage[, synthetic := 1:.N]
    for(short.type in c("real", "synthetic")){
      some.cov <- data.table(Mono27ac$coverage)
      some.cov$count <- some.cov[[short.type]]
      for(offset in 0:9){
        df <- some.cov[(1:N.data)+offset, .(chrom, chromStart, chromEnd, count)]
        fit <- PeakSegDisk::PeakSegFPOP_df(df, penalty)
        loss.list[[paste(penalty, N.data, short.type, offset)]] <- data.table(
          N.data,
          offset,
          short.type,
          fit$loss)
      }
    }
  }
}
(loss <- do.call(rbind, loss.list))[, .(
  penalty, short.type, N.data,
  mean.intervals, max.intervals,
  megabytes, seconds)][order(penalty, short.type, N.data)]

(worst.dt <- data.table(
  N.data=N.data.vec,
  short.type="theoretical",
  mean.intervals=(N.data.vec+1)/2))
measured.dt <- loss[, .(
  mean=mean(mean.intervals),
  sd=sd(mean.intervals),
  count=.N
), by=.(penalty, short.type, N.data)]

library(ggplot2)
one <- function(short.type, data.type, color){
  data.table(short.type, data.type, color)
}
type.dt <- rbind(
  one("theoretical", "Theoretical\nworst case", "grey50"),
  one("synthetic", "Synthetic\nIncreasing", "red"),
  one("real", "Real ChIP-seq", "blue"))
measured.types <- type.dt[measured.dt, on=list(short.type)]
worst.types <- type.dt[worst.dt, on=list(short.type)]
(type.colors <- type.dt[, structure(color, names=data.type)])
gg <- ggplot()+
  guides(
    color=guide_legend(keyheight=3)
  )+
  geom_blank(aes(
    N.data, 1),
    data=data.table(N.data=c(5, 2000)))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ penalty, labeller=label_both)+
  geom_line(aes(
    N.data, mean.intervals, color=data.type),
    size=1,
    data=worst.types)+
  scale_color_manual(
    values=type.colors,
    breaks=names(type.colors))+
  scale_fill_manual(
    values=type.colors,
    guide=FALSE,
    breaks=names(type.colors))+
  geom_ribbon(aes(
    N.data, ymin=mean-sd, ymax=mean+sd, fill=data.type),
    alpha=0.5,
    data=measured.types)+
  geom_line(aes(
    N.data, mean, color=data.type),
    size=1,
    data=measured.types)+
  scale_x_log10("N data")+
  scale_y_log10("Mean intervals (candidate changepoints)")
pdf("figure-worst-case.pdf", width=7, height=4)
print(gg)
dev.off()

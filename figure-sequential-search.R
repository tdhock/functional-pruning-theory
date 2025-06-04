data(Mono27ac, package="PeakSegDisk")
library(data.table)
Mono27ac$coverage
library(ggplot2)
ggplot()+
  theme_bw()+
  geom_step(aes(
    chromStart/1e3, count),
    color="grey50",
    data=Mono27ac$coverage)
data.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11-60000-580000")
dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
future::plan("multisession")
fit10 <- PeakSegDisk::sequentialSearch_dir(data.dir, 10L, verbose=1)
plot(fit10)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=Mono27ac$coverage)

target.min.peaks <- 100
target.max.peaks <- 100
pen.vec <- c(0,Inf)
loss.dt.list <- list()
for(pen in pen.vec){
  fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, pen)
  loss.dt.list[[paste(pen)]] <- fit$loss[, .(penalty, peaks, total.loss)]
}
(loss.dt <- rbindlist(loss.dt.list))
get_selection <- function(loss){
  selection.df <- penaltyLearning::modelSelection(loss, "total.loss", "peaks")
  with(selection.df, data.table(
    total.loss,
    min.lambda,
    penalty,
    max.lambda,
    peaks_after=c(peaks[-1],NA),
    peaks
  ))[
  , helpful := peaks_after<target.max.peaks & peaks>target.min.peaks &
      peaks_after+1 < peaks & !max.lambda %in% names(loss.dt.list)
  ][]
}
(selection.dt <- get_selection(loss.dt))
get_gg <- function(){
  all.dot.dt <- rbind(
    data.table(loss.dt[, .(
      penalty, peaks, total.loss, helpful=FALSE)], dot="evaluated"),
    data.table(selection.dt[-.N, .(
      penalty=max.lambda, peaks, total.loss, helpful)], dot="candidate"))
  dot.dt <- all.dot.dt[order(penalty)][{
    i <- which(if(any(helpful))helpful==TRUE else peaks==target.min.peaks)
    seq(max(1, i-2), min(.N, i+2))
  }]
  ggplot()+
    geom_abline(aes(
      slope=peaks,
      intercept=total.loss),
      data=loss.dt)+
    geom_point(aes(
      penalty,
      total.loss+ifelse(peaks==0, 0, penalty*peaks),
      color=helpful),
      size=5,
      data=dot.dt)+
    geom_point(aes(
      penalty,
      total.loss+ifelse(peaks==0, 0, penalty*peaks),
      color=helpful,
      fill=dot),
      shape=21,
      size=4,
      data=dot.dt)+
    geom_label(aes(
      penalty,
      total.loss+ifelse(peaks==0, 0, penalty*peaks),
      hjust=fcase(
        penalty==Inf, 1,
        penalty==0, -0.5,
        default=0.5),
      vjust=fcase(
        penalty==0, 0.5,
        default=1.5),
      label=peaks),
      alpha=0.5,
      data=dot.dt[dot=="evaluated"])+
    scale_fill_manual(values=c(evaluated="white", candidate="black"))+
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"))+
    scale_y_continuous("Cost = Loss + penalty * (model size)")
}
get_gg()
iteration <- 1
slide.list <- character()
while({
  iteration <- iteration+1
  it.png <- sprintf("figure-sequential-search-%d.png", iteration)
  png(it.png, width=6, height=4, units="in", res=200)
  print(get_gg()+ggtitle(paste0(
    "Goal size=", target.min.peaks,
    ", iteration=", iteration,
    ", model sizes computed so far:\n",
    paste(loss.dt$peaks, collapse=", "))))
  dev.off()
  slide.list[[paste(iteration)]] <- sprintf("
\\begin{frame}
  \\includegraphics[width=\\textwidth]{%s}
\\end{frame}
", it.png)
  sum(selection.dt$helpful)==1
}){
  pen <- selection.dt[helpful==TRUE, max.lambda]
  fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, pen)
  print(loss.dt.list[[paste(pen)]] <- fit$loss[, .(penalty, peaks, total.loss)])
  loss.dt <- rbindlist(loss.dt.list)
  selection.dt <- get_selection(loss.dt)
}
cat(slide.list, file="figure-sequential-search.tex")

pen.list <- list()
seq_it <- function(){
  pen.vec <- selection.dt[helpful==TRUE, max.lambda]
  pen.list[[length(pen.list)+1]] <<- pen.vec
  loss.dt.list[paste(pen.vec)] <<- future.apply::future_lapply(pen.vec, function(pen){
    fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, pen)
    fit$loss[, .(penalty, peaks, total.loss)]
  })
  loss.dt <<- rbindlist(loss.dt.list)
  selection.dt <<- get_selection(loss.dt)
  get_gg()
}
seq_it()

seq_it()

while(any(selection.dt$helpful)){
  gg <- seq_it()
}
gg
selection.dt

central_launch <- function(target.min.peaks, target.max.peaks, LAPPLY=future.apply::future_lapply){
  data.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11-60000-580000")
  dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
  write.table(
    Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  pen.vec <- c(0,Inf)
  loss.dt.list <- list()
  cand.dt.list <- list()
  iteration <- 1L
  while(length(pen.vec)){
    loss.dt.list[paste(pen.vec)] <- LAPPLY(
      pen.vec, function(pen){
        start.time <- Sys.time()
        fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, pen)
        fit$loss[, .(
          penalty, peaks, total.loss, iteration,
          process=factor(Sys.getpid()), start.time, end.time=Sys.time())]
      }
    )
    start.time <- Sys.time()
    loss.dt <- rbindlist(loss.dt.list)
    selection.df <- penaltyLearning::modelSelection(loss.dt, "total.loss", "peaks")
    selection.dt <- with(selection.df, data.table(
      total.loss,
      min.lambda,
      penalty,
      max.lambda,
      peaks_after=c(peaks[-1],NA),
      peaks
    ))[
    , helpful := peaks_after<target.max.peaks & peaks>target.min.peaks &
        peaks_after+1 < peaks & !max.lambda %in% names(loss.dt.list)
    ][]
    pen.vec <- selection.dt[helpful==TRUE, max.lambda]
    cand.dt.list[[length(cand.dt.list)+1]] <- data.table(
      iteration,
      process=factor(Sys.getpid()), start.time, end.time=Sys.time())
    iteration <- iteration+1L
  }
  list(loss=loss.dt, candidates=rbindlist(cand.dt.list))
}
(loss_5_15 <- central_launch(5,15))

viz_workers <- function(L){
  min.time <- min(L$loss$start.time)
  for(N in names(L)){
    L[[N]] <- data.table(L[[N]], computation=N)
    for(pos in c("start", "end")){
      time.vec <- L[[N]][[paste0(pos,".time")]]
      set(L[[N]], j=paste0(pos,".seconds"), value=time.vec-min.time)
    }
  }
  ##L$cand$process <- L$loss[, .(count=.N), by=process][order(-count)][1, process]
  biggest.diff <- L$loss[,as.numeric(diff(c(min(start.seconds),max(end.seconds))))]
  n.proc <- length(unique(L$loss$process))
  best.seconds <- n.proc*biggest.diff
  total.seconds <- L$loss[,sum(end.seconds-start.seconds)]
  efficiency <- total.seconds/best.seconds
  seg.dt <- rbind(L$cand,L$loss[,names(L$cand),with=FALSE])
  proc.levs <- unique(c(levels(L$cand$process), levels(L$loss$process)))
  text.dt <- L$loss[,.(
    mid.seconds=as.POSIXct((as.numeric(min(start.seconds))+as.numeric(max(end.seconds)))/2),
    label=paste(peaks,collapse=",")
  ),by=.(process,iteration)]
  it.dt <- seg.dt[, .(min.seconds=min(start.seconds), max.seconds=max(end.seconds)), by=iteration]
  comp.colors <- c(
    loss="black",
    candidates="red")
  ggplot()+
    theme_bw()+
    ggtitle(sprintf("Worker efficiency = %.1f%%", 100*total.seconds/best.seconds))+
    geom_rect(aes(
      xmin=min.seconds, xmax=max.seconds,
      ymin=-Inf, ymax=Inf),
      fill="grey",
      color="white",
      alpha=0.5,
      data=it.dt)+
    geom_text(aes(
      x=(min.seconds+max.seconds)/2, Inf,
      hjust=ifelse(iteration==1, 1, 0.5),
      label=paste0(ifelse(iteration==1, "it=", ""), iteration)),
      data=it.dt,
      vjust=1)+
    directlabels::geom_dl(aes(
      start.seconds, process,
      label.group=paste(iteration, process),
      label=peaks),
      data=L$loss,
      method=polygon.mine("bottom", offset.cm=0.2, custom.colors=list(box.color="white")))+
    geom_segment(aes(
      start.seconds, process,
      color=computation,
      xend=end.seconds, yend=process),
      data=seg.dt)+
    geom_point(aes(
      start.seconds, process,
      color=computation),
      shape=1,
      data=seg.dt)+
    scale_color_manual(
      values=comp.colors,
      limits=names(comp.colors))+
    scale_y_discrete(
      limits=proc.levs, # necessary for red dot on bottom.
      name="process")+
    scale_x_continuous("time (seconds)")
}
polygon.mine <- function
### Make a Positioning Method that places non-overlapping speech
### polygons at the first or last points.
(top.bottom.left.right,
### Character string indicating what side of the plot to label.
  offset.cm=0.1,
### Offset from the polygon to the most extreme data point.
  padding.cm=0.05,
### Padding inside the polygon.
  custom.colors=NULL
### Positioning method applied just before draw.polygons, can set
### box.color and text.color for custom colors.
){
  if(is.null(custom.colors)){
    custom.colors <- directlabels::gapply.fun({
      rgb.mat <- col2rgb(d[["colour"]])
      d$text.color <- with(data.frame(t(rgb.mat)), {
        gray <- 0.3*red + 0.59*green + 0.11*blue
        ifelse(gray/255 < 0.5, "white", "black")
      })
      d
    })
  }
  opposite.side <- c(
    left="right",
    right="left",
    top="bottom",
    bottom="top")[[top.bottom.left.right]]
  direction <- if(
    top.bottom.left.right %in% c("bottom", "left")
  ) -1 else 1
  min.or.max <- if(
    top.bottom.left.right %in% c("top", "right")
  ) max else min
  if(top.bottom.left.right %in% c("left", "right")){
    min.or.max.xy <- "x"
    qp.target <- "y"
    qp.max <- "top"
    qp.min <- "bottom"
    padding.h.factor <- 2
    padding.w.factor <- 1
    limits.fun <- ylimits
    reduce.method <- "reduce.cex.lr"
  }else{
    min.or.max.xy <- "y"
    qp.target <- "x"
    qp.max <- "right"
    qp.min <- "left"
    padding.h.factor <- 1
    padding.w.factor <- 2
    limits.fun <- directlabels::xlimits
    reduce.method <- "reduce.cex.tb"
  }
  list(
    hjust=0.5, vjust=1,
    function(d,...){
      ## set the end of the speech polygon to the original data point.
      for(xy in c("x", "y")){
        extra.coord <- sprintf(# e.g. left.x
          "%s.%s", opposite.side, xy)
        d[[extra.coord]] <- d[[xy]]
      }
      ## offset positions but do NOT set the speech polygon position
      ## to the min or max.
      d[[min.or.max.xy]] <- d[[min.or.max.xy]] + offset.cm*direction
      d
    },
    "calc.boxes",
    reduce.method,
    function(d, ...){
      d$h <- d$h + padding.cm * padding.h.factor
      d$w <- d$w + padding.cm * padding.w.factor
      d
    },
    "calc.borders",
    function(d,...){
      do.call(rbind, lapply(split(d, d$y), function(d){
        directlabels::apply.method(directlabels::qp.labels(
          qp.target,
          qp.min,
          qp.max,
          directlabels::make.tiebreaker(min.or.max.xy, qp.target),
          limits.fun), d)
      }))
    },
    "calc.borders",
    custom.colors,
    "draw.polygons")
}
viz_workers(loss_5_15)
```

The figure above shows how the computation proceeds over time (X axis)
in the different parallel processes (Y axis). Iterations are shown in
grey rectangles, with iteration numbers at the top of the plot. In
each iteration, There are at most four parallel workers (black) at any
given time. We see that before starting a new iteration, all workers
must complete, and wait for the central process to compute new
candidates. The worker efficiency reported at the top of the plot is
the amount of time taken computing models (black segments), divided by
the max time that could be taken by that number of workers (if black
line segments were from start to end of the plot for each
worker). There are two sources of inefficiency:

* The overhead for starting a new future job is about 0.2 seconds, and
  this overhead could be reduced by using mirai package instead. But
  for long running computations (big data or complex models), this
  overhead is not the important bottleneck.
* Because of the centralized method for computing new candidates, the
  first process that finishes in an iteration must wait for the last
  process, before it starts working again. This can be fixed by using
  rush package instead (as we show in the next section).

To see these patterns more clearly, below we run the central launch
algorithm with a larger number of models:

```{r central-1-100-future_lapply}
results_1_100 <- list()
(results_1_100$future_lapply <- central_launch(1,100))
viz_workers(results_1_100$future_lapply)
```

The figure above shows more iterations and more processes. In some
iterations (9-11), the number of candidates is greater than 14, which
is the number of CPUs on my machine (and the max number of future
workers). In those iterations, we see calculation of either 1 or 2
models in each worker process.

### Centralized launching, no parallelization

A baseline to compare is no parallel computation (everything in one
process), as coded below.

```{r central-1-100-lapply}
results_1_100$lapply <- central_launch(1,100,lapply)
viz_workers(results_1_100$lapply)
```

The result above is almost 100% efficient (because only one CPU is
used instead of 14), but it takes longer overall.

### Centralized launching, mirai

Another comparison is mirai, which offers lower overhead than `future`.

```{r central-1-100-mirai_map}
if(mirai::daemons()$connections==0){
  mirai::daemons(future::availableCores())
}
results_1_100$mirai_map <- central_launch(1,100,function(...){
  mirai::mirai_map(...)[]
})
viz_workers(results_1_100$mirai_map)
```

The figure above shows that there is very little overhead for
launching mirai parallel tasks. In each iteration, we can see that all
14 processes start almost at the same time. This results in a shorter
overall comptutation time. However, we can see still some
overhead/inefficiency that comes from the centralized computation of
target penalties. That is, at each iteration, some processes are
short, and must wait until the longest process in the iteration
finishes, before receiving a new penalty to compute. This observation
motivates the de-centralized parallelized model that should be
possible using rush.

### Comparison of central launching methods

In the code below, we combine the results from the three methods in the previous sections.

```{r compare-times}
loss_1_100_list <- list()
fun_levs <- c("lapply", "future_lapply", "mirai_map")
for(fun in names(results_1_100)){
  loss_fun <- results_1_100[[fun]]$loss
  min.time <- min(loss_fun$start.time)
  loss_1_100_list[[fun]] <- data.table(fun=factor(fun, fun_levs), loss_fun)[, let(
    process_i = as.integer(process),
    start.time=start.time-min.time,
    end.time=end.time-min.time)]
}
(loss_1_100 <- rbindlist(loss_1_100_list))
ggplot()+
  geom_segment(aes(
    start.time, process_i,
    xend=end.time, yend=process_i),
    data=loss_1_100)+
  geom_point(aes(
    start.time, process_i),
    shape=1,
    data=loss_1_100)+
  facet_grid(fun ~ ., scales="free", labeller=label_both)+
  scale_y_continuous(breaks=seq(1, parallel::detectCores()))+
  scale_x_continuous("Time from start of computation (seconds)")
```

In the figure above, we can see that `lapply` and `future_lapply` take
the most time, and `mirai_map` is much faster (about 5 seconds). This
is an example when `mirai_map` is particularly advantageous:

* centralized process does not take much time to assign new tasks.
* computation time of each task is relatively small, so overhead of `future_lapply` launching is relevant.

For longer running computations, for example several minutes or
seconds, there should be smaller differences between `future_lapply`
and `mirai_map`.

## Parallelization without a central launcher

In this section, we explore the speedups that we can get by using a
de-centralized parallel computation, in which each process decides
what penalty to compute (there is no central launcher, so one less bottleneck).

### Declare some functions

We coordinate the workers using a lock file, in the functions below.
The function below creates a new data set directory, and returns the path where we can save model info in an RDS file.

```{r}
new_data_dir <- function(){
  data.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11-60000-580000")
  dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
  data(Mono27ac, package="PeakSegDisk")
  write.table(
    Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  file.path(data.dir, "summary.rds")
}
```

The function below attempts to run one new/helpful penalty, based on the models in the RDS file.
Note that we use the [filelock](https://github.com/r-lib/filelock) R package to
manage a lock file that can guarantee only one process reading/writing this RDS file at a time.
Below we use `lock()` and `unlock()` to protect parts of the code where we need to read from the RDS and then write it again.

* in the first lock/unlock block, we read RDS, decide on a new `pen` value to compute, then add a row with that penalty and `peaks=NA` to the RDS table.
* at this point, if another process looks at the RDS, it will not decide to compute the same `pen` value, because it already exists in the RDS table.
* then we compute the model for the new `pen` value via `PeakSegFPOP_dir()`.
* then in the second lock/unlock block, we read RDS, and save the new `peaks` value to RDS.

```{r}
run_penalty <- function(){
  library(data.table)
  data.dir <- dirname(summary.rds)
  lock.file <- paste0(summary.rds, ".LOCK")
  target.max.peaks <- 100
  target.min.peaks <- 1
  lock.obj <- filelock::lock(lock.file)
  start.time.candidates <- Sys.time()
  if(file.exists(summary.rds)){
    task_dt <- readRDS(summary.rds)
    finite_dt <- task_dt[is.finite(peaks)]
  }else{
    finite_dt <- task_dt <- data.table()
  }
  candidates <- if(nrow(finite_dt) < 2){
    c(0, Inf)
  }else{
    selection.df <- penaltyLearning::modelSelection(finite_dt, "total.loss", "peaks")
    selection.dt <- with(selection.df, data.table(
      total.loss,
      min.lambda,
      penalty,
      max.lambda,
      peaks_after=c(peaks[-1],NA),
      peaks
    ))[
      peaks_after<target.max.peaks & peaks>target.min.peaks & peaks_after+1 < peaks,
      max.lambda
    ]
  }
  new.cand <- setdiff(candidates, task_dt$penalty)
  if(length(new.cand)){
    pen <- new.cand[1]
    new_task_dt <- rbind(
      task_dt,
      data.table(
        penalty=pen,
        total.loss=NA_real_,
        process=factor(Sys.getpid()),
        peaks=NA_integer_,
        start.time.candidates,
        end.time.candidates=Sys.time(),
        start.time.model=NA_real_,
        end.time.model=NA_real_))
    saveRDS(new_task_dt, summary.rds)
  }else{
    pen <- NULL
  }
  filelock::unlock(lock.obj)
  if(is.numeric(pen)){
    start.time <- Sys.time()
    fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, pen)
    end.time <- Sys.time()
    lock.obj <- filelock::lock(lock.file)
    task_dt <- readRDS(summary.rds)
    task_dt[penalty==pen, let(
      total.loss=fit$loss$total.loss,
      peaks=fit$loss$peaks,
      start.time.model=start.time,
      end.time.model=end.time
    )]
    saveRDS(task_dt, summary.rds)
    filelock::unlock(lock.obj)
  }
  pen
}
```

The code below is the main worker loop of each parallel process. The
function includes a while loop that repeatedly does `run_penalty()`.
If that returns numeric (indicating a model was computed for a new
penalty), then the time is saved. If more than 5 seconds has elapsed
since the last model was computed, then we are done.

```{r}
run_penalties <- function(seconds.thresh=5){
  done <- FALSE
  last.time <- Sys.time()
  while(!done){
    pen <- run_penalty()
    if(is.numeric(pen)){
      last.time <- Sys.time()
    }
    if(Sys.time()-last.time > seconds.thresh){
      done <- TRUE
    }
  }
}
```

The function below computes a summary table based on the data in the RDS file.

```{r}
reshape_summary <- function(){
  (summary_dt <- readRDS(summary.rds))
  nc::capture_melt_multiple(
    summary_dt[, let(
      row = .I,
      process_i = as.integer(process)
    )],
    column=".*", "[.]",
    computation="model|candidates"
  )[, let(
    start.seconds=start.time-min(start.time),
    end.seconds=end.time-min(start.time)
  )][]
}
```

The function below plots the summary table.

```{r}
comp.colors <- c(
  model="black",
  candidates="red")
gg_summary <- function(summary_long){
  library(ggplot2)
  ggplot()+
    theme_bw()+
    geom_segment(aes(
      start.seconds, process,
      color=computation,
      xend=end.seconds, yend=process),
      data=summary_long)+
    geom_point(aes(
      start.seconds, process,
      size=computation,
      color=computation),
      shape=1,
      data=summary_long)+
    scale_color_manual(values=comp.colors)+
    scale_size_manual(values=c(
      candidates=3,
      model=2))+
    directlabels::geom_dl(aes(
      start.seconds, process,
      label.group=row,
      label=peaks),
      data=summary_long[computation=="model"],
      method=polygon.mine("bottom", offset.cm=0.2, custom.colors=list(box.color="white")))+
    scale_x_continuous("Time from start of computation (seconds)")
}
```

### Function demo

Below we create a new data directory, then run three penalties.

```{r}
summary.rds <- new_data_dir()
run_penalty()
readRDS(summary.rds)
run_penalty()
readRDS(summary.rds)
run_penalty()
readRDS(summary.rds)
```

The tables above show that each call to `run_penalty()` adds a row with a new penalty value to the RDS file.
The idea is to repeated call this function in many parallel processes.

### De-centralized parallelization, lock file, future

Below we use `future_lapply()` to implement the de-centralized parallel computation.

```{r filelock-future}
future::plan("multisession")
summary.rds <- new_data_dir()
null.list <- future.apply::future_lapply(seq(1, parallel::detectCores()), function(i)run_penalties())
summary_future <- reshape_summary()
gg_summary(summary_future)
```

The plot above shows how the parallel computation proceeds over time using `future_lapply()`.

* The first process starts with 3199 peaks at the bottom.
* There is a delay of about 0.2 seconds before starting computation in each other process, as can be seen via the stair step pattern going up to the right.
* With the previous centralized launcher result, when a process
  finished, it had to wait for the other processes to finish. Then the
  launcher would assign new work, at which case the process could
  resume work. Here we see a different pattern: when a computation
  finishes, a new computation starts almost immediately, with no
  visible periods where processes have to wait for each other. This
  reduction of idle time is the main advantage of the de-centralized
  method of parallelization.
* The computation finishes in about 4 seconds, after each process has
  decided that there is no more work to do.

### De-centralized parallelization, lock file, mirai

Below we use `mirai_map()` to implement the de-centralized parallel computation.

```{r filelock-mirai}
if(mirai::daemons()$connections==0){
  mirai::daemons(parallel::detectCores())
}
summary.rds <- new_data_dir()
list.of.mirai <- mirai::mirai_map(
  1:parallel::detectCores(),
  function(i)run_penalties(),
  run_penalties=run_penalties,
  run_penalty=run_penalty,
  summary.rds=summary.rds
)[]
summary_mirai <- reshape_summary()
gg_summary(summary_mirai)
```

The figure above shows similar patterns as the previous section.

* A less pronounced stair step pattern is evident, because there is a smaller delay between creation of processes.
* There is a somewhat large gap in the bottom left, which is due to the lack of obviously helpful penalties to explore in the beginning of the algorithm.

### Comparing de-centralized parallelization packages

The figure below compares the two R functions that we used for de-centralized parallelization.

```{r filelock-compare}
summary_compare <- rbind(
  data.table(fun="mirai_map", summary_mirai),
  data.table(fun="future_lapply", summary_future))
ggplot()+
  theme_bw()+
  geom_segment(aes(
    start.seconds, process_i,
    color=computation,
    xend=end.seconds, yend=process_i),
    data=summary_compare)+
  geom_point(aes(
    start.seconds, process_i,
    size=computation,
    color=computation),
    shape=1,
    data=summary_compare)+
  scale_color_manual(values=comp.colors)+
  scale_size_manual(values=c(
    candidates=3,
    model=2))+
  scale_y_continuous(breaks=seq(1, parallel::detectCores()))+
  scale_x_continuous("Time from start of computation (seconds)")+
  facet_grid(fun ~ ., labeller=label_both, scales="free", space="free")
```

The figure above shows timings on my laptop: `mirai_map()` is actually a bit slower
than `future_lapply()` in this case, although the difference is not
large (about 1 second).

Figures below are the same analysis on two Alliance Canada Clusters.

![mammouth-projects](/assets/img/2025-05-15-rush-change-point/mammouth-projects.png)

![beluga-scratch](/assets/img/2025-05-15-rush-change-point/beluga-scratch.png)

The results above show that:

* the file lock mechanism works on both file systems (/projects and /scratch).
* mirai is much faster than future on
  mammouth. `parallel::detectCores()` returns 24 on the login node,
  and this is the number of CPUs used by mirai, but not future.
* future is slightly faster than mirai on
  beluga. `paralllel::detectCores()` returns 40 on the login node, but
  both frameworks only use 6 CPUs.

### Compare everything

The code below combines all timings into a single visualization.

```{r all-compare, fig.height=8}
loss_1_100[, let(
  start.seconds=start.time-min(start.time),
  end.seconds=end.time-min(start.time)
)]
common.names <- intersect(names(loss_1_100), names(summary_compare))
all_compare <- rbind(
  data.table(design="de-centralized", summary_compare[, ..common.names]),
  data.table(design="centralized", loss_1_100[, ..common.names])
)[, Fun := factor(fun, fun_levs)][]
ggplot()+
  theme_bw()+
  geom_segment(aes(
    start.seconds, process_i,
    color=design,
    xend=end.seconds, yend=process_i),
    data=all_compare)+
  geom_point(aes(
    start.seconds, process_i, color=design),
    shape=1,
    data=all_compare)+
  facet_grid(Fun + design ~ ., scales="free", labeller=label_both)+
  scale_y_continuous(breaks=seq(1, parallel::detectCores()))+
  scale_x_continuous("Time from start of computation (seconds)")
```

The figure above compares the different parallelization methods we have explored.

* The two methods using `mirai_map` are fast because the delay for launch new processes is very small.
* The two methods using de-centralized parallelization are fast because there is no need to wait for a central launcher to assign work. This difference is especially evident when comparing the two `future_lapply` methods.

## Conclusions

We have discussed how to implement Change-points for a Range Of
ComplexitieS (CROCS), an algorithm that returns a range of peak
models. We used a real genomic data set, and computed models with from
1 to 100 peaks. We observed the differences between several computation methods:

* `lapply` sequential computation is relatively slow.
* `future.apply::future_lapply` results in somewhat faster computation
  using 14 CPUs on my laptop. The computation time of each job was on
  the same order of magnitude as the overhead of launching parallel
  jobs (about 0.2 seconds). For much longer jobs, this overhead is not
  significant, but for these small jobs, it does result in noticeable
  slow-downs.
* `mirai::mirai_map` results in even faster computation, because its
  overhead is much smaller.
* Both functions can be used with de-centralized parallelization,
  which greatly reduces computation time using `future_lapply` in this
  case.

## Session Info

```{r}
sessionInfo()
```

## Not working yet

For future work, it would be interesting to compare with these other frameworks.

### De-centralized candidate computation, rush

TODO, this is not working, but I asked for help <https://github.com/mlr-org/rush/issues/44>

```{r}
data.rush <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11-60000-580000")
dir.create(data.rush, recursive=TRUE, showWarnings=FALSE)
data(Mono27ac,package="PeakSegDisk")
write.table(
  Mono27ac$coverage, file.path(data.rush, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## rush = rush::RushWorker$new(network_id = "test", config=config, remote=FALSE)
## key = rush$push_running_tasks(list(list(penalty=0)))
## rush$fetch_tasks()
## rush$push_results(key, list(list(peaks=5L, total.loss=5.3)))
## rush$fetch_tasks()
run_penalty <- function(rush){
  target.max.peaks <- 15
  target.min.peaks <- 5
  get_tasks <- function(){
    task_dt <- rush$fetch_tasks()
    if(is.null(task_dt$penalty)){
      task_dt$penalty <- NA_real_
    }
    if(is.null(task_dt$peaks)){
      task_dt$peaks <- NA_integer_
    }
    task_dt
  }
  task_dt <- get_tasks()
  start.time.cand <- Sys.time()
  first_pen_cand <- c(0, Inf)
  done <- first_pen_cand %in% task_dt$penalty
  pen <- if(any(!done)){
    first_pen_cand[!done][1]
  }else{
    print(task_dt)
    while(nrow(task_dt[!is.na(peaks)])<2){
      task_dt <- get_tasks()
    }
    print(task_dt)
    selection.df <- penaltyLearning::modelSelection(task_dt, "total.loss", "peaks")
    selection.dt <- with(selection.df, data.table(
      total.loss,
      min.lambda,
      penalty,
      max.lambda,
      peaks_after=c(peaks[-1],NA),
      peaks
    ))[
      peaks_after<target.max.peaks & peaks>target.min.peaks &
        peaks_after+1 < peaks & !max.lambda %in% task_dt$penalty
    ]
    if(nrow(selection.dt)){
      selection.dt[1, max.lambda]
    }
  }
  if(is.numeric(pen)){
    key = rush$push_running_tasks(xss=list(list(
      penalty=pen,
      start.time.cand=start.time.cand,
      end.time.cand=Sys.time())))
    start.time.model <- Sys.time()
    fit <- PeakSegDisk::PeakSegFPOP_dir(data.rush, pen)
    rush$push_results(key, yss=list(list(
      peaks=fit$loss$peaks,
      total.loss=fit$loss$total.loss,
      start.time.model=start.time.model,
      end.time.model=Sys.time())))
    pen
  }
}
wl_penalty <- function(rush){
  done <- FALSE
  while(!done){
    pen <- run_penalty(rush)
    done <- is.null(pen)
  }
}

if(FALSE){
  redux::hiredis()$pipeline("FLUSHDB")
  config = redux::redis_config()
  rush = rush::Rush$new(network_id = "test", config=config)
  wl_penalty(rush)
  devtools::install_github("mlr-org/rush@mirai")
}


wl_random_search = function(rush) {
  # stop optimization after 100 tasks
  while(rush$n_finished_tasks < 100) {
    # draw new task
    xs = list(x1 = runif(1, -5, 10), x2 = runif(1, 0, 15))
    # mark task as running
    key = rush$push_running_tasks(xss = list(xs))
    # evaluate task
    ys = list(y = branin(xs$x1, xs$x2))
    # push result
    rush$push_results(key, yss = list(ys))
  }
}
branin = function(x1, x2) {
  (x2 - 5.1 / (4 * pi^2) * x1^2 + 5 / pi * x1 - 6)^2 + 10 * (1 - 1 / (8 * pi)) * cos(x1) + 10
}


library(data.table)

redux::hiredis()$pipeline("FLUSHDB")
config = redux::redis_config()
rush = rush::Rush$new(network_id = "test", config=config)
rush$fetch_tasks()
rush$start_local_workers(
  worker_loop = wl_random_search,
  n_workers = 4,
  globals = "branin")
rush
rush$fetch_tasks()

redux::hiredis()$pipeline("FLUSHDB")
config = redux::redis_config()
rush = rush::Rush$new(network_id = "test", config=config)
rush$fetch_tasks()
rush$start_local_workers(
  worker_loop = wl_penalty,
  n_workers = 4,
  globals = c("run_penalty","data.rush"))
rush
task_dt <- rush$fetch_tasks()
task_dt[order(peaks)]
```

### redis list

TODO code based on redux and
[transactions](https://redis.io/docs/latest/develop/interact/transactions/)
or [locks](https://redis.io/docs/latest/develop/use/patterns/distributed-locks/)?

```{r}
r <- redux::hiredis()
r$FLUSHDB()
r$LPUSH()
```


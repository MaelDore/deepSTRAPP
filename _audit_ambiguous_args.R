###########################################################################################
##   Audit of ambiguous / partially-matched argument names                               ##
##   Author: Mael Dore                                                                   ##
##   For every function holding a "..." argument, reports the formals declared BEFORE it  ##
##   that a graphical parameter passed through "..." would partially match, and therefore ##
##   silently capture. Run from the root of the package directory.                        ##
###########################################################################################

library(ape)
## ---- the set of names a user may plausibly pass through '...' -------------------
pdf(NULL); par_names <- names(par(no.readonly = TRUE)); dev.off()
plot_args <- c("main","sub","xlab","ylab","xlim","ylim","type","add","axes","ann","log","asp",
               "border","angle","density","labels","legend","offset","direction","edge.width",
               "show.tip.label","tip.color","label.offset","align.tip.label","node.pos","root.edge",
               "fsize","ftype","outline","sig","hold","leg.txt","plot","x.lim","y.lim","lend","srt")
ape_args <- names(formals(ape::plot.phylo))
## BAMMtools formals, parsed from source (not installed here)
bt <- character(0)
for (f in list.files("/tmp/bt/BAMMtools/R", pattern="\\.R$", full.names=TRUE)) {
  e <- tryCatch(parse(f), error=function(e) NULL)
  for (x in e) if (is.call(x) && as.character(x[[1]])[1] %in% c("<-","=") &&
                   is.call(x[[3]]) && as.character(x[[3]][[1]])[1]=="function" &&
                   as.character(x[[2]]) %in% c("plot.bammdata","addBAMMshifts","addBAMMlegend","dtRates"))
    bt <- c(bt, names(as.list(x[[3]][[2]])))
}
risk_names <- setdiff(unique(c(par_names, plot_args, ape_args, bt)), c("...",""))

## ---- formals of every deepSTRAPP function ---------------------------------------
fun <- list(); where <- list()
for (f in list.files("R", pattern="\\.R$", full.names=TRUE)) {
  e <- tryCatch(parse(f), error=function(e) NULL)
  for (x in e) if (is.call(x) && as.character(x[[1]])[1] %in% c("<-","=") &&
                   is.call(x[[3]]) && as.character(x[[3]][[1]])[1]=="function") {
    nm <- as.character(x[[2]]); fun[[nm]] <- names(as.list(x[[3]][[2]])); where[[nm]] <- basename(f) }
}
exported <- sub("^export\\((.*)\\)$","\\1", grep("^export\\(", readLines("NAMESPACE"), value=TRUE))

cat(sprintf("%d risk names; %d functions; checking those with '...'\n\n", length(risk_names), length(fun)))
res <- NULL
for (nm in names(fun)) {
  fm <- fun[[nm]]
  if (!("..." %in% fm)) next
  before <- fm[seq_len(match("...", fm) - 1L)]
  for (p in risk_names) {
    if (p %in% before) next                       # exact formal: user cannot reach '...' with it anyway
    hit <- before[startsWith(before, p) & before != p]
    if (length(hit) == 1L)
      res <- rbind(res, data.frame(fn=nm, file=where[[nm]], passed=p, effect="SILENT CAPTURE",
                                   detail=hit, exported=nm %in% exported, stringsAsFactors=FALSE))
    else if (length(hit) > 1L)
      res <- rbind(res, data.frame(fn=nm, file=where[[nm]], passed=p, effect="hard error",
                                   detail=paste(hit, collapse=" / "), exported=nm %in% exported, stringsAsFactors=FALSE))
  }
}
res <- res[order(res$effect, res$fn, res$passed), ]
cat("======== SILENT CAPTURE (wrong behaviour, no error) ========\n")
s <- res[res$effect=="SILENT CAPTURE",]
for (i in seq_len(nrow(s))) cat(sprintf("  %-52s %-14s -> %-26s %s\n", s$fn[i], s$passed[i], s$detail[i], ifelse(s$exported[i],"[exported]","")))
cat("\n======== hard error if passed ========\n")
h <- res[res$effect=="hard error",]
for (i in seq_len(nrow(h))) cat(sprintf("  %-52s %-14s -> %s\n", h$fn[i], h$passed[i], h$detail[i]))
cat(sprintf("\n%d silent captures, %d hard errors\n", nrow(s), nrow(h)))

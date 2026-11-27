# Exports the random forest of a CHORD model (.rds, as chord_predict.R reads it) into the text resource
# that chord-rs embeds (chord-rs/resources/chord_forest.tsv), so that the port predicts without R.
#
# Usage: Rscript export_forest.R <CHORD.rds> <output.tsv>
# Needs the randomForest package only to print its version (readRDS reads the object without it).
#
# Format (tab-separated, '#' lines are comments):
#   type        classification
#   classes     <class 1> ... <class k>          (object$classes, the order of the vote columns)
#   variables   <var 1> ... <var m>              (rownames(object$importance): the order predict uses)
#   ncat        <1 per variable>                 (1 = numeric; other values are not supported)
#   ntree       <n>
#   tree        <tree number> <node count>       (object$forest$ndbigtree[tree])
#   then one line per node of that tree, in node order:
#   <nodestatus> <bestvar> <left daughter> <right daughter> <nodepred> <xbestsplit>
# xbestsplit is written as a C99 hexadecimal float (sprintf "%a", e.g. 0x1.881723c1f2245p-4), which
# holds the double exactly; the script checks that R reads each one back to the same double. (Decimal
# with 17 digits would hold it too, but R's own string-to-double conversion is not correctly rounded, so
# the check could not be made in R.)

args <- commandArgs(TRUE)
if (length(args) != 2) stop("usage: Rscript export_forest.R <CHORD.rds> <output.tsv>")

m <- readRDS(args[1])
if (!inherits(m, "randomForest")) stop("not a randomForest object")
if (m$type != "classification") stop("only classification forests are supported")
f <- m$forest
if (any(f$ncat != 1)) stop("categorical predictors are not supported")

vars <- if (is.null(dim(m$importance))) names(m$importance) else rownames(m$importance)
if (!identical(names(f$ncat), vars)) stop("ncat names differ from the importance rownames")

fmt_split <- function(x) sprintf("%a", x)

lines <- c(
  sprintf("# CHORD random forest exported by export_forest.R from %s", basename(args[1])),
  sprintf("# R %s, randomForest %s", paste(R.version$major, R.version$minor, sep = "."),
          tryCatch(as.character(packageVersion("randomForest")), error = function(e) "not installed")),
  paste("type", m$type, sep = "\t"),
  paste(c("classes", m$classes), collapse = "\t"),
  paste(c("variables", vars), collapse = "\t"),
  paste(c("ncat", as.integer(f$ncat)), collapse = "\t"),
  paste("ntree", as.integer(f$ntree), sep = "\t")
)

out <- file(args[2], "w")
writeLines(lines, out)
max_diff <- 0
for (t in seq_len(f$ntree)) {
  n <- f$ndbigtree[t]
  split <- f$xbestsplit[1:n, t]
  split_text <- fmt_split(split)
  if (!identical(as.numeric(split_text), as.numeric(split))) stop(sprintf("split values of tree %d do not round-trip", t))
  writeLines(sprintf("tree\t%d\t%d", t, n), out)
  writeLines(paste(f$nodestatus[1:n, t], f$bestvar[1:n, t], f$treemap[1:n, 1, t], f$treemap[1:n, 2, t],
                   f$nodepred[1:n, t], split_text, sep = "\t"), out)
}
close(out)
cat(sprintf("wrote %s: %d trees, %d nodes\n", args[2], f$ntree, sum(f$ndbigtree)))

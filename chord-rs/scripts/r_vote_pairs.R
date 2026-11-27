# Every pair of vote fractions the model can give (k1/500, k2/500 with k1 + k2 <= 500) and their sum,
# as chord_predict.R computes them (rowSums) and write.table writes them. Usage: Rscript r_vote_pairs.R <out.tsv>
args <- commandArgs(TRUE)
g <- expand.grid(k1 = 0:500, k2 = 0:500)
g <- g[g$k1 + g$k2 <= 500, ]
probs <- data.frame(BRCA1 = g$k1 / 500, BRCA2 = g$k2 / 500)
out <- data.frame(k1 = g$k1, k2 = g$k2, p_BRCA1 = probs$BRCA1, p_BRCA2 = probs$BRCA2, p_hrd = rowSums(probs))
write.table(out, args[1], sep = "\t", quote = FALSE, row.names = FALSE)

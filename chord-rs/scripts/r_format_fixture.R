# The fixture of chord-rs's write.table number test: each value as an exact hex float, with the text
# write.table writes for it. Usage: Rscript r_format_fixture.R <out.tsv>
args <- commandArgs(TRUE)
set.seed(42)
special <- c(0, 1, -1, 0.5, -0.25, 0.004, 0.002, 0.006, 0.586, 1/3, 2/3, 1e5, 1e-5, 1e-4, 0.00012, 123456, 1234567,
             1e15, 1e16, 1e-99, 1e-100, 1e99, 1e100, 1e300, 1e-300, 9.9999999999999995, 99999.5, 0.1 + 0.2,
             123456789.123, -123456789.123, 2^-1074, .Machine$double.xmax, 5e-324 * 3, 999999999999999, 9999999999999999,
             0.000099999, 1 - 2^-53, 1 + 2^-52, 100000, 1e-15, 314159.26535)
votes <- (0:500) / 500
sums <- (0:500) / 500 + rev((0:500) / 500) * 0.5
mags <- 10^runif(400, -320, 308) * sample(c(-1, 1), 400, TRUE)
digits <- round(runif(300) * 10^sample(0:15, 300, TRUE)) / 10^sample(0:15, 300, TRUE)
x <- c(special, votes, sums, mags, digits)
x <- x[is.finite(x)]
df <- data.frame(hex = sprintf("%a", x), value = x, stringsAsFactors = FALSE)
write.table(df, args[1], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bw_pair_cnts <- function(x, nb, binned)
{
  if(binned) {
    r <- range(x)
    d <- diff(r) * 1.01/nb
    xx <- trunc(abs(x)/d) *sign(x)
    xx <- xx - min(xx) + 1
    xxx <- tabulate(xx, nb)
    list(d, bw_den_binned(xxx))
  } else bw_den(nb, x)
}
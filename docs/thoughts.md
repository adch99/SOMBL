# Thoughts and Questions to Ask

1. Should I be checking first for the optimal number of
    eigenvalues or should I directly use the total number
    of states? It may not be worth optimizing for though.
    The lapack docs tell to do a double call to zheev
    to get optimal `LWORK` and then set the `WORK` param.

2. Should I average over the localization lengths
    $ξ(E_i)$ or over the lyapunov exponents $γ(E_i)$?

3. What library should I use for sparse matrix 
    calculations?
# Thoughts and Questions to Ask

> Turns out I am an idiot and didn't realize that the
    eigenvalue method worked only for 1D open boundary
    systems. Now I need to find other methods to deal
    with this. Most will involve getting eigenfunctions
    which is obviously not great computationally. 

- Should I average over the localization lengths
    $ξ(E_i)$ or over the lyapunov exponents $γ(E_i)$?

- What library should I use for sparse matrix 
    calculations?

- I keep getting this feeling that the Lyapunov
    exponent calculation is getting messed up
    somewhere. It keeps giving negative localization
    lengths. Answer: It was getting messed up because
    of passing the wrong parameter.

## Done

- Should I be checking first for the optimal number of
    eigenvalues or should I directly use the total number
    of states? It may not be worth optimizing for though.
    The lapack docs tell to do a double call to zheev
    to get optimal `LWORK` and then set the `WORK` param.
    Answer: The LAPACKE functions set the workspace
    variables on their own unless we want to do it
    explicitly. So this is not needed.


## Parts of the Program

1. Initializing parameters
2. get_neighbour_lists - Tested - OK
3. hamiltonian_nospin - Tested - OK
4. utils_get_eigvalsh - Tested - OK
5. analysis - Not Tested
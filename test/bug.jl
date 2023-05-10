function test(ztol = Base.rtoldefault(Float64))
    X = [0.5459627556242905, 1.7288950489507429, 0.7167681447476535]
    Y = [-54.06002080721971, 173.77393714162503, 71.48154370522498]
    n = 3
    @polyvar W[1:n] α β
    I = @set(
        sum([W[i] * X[i] * (β * X[i] + α - Y[i]) for i in 1:n]) == 0,
        library = Buchberger(ztol)
    )
    for i in 1:n
        addequality!(I, W[i] - W[i]^2)
    end
    SemialgebraicSets.computegröbnerbasis!(I.I)
    return I
end

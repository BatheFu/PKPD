using CairoMakie

using AlgebraOfGraphics

logistic_fun(x) = 1 / (1 + exp(-x));

data((; x = -10:0.1:10, y = logistic_fun.(-10:0.1:10))) *
mapping(:x, :y => "logistic(x)") *
visual(Lines; color = :steelblue, linewidth = 3) |> draw


using Distributions

x = [0, 1]
probs = [0.5, 0.25, 0.75]
dists = [Bernoulli(p) for p in probs]
pdfs = mapreduce(d -> pdf.(d, x), vcat, dists)
plt =
    data((; probs = repeat(probs; inner = 2), x = repeat(x; outer = 3), pdfs)) *
    mapping(
        :x => "outcome",
        :pdfs => "PMF";
        color = :probs => nonnumeric => "parameter p",
        col = :probs => nonnumeric,
    ) *
    visual(BarPlot)
draw(plt; axis = (; xticks = [0, 1]))

using  Statistics, DataFrames, Chain, BenchmarkTools

df = DataFrame(
    x = rand(["A", "B", "C", "D"], 10_000),
    y = rand(10_000),
    z = rand(10_000),
)

benchResults = @benchmark @chain df begin
    groupby(:x)
    combine(:y => median, :z => mean)
end
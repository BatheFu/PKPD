using Pumas
using PumasUtilities
using NCA
using NCAUtilities
using GLM: lm, @formula
using Random
using CSV
using DataFramesMeta
using CairoMakie
using PharmaDatasets

pkpain_df = dataset("pk_painrelief")
first(pkpain_df, 5)

pkpain_noplb_df = @rsubset pkpain_df :Dose != "Placebo";
first(pkpain_noplb_df, 5)

@rtransform! pkpain_noplb_df begin
    :route = "ev"
    :Dose = parse(Int, chop(:Dose; tail = 3))
end

@rtransform! pkpain_noplb_df :amt = :Time == 0 ? :Dose : missing

pkpain_nca = read_nca(
    pkpain_noplb_df;
    id = :Subject,
    time = :Time,
    amt = :amt,
    observations = :Conc,
    group = [:Dose],
    route = :route,
)

f = observations_vs_time(
    pkpain_nca;
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "CTMNoPain Concentration (ng/mL)"),
    facet = (; combinelabels = true),
)
f[1]

summary_observations_vs_time(
    pkpain_nca,
    figure = (; fontsize = 22, resolution = (800, 1000)),
    color = "black",
    linewidth = 3,
    axis = (; xlabel = "Time (hr)", ylabel = "CTMX Concentration (μg/mL)"),
    facet = (; combinelabels = true, linkaxes = true),
)

pk_nca = run_nca(pkpain_nca; sigdigits = 3)

f = subject_fits(
    pk_nca,
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "CTMX Concentration (μg/mL)"),
    facet = (; combinelabels = true, linkaxes = true),
)
f[1]


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

strata = [:Dose]
params = [:cmax, :aucinf_obs]

pk_nca = run_nca(pkpain_nca; sigdigits = 3)

output = summarize(pk_nca; stratify_by = strata, parameters = params)

parameters_vs_group(
    pk_nca,
    parameter = :cmax,
    axis = (; xlabel = "Dose (mg)", ylabel = "Cₘₐₓ (ng/mL)"),
    figure = (; fontsize = 18),
)

@rtransform! pkpain_noplb_df begin
    :evid = :Time == 0 ? 1 : 0
    :cmt = :Time == 0 ? 1 : 2
    :cmt2 = 1
end

@rtransform! pkpain_noplb_df :Conc = :evid == 1 ? missing : :Conc

pkpain_noplb = read_pumas(
    pkpain_noplb_df;
    id = :Subject,
    time = :Time,
    amt = :amt,
    observations = [:Conc],
    covariates = [:Dose],
    evid = :evid,
    cmt = :cmt,
)

pk_1cmp = @model begin

    @metadata begin
        desc = "One Compartment Model"
        timeu = u"hr"
    end

    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init = 3.2)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 16.4)
        """
        Absorption rate constant (h-1)
        """
        tvka ∈ RealDomain(; lower = 0, init = 3.8)
        """
          - ΩCL
          - ΩVc
          - ΩKa
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001, init = 0.2)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates begin
        """
        Dose (mg)
        """
        Dose
    end

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        Ka = tvka * exp(η[3])
    end

    @dynamics Depots1Central1

    @derived begin
        cp := @. Central / Vc
        """
        CTMx Concentration (ng/mL)
        """
        Conc ~ @. Normal(cp, abs(cp) * σ_p)
    end

end

# zero out the random effects
etas = zero_randeffs(pk_1cmp, init_params(pk_1cmp), pkpain_noplb)

simpk_iparams = simobs(pk_1cmp, pkpain_noplb, init_params(pk_1cmp), etas)

sim_plot(
    pk_1cmp,
    simpk_iparams;
    observations = [:Conc],
    figure = (; fontsize = 18),
    axis = (;
        xlabel = "Time (hr)",
        ylabel = "Observed/Predicted \n CTMx Concentration (ng/mL)",
    ),
)
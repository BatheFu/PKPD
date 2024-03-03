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

pkpain_noplb_df = @rsubset pkpain_df :Dose != "Placebo";

@rtransform! pkpain_noplb_df begin
  :route = "ev"
  :Dose = parse(Int, chop(:Dose; tail = 3))
end;

@rtransform! pkpain_noplb_df :amt = :Time == 0 ? :Dose : missing;

pkpain_nca = read_nca(
  pkpain_noplb_df;
  id = :Subject,
  time = :Time,
  amt = :amt,
  observations = :Conc,
  group = [:Dose],
  route = :route,
)

let
  f = observations_vs_time(
    pkpain_nca;
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "CTMNoPain Concentration (ng/mL)"),
    facet = (; combinelabels = true),
  )
  f[1]
end

summary_observations_vs_time(
  pkpain_nca,
  figure = (; fontsize = 22, resolution = (800, 1000)),
  color = "black",
  linewidth = 3,
  axis = (; xlabel = "Time (hr)", ylabel = "CTMX Concentration (μg/mL)"),
  facet = (; combinelabels = true, linkaxes = true),
)

pk_nca = run_nca(pkpain_nca; sigdigits = 3)

let
  f = subject_fits(
    pk_nca,
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "CTMX Concentration (μg/mL)"),
    facet = (; combinelabels = true, linkaxes = true),
  )
  f[1]
end

strata = [:Dose]

parms = [:cmax, :aucinf_obs]

output = summarize(pk_nca; stratify_by = strata, parameters = parms)

parameters_vs_group(
  pk_nca,
  parameter = :cmax,
  axis = (; xlabel = "Dose (mg)", ylabel = "Cₘₐₓ (ng/mL)"),
  figure = (; fontsize = 18),
)

dp = NCA.DoseLinearityPowerModel(pk_nca, :cmax; level = 0.9)

power_model(dp)

dose_vs_dose_normalized(pk_nca, :cmax)

dose_vs_dose_normalized(pk_nca, :aucinf_obs)

@rtransform! pkpain_noplb_df begin
  :evid = :Time == 0 ? 1 : 0
  :cmt = :Time == 0 ? 1 : 2
  :cmt2 = 1 # for zero order absorption
end;

@rtransform! pkpain_noplb_df :Conc = :evid == 1 ? missing : :Conc;

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
    "Clearance (L/hr)"
    tvcl ∈ RealDomain(; lower = 0, init = 3.2)
    "Volume (L)"
    tvv ∈ RealDomain(; lower = 0, init = 16.4)
    "Absorption rate constant (h-1)"
    tvka ∈ RealDomain(; lower = 0, init = 3.8)
    """
    - ΩCL
    - ΩVc
    - ΩKa
    """
    Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04])
    "Proportional RUV"
    σ_p ∈ RealDomain(; lower = 0.0001, init = 0.2)
  end
  @random begin
    η ~ MvNormal(Ω)
  end
  @covariates begin
    "Dose (mg)"
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
  observatins = [:Conc],
  figure = (; fontsize = 18),
  axis = (;
    xlabel = "Time (hr)",
    ylabel = "Observed/Predicted \n CTMx Concentration (ng/mL)",
  ),
)

pkparam = (; init_params(pk_1cmp)..., tvka = 2, tvv = 10)

simpk_changedpars = simobs(pk_1cmp, pkpain_noplb, pkparam, etas)

sim_plot(
  pk_1cmp,
  simpk_changedpars;
  observations = [:Conc],
  figure = (; fontsize = 18),
  axis = (
    xlabel = "Time (hr)",
    ylabel = "Observed/Predicted \n CTMx Concentration (ng/mL)",
  ),
)

pkfit_np = fit(pk_1cmp, pkpain_noplb, init_params(pk_1cmp), NaivePooled(); omegas = (:Ω,))

coefficients_table(pkfit_np)

let
  lls = []
  for subj in pkpain_noplb
    push!(lls, loglikelihood(pk_1cmp, subj, pkparam, FOCE()))
  end
  # the plot below is using native CairoMakie `hist`
  hist(lls; bins = 10, normalization = :none, color = (:black, 0.5), x_gap = 0)
end

influential_subjects = findinfluential(pk_1cmp, pkpain_noplb, pkparam, FOCE())

pkfit_1cmp = fit(pk_1cmp, pkpain_noplb, pkparam, FOCE(); constantcoef = (; tvka = 2))

infer(pkfit_1cmp)

pk_2cmp = @model begin
  @param begin
    "Clearance (L/hr)"
    tvcl ∈ RealDomain(; lower = 0, init = 3.2)
    "Central Volume (L)"
    tvv ∈ RealDomain(; lower = 0, init = 16.4)
    "Peripheral Volume (L)"
    tvvp ∈ RealDomain(; lower = 0, init = 10)
    "Distributional Clearance (L/hr)"
    tvq ∈ RealDomain(; lower = 0, init = 2)
    "Absorption rate constant (h-1)"
    tvka ∈ RealDomain(; lower = 0, init = 1.3)
    """
    - ΩCL
    - ΩVc
    - ΩKa
    - ΩVp
    - ΩQ
    """
    Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04])
    "Proportional RUV"
    σ_p ∈ RealDomain(; lower = 0.0001, init = 0.2)
  end
  @random begin
    η ~ MvNormal(Ω)
  end
  @covariates begin
    "Dose (mg)"
    Dose
  end
  @pre begin
    CL = tvcl * exp(η[1])
    Vc = tvv * exp(η[2])
    Ka = tvka * exp(η[3])
    Vp = tvvp * exp(η[4])
    Q = tvq * exp(η[5])
  end
  @dynamics Depots1Central1Periph1
  @derived begin
    cp := @. Central / Vc
    """
    CTMx Concentration (ng/mL)
    """
    Conc ~ @. Normal(cp, cp * σ_p)
  end
end

pkfit_2cmp =
  fit(pk_2cmp, pkpain_noplb, init_params(pk_2cmp), FOCE(); constantcoef = (; tvka = 2))

compare_estimates(; pkfit_1cmp, pkfit_2cmp)

lrtest(pkfit_1cmp, pkfit_2cmp)

@chain metrics_table(pkfit_2cmp) begin
  leftjoin(metrics_table(pkfit_1cmp); on = :Metric, makeunique = true)
  rename!(:Value => :pk2cmp, :Value_1 => :pk1cmp)
end

res_inspect_1cmp = inspect(pkfit_1cmp)

res_inspect_2cmp = inspect(pkfit_2cmp)

gof_1cmp = goodness_of_fit(res_inspect_1cmp; figure = (; fontsize = 12))

gof_2cmp = goodness_of_fit(res_inspect_2cmp; figure = (; fontsize = 12))

fig_subject_fits = subject_fits(
  res_inspect_2cmp;
  separate = true,
  paginate = true,
  facet = (; combinelabels = true),
  figure = (; fontsize = 18),
  axis = (; xlabel = "Time (hr)", ylabel = "CTMx Concentration (ng/mL)"),
)

fig_subject_fits[1]

empirical_bayes_dist(res_inspect_2cmp; zeroline_color = :red)

empirical_bayes_vs_covariates(
  res_inspect_2cmp;
  categorical = [:Dose],
  figure = (; resolution = (600, 800)),
)

pkfit_2cmp_unfix_ka = fit(pk_2cmp, pkpain_noplb, init_params(pk_2cmp), FOCE())

compare_estimates(; pkfit_2cmp, pkfit_2cmp_unfix_ka)

res_inspect_2cmp_unfix_ka = inspect(pkfit_2cmp_unfix_ka)

goodness_of_fit(res_inspect_2cmp_unfix_ka; figure = (; fontsize = 12))

empirical_bayes_vs_covariates(
  res_inspect_2cmp_unfix_ka;
  categorical = [:Dose],
  ebes = [:η₃],
  figure = (; resolution = (600, 800)),
)

fig_subject_fits2 = subject_fits(
  res_inspect_2cmp_unfix_ka;
  separate = true,
  paginate = true,
  facet = (; combinelabels = true, linkyaxes = false),
  figure = (; fontsize = 18),
  axis = (; xlabel = "Time (hr)", ylabel = "CTMx Concentration (ng/mL)"),
)

fig_subject_fits2[6]

pk_vpc = vpc(
  pkfit_2cmp_unfix_ka,
  200;
  observations = [:Conc],
  stratify_by = [:Dose],
  ensemblealg = EnsembleThreads(), # multi-threading
)

vpc_plot(
  pk_2cmp,
  pk_vpc;
  rows = 1,
  columns = 3,
  figure = (; resolution = (1400, 1000), fontsize = 22),
  axis = (;
    xlabel = "Time (hr)",
    ylabel = "Observed/Predicted\n CTMx Concentration (ng/mL)",
  ),
  facet = (; combinelabels = true),
)
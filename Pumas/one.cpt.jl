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

path = "C:\\Users\\m1825\\Documents\\data\\"
dat = CSV.read(path*"501.csv", DataFrame)

first(dat, 5)

@rtransform! dat begin
    :EVID = :AMT == 0 ? 0 : 1
    :CMT = :AMT > 0 ? 1 : missing
end

@select! dat $(Not(:C))

@rtransform! dat begin
    :AMT = :AMT == 0 ? missing : :AMT
    :DV = :DV == 0 ? missing : :DV
end

first(dat, 3)

pk501 = read_pumas(
    dat;
    id = :ID,
    time = :TIME,
    amt = :AMT,
    observations = [:DV],
    covariates = [:WT, :AGE, :SEX],
    cmt = :CMT,
    evid = :EVID
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
        tvcl ∈ RealDomain(; lower = 0, init = 5)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 10)
        """
          - ΩCL
          - ΩVc
        """
        Ω ∈ PDiagDomain(init = [0.16, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001, init = 0.1)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates begin

    end

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
    end

    @dynamics begin
        Central' = -(CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        """
        CTMx Concentration (ng/mL)
        """
        Conc ~ @. Normal(cp, abs(cp) * σ_p)
    end

end

pkparam = (; init_params(pk_1cmp)..., tvv = 10)

pkfit_1cmp = fit(pk_1cmp, pk501, pkparam, FOCE();)








pkpain_df = dataset("pk_painrelief")
first(pkpain_df, 5)
pkpain_noplb_df = @rsubset pkpain_df :Dose != "Placebo";
first(pkpain_noplb_df, 5)
@rtransform! pkpain_noplb_df begin
    :route = "ev"
    :Dose = parse(Int, chop(:Dose; tail = 3))
end
@rtransform! pkpain_noplb_df :amt = :Time == 0 ? :Dose : missing
@rtransform! pkpain_noplb_df begin
    :evid = :Time == 0 ? 1 : 0
    :cmt = :Time == 0 ? 1 : 2
    :cmt2 = 1 # for zero order absorption
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

pkparam = (; init_params(pk_1cmp)..., tvka = 2, tvv = 10)

pkfit_1cmp = fit(pk_1cmp, pkpain_noplb, pkparam, FOCE(); constantcoef = (; tvka = 2))

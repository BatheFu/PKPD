using Pumas
using PharmaDatasets

df = dataset("po_sad_1")
first(df,5)

describe(df)

#parse df into a Population with read_pumas:
population = read_pumas(df;
    observations = [:dv], 
    covariates = [:wt, :isPM, :isfed],
    route = :route)

#a 2-compartment oral absorption base model with no covariate effects:
base_model = @model begin
    @metadata begin
        desc = "base model"
        timeu = u"hr"
    end
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0)
        """
        Central Volume (L)
        """
        tvvc ∈ RealDomain(; lower = 0)
        """
        Peripheral Volume (L)
        """
        tvvp ∈ RealDomain(; lower = 0)
        """
        Distributional Clearance (L/hr)
        """
        tvq ∈ RealDomain(; lower = 0)
        """
        Absorption rate constant (1/h)
        """
        tvka ∈ RealDomain(; lower = 0)
        """
          - ΩCL
          - ΩVc
          - ΩKa
          - ΩVp
          - ΩQ
        """
        Ω ∈ PDiagDomain(5)
        """
        Proportional RUV (SD scale)
        """
        σₚ ∈ RealDomain(; lower = 0)
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
        Ka = tvka * exp(η[3])
        Q = tvq * exp(η[4])
        Vp = tvvp * exp(η[5])
    end

    @dynamics Depots1Central1Periph1

    @derived begin
        cp := @. 1000 * (Central / Vc)
        """
        Drug Concentration (ng/mL)
        """
        dv ~ @. Normal(cp, cp * σₚ)
    end

end

iparams = (;
    tvka = 0.4,
    tvcl = 4.0,
    tvvc = 70.0,
    tvq = 4.0,
    tvvp = 50.0,
    Ω = Diagonal(fill(0.04, 5)),
    σₚ = 0.1,
)

base_fit = fit(base_model, population, iparams, FOCE())


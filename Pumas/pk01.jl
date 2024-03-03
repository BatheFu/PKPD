using Dates


using Pumas
using PumasUtilities
using Random
using CairoMakie


pk_01        = @model begin
  @metadata begin
    desc = "One Compartment Model"
    timeu = u"minute"
  end
  @param begin
    "Clearance (L/hr)"
    tvcl     ∈ RealDomain(lower=0)
    "Volume (L)"
    tvvc     ∈ RealDomain(lower=0)
    Ω        ∈ PDiagDomain(2)
    "Additive RUV"
    σ        ∈ RealDomain(lower=0)
  end

  @random begin
    η        ~ MvNormal(Ω)
   end

  @pre begin
    Cl       = tvcl * exp(η[1])
    Vc       = tvvc * exp(η[2])
   end

  @dynamics begin
    Central' = -(Cl/Vc)*Central
   end

  @derived begin
    """
    PK01 Concentrations (ug/L)
    """
    cp       = @. 1000*(Central/Vc)
    """
    PK01 Concentrations (ug/L)
    """
    dv       ~ @. Normal(cp, σ)
   end
end


param = [
          (tvcl = 0.10,
          tvvc = 9.98,
          Ω    = Diagonal([0.00,0.00]),
          σ    = 20.80), 

          (tvcl = 0.20,
          tvvc = 9.82,
          Ω    = Diagonal([0.00,0.00]),
          σ    = 27.46),

          (tvcl = 0.20,
          tvvc = 10.22,
          Ω    = Diagonal([0.00,0.00]),
          σ    = 8.78),

          (tvcl = 0.20,
          tvvc = 19.95,
          Ω    = Diagonal([0.00,0.00]),
          σ    = 8.50)
        ]


ev1   = DosageRegimen(10,time=0,cmt=1)

ids   = ["1: CL = 0.10, V = 9.98", "2: CL = 0.20, V = 9.82",
         "3: CL = 0.20, V = 10.22", "4: CL = 0.20, V = 19.95"]

pop = map(i -> Subject(id = ids[i], events = ev1, observations = (cp = nothing,)), 1:length(ids))


Random.seed!(123)

sim = map(zip(pop, param)) do (subj, p)
  return simobs(pk_01, subj, p, obstimes = [10,20,30,40,50,60,70,90,110,150])
end


f, a, p = sim_plot(pk_01, sim, 
        observations = :cp, 
        color = :redsblues,
        linewidth = 4,
        axis = (xlabel = "Time (minute)", 
                ylabel = "PK01 Concentrations (ug/L)",
                xticks = 0:20:160))
axislegend(a) 
f


## Generation of Population - Dataset
par1 = (tvcl = 0.10,
        tvvc = 9.98,
        Ω    = Diagonal([0.04,0.09]),
        σ    = 20.80)

par2 = (tvcl = 0.20,
        tvvc = 9.82,
        Ω    = Diagonal([0.15,0.0225]),
        σ    = 27.46)

par3 = (tvcl = 0.20,
        tvvc = 10.22,
        Ω    = Diagonal([0.04,0.09]),
        σ    = 8.78)

par4 = (tvcl = 0.20,
        tvvc = 19.95,
        Ω    = Diagonal([0.15,0.0225]),
        σ    = 8.50)

evs2 = DosageRegimen(10, time=0, cmt=1)

pop1 = map(i -> Subject(id=i, events=evs2, covariates=(sex="Male",)), 1:20)
pop2 = map(i -> Subject(id=i, events=evs2, covariates=(sex="Male",)), 21:40)
pop3 = map(i -> Subject(id=i, events=evs2, covariates=(sex="Female",)), 41:60)
pop4 = map(i -> Subject(id=i, events=evs2, covariates=(sex="Female",)), 61:80)

Random.seed!(1234)
pop_sim1 = simobs(pk_01, pop1, par1, obstimes=[10,20,30,40,50,60,70,90,110,150])
Random.seed!(1234)
pop_sim2 = simobs(pk_01, pop2, par2, obstimes=[10,20,30,40,50,60,70,90,110,150])
Random.seed!(1234)
pop_sim3 = simobs(pk_01, pop3, par3, obstimes=[10,20,30,40,50,60,70,90,110,150])
Random.seed!(1234)
pop_sim4 = simobs(pk_01, pop4, par4, obstimes=[10,20,30,40,50,60,70,90,110,150])

sim_pop_all = [pop_sim1;pop_sim2;pop_sim3;pop_sim4]
df_sim = DataFrame(sim_pop_all)

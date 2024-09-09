using Pumas

model = @model begin
    @param begin
        θ ∈ VectorDomain(lower = zeros(5))
        σ ∈ RealDomain(lower = 0.0)
    end
    @pre begin
        CL = θ[1]
        Vc = θ[2]
        Ka = θ[3]
        Vp = θ[4]
        Q = θ[5]
    end
    @dynamics begin
        Depot' = -Ka * Depot
        Central' = Ka * Depot - (CL + Q) / Vc * Central + Q / Vp * Peripheral
        Peripheral' = Q / Vc * Central - Q / Vp * Peripheral
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ + 1e-6)
    end
end

params = (θ = [35, 100, 0.5, 210, 30], σ = 0.1)

skeleton = Subject(
    id = 1,
    time = 0.0:0.5:30.0,
    events = DosageRegimen(3000, time = 0.0, cmt = 1),
    observations = (; dv = nothing),
)

using OptimalDesign

times = OptimalDesign.ObsTimes(skeleton.time)
F = OptimalDesign.fim(model, [skeleton], params, [times], FO())

E = eigen(F)
E.values[1]

using Random

rng = Random.default_rng()
Random.seed!(rng, 12345)
pop = [Subject(simobs(model, skeleton, params; rng))]

ll0 = loglikelihood(model, pop, params, NaivePooled())

d = E.vectors[:, 1]

params1 = (θ = [35, 100, 0.5, 210, 30], σ = 0.1)
fpm1 = fit(model, pop, params1, NaivePooled())

params2 = (θ = [10.0, 300, 1.0, 10.0, 5], σ = 0.2)
fpm2 = fit(model, pop, params2, NaivePooled())

hcat(coef(fpm1).θ, coef(fpm2).θ)

rng = Random.default_rng()
Random.seed!(rng, 12345)
newparams = (; θ = params.θ, σ = 0.2)
newpop = [Subject(simobs(model, skeleton, newparams; rng))]

fpm2 = fit(model, newpop, newparams, NaivePooled())


using Pumas
using PumasUtilities
using Random
using CairoMakie


pk_01 = @model begin
    @metadata begin
        desc = "One Compartment Model"
        timeu = u"minite"
    end
    @param begin
    "Clearance (L/hr)"
    tvcl = RealDomain(lower=0)
    "Volumn (L)"
    tvvc = Realdomain(lower=0)
    Ω = PDiagDomain(2)
    end
end
module Monotones

export relative_entropy
import RTFOOL.Context

function relative_entropy(ps::AbstractVector{Float64}, qs::AbstractVector{Float64})
    if length(ps) != length(qs)
        error("states have different sizes")
    elseif !(sum(ps) ≈ 1.0) || !(sum(qs) ≈ 1.0)
        error("state does not sum to 1.0")
    end
    sum(map((p,q) -> (p > 0.0) ? p*log(p/q) : 0.0, ps, qs))
end
relative_entropy(ctx::Context) = relative_entropy(ctx.system_state, ctx.bath_state)

end

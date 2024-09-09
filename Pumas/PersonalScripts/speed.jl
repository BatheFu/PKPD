using BenchmarkTools

function testp(x)::Int64
	for i=2:x-1
		if x%i==0
			return 0
		end
	end
	return 1
end

function foo()
    for i=1:500000
	testp(i)
    end
end

@benchmark foo()

function testp(x::Int32)
    for i::Int32 = 2:x-1
        if x % i == 0
            return 0
        end
    end
    return 1
end

c = Vector{Int32}(undef, 500000)
@time begin
    for i::Int32 = 1:500000
        c[i] = testp(i)
    end
end
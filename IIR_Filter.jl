# Infinite Impulse Responce Digital Filter (IIR)
# when multiplier_backward = zeros(), Finite Impulse Responce Digital Filter (FIR)
struct IIR_Filter
    multiplier_backward::Array{Float64}
    multiplier_forward::Array{Float64}
end

function processing_IIR(signal, filter::IIR_Filter)
    sig_len = length(signal)
    N = length(filter.multiplier_backward)

    vec = vcat(zeros(N), signal, zeros(N))      # signal with mergin
    for i = N + 1 : sig_len + 2 * N 
        for j = 1 : N 
            vec[i] += filter.multiplier_backward[j] * vec[i - j]
        end
    end

    responce = zeros(sig_len + N) 
    for i = 1 : sig_len + N 
        for j = 0 : N 
            responce[i] += filter.multiplier_forward[j + 1] * vec[i + N - j]
        end
    end

    return responce
end

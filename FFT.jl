# Fast Fourier Transform algorithm

function fft(N, data, window; is_inverse=false)
    N = UInt64(N)
    @assert N <= length(data) "N is larger than data length"
    
    n = 63 - leading_zeros(N)   # n = floor(log2(N))
    @assert N == 1 << n "N must be equal to 2^n"

    formatted = _data_formatting(N, data, window)
    reversed = _bit_reverse_data(N, formatted)
    return _fft!(N, reversed, is_inverse)
end

function inverse_fft(N, data, window)
    return fft(N, data, window; is_inverse=true)
end

# subroutines
function _data_formatting(N, data, window)
    # select window function
    w = if window == :Hanning 
        n -> 0.5 * (1.0 - cos(2π * n / N))
    elseif window == :Hamming
        n -> 0.54 - 0.46 * cos(2π * n / N)
    elseif window == :Blackman
        n -> 0.42 - 0.5 * cos(2π * n / N) + 0.08 * cos(4π * n / N)
    elseif window == :Square
        n -> 1.0
    else
        error("invalid window function. valid functions are:\nwindow = :Hanning, :Hamming, :Blackman, :Square")
    end

    # apply window function
    newdata = zeros(ComplexF64, N)
    for i = 1 : N
        newdata[i] = data[i] * w(i)
    end
    return newdata
end

function _bit_reverse_increment(x, N)
    m = N >> 1
    y = x ⊻ m 
    if y < m 
        while true 
            m = m >> 1
            y = y ⊻ m
            if y >= m 
                break
            end
        end
    end
    return y
end

function _bit_reverse_data(N, data)
    rev = zeros(ComplexF64, N)
    x = 0
    for i = 1 : N
        rev[i] = data[x + 1]
        x = _bit_reverse_increment(x, N)
    end
    return rev
end

function _fft!(N, data, is_inverse)
    ω = is_inverse ? (2π / N) : -(2π / N) 
    W = ones(ComplexF64, N)
    for j = 1 : N - 1
        W[j + 1] = cos(ω * j) + sin(ω * j) * im
    end

    # butterfly 
    n = 63 - leading_zeros(N)   # n = floor(log2(N))
    for step = 1 : n 
        pair_distance = 1 << (step - 1)
        rot_factor = N >> step
        for j = 0 : N - 1 
            if pair_distance & j != 0
                continue
            end

            tmp = data[j + 1]
            Δ = data[j + pair_distance + 1] * W[((j * rot_factor) & (N - 1)) + 1]
            data[j + 1] = tmp + Δ
            data[j + pair_distance + 1] = tmp - Δ 
        end
    end

    if is_inverse 
        data = data ./ N
    end

    return data 
end

# test

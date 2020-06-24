@inline function CUDA.exp(x::ComplexF64)
    scale = CUDA.exp( x.re )
    return ComplexF64( scale * CUDA.cos(x.im), scale * CUDA.sin(x.im) )
end

function kernel_initial_psi(a,N,k_in)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    l = (blockIdx().z - 1) * blockDim().z + threadIdx().z
    if j <= size(a,1)
        a[j , k, l] = CUDA.exp(-im*((j-N/2)*k_in[1]+ (k-N/2)*k_in[2]+(l-N/2)*k_in[3]))
    end
    return nothing
end

function kernel_laplace(a,b)
    len = size(a,1)
    wid = size(a,2)
    ti = threadIdx().x
    tj = threadIdx().y
    i = blockIdx().x * blockDim().x + ti
    j = blockIdx().y * blockDim().y + tj
    
    if i < len && i > 1 && j < wid && j > 1
        a[i + len*j] = 1.0/6.0 * (4.0 * (b[i+1 + len*j] + b[i-1 + len*j] + b[i + len*(j+1)] + b[i + len*(j-1)])
                                 + 1.0 * (b[i+1 + len*(j+1)] + b[i+1 + len*(j-1)] + b[i-1 + len*(j+1)] + b[i-1 + len*(j-1)])
                                 - 20.0 * b[i + len*j])
    end
    return nothing
end

function kernel_laplace2(a,b, T=Float32)
    len = size(a,1)
    wid = size(a,2)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1) * blockDim().y + threadIdx().y
    
    if i > 1 && i < len && j > 1 && j < wid
        a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                                 + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[i+1, j-1])
                                 - T(20.0) * b[i, j])
    else
        a[i,j] = b[i,j]
    end
    return nothing
end

function kernel_laplace3(a,b, T=Float32)
    len = size(a,1)
    wid = size(a,2)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1) * blockDim().y + threadIdx().y
    for i = 1:len, j = 1:wid
        if i > 1 && i < len && j > 1 && j < wid
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                                 + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[i+1, j-1])
                                 - T(20.0) * b[i, j])
        elseif i == len
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                        + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[i+1, j-1])
                        - T(20.0) * b[i, j])
        end
    end
    return nothing
end


function kernel_laplace4(a,b, T=Float32)
    # working good! :)
    len = size(a,1)
    wid = size(a,2)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1) * blockDim().y + threadIdx().y
    for i = 1:len, j = 1:wid
        # within the body
        if i > 1 && i < len && j > 1 && j < wid
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                                 + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[i+1, j-1])
                                 - T(20.0) * b[i, j])
        # at the right horizontal edge
        elseif i == len && j > 1 && j < wid
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                        + T(1.0) * (b[1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[1, j-1])
                        - T(20.0) * b[i, j])
        # at the left horizontal edge
        elseif i == 1 && j > 1 && j < wid
        a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[len, j] + b[i, j+1] + b[i, j-1])
                    + T(1.0) * (b[i+1, j+1] + b[len, j+1] + b[len, j-1] + b[i+1, j-1])
                    - T(20.0) * b[i, j])
        # at the bottom
        elseif j == 1 && i > 1 && i < wid
        a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, wid])
                    + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, wid] + b[i+1, wid])
                    - T(20.0) * b[i, j])

        # at the top
        elseif j == wid && i > 1 && i < wid
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, 1] + b[i, j-1])
                        + T(1.0) * (b[i+1, 1] + b[i-1, 1] + b[i-1, j-1] + b[i+1, j-1])
                        - T(20.0) * b[i, j])
        end
        # Now the four corners
        a[1,1] = T(1.0/6.0) * (T(4.0) * (b[2, 1] + b[len, 1] + b[1, 2] + b[1, wid])
               + T(1.0) * (b[2, 2] + b[len, 2] + b[len, wid] + b[2, wid])
               - T(20.0) * b[1, 1])
        a[len,1] = T(1.0/6.0) * (T(4.0) * (b[1, 1] + b[len, 2] + b[len-1, 1] + b[len, wid])
               + T(1.0) * (b[1, 2] + b[1, wid] + b[len-1, wid] + b[len-1, 2])
               - T(20.0) * b[len, 1])
        a[1,wid] = T(1.0/6.0) * (T(4.0) * (b[2, wid] + b[1, 1] + b[1, wid-1] + b[len, wid])
               + T(1.0) * (b[2, 1] + b[2, wid-1] + b[len, wid-1] + b[len, 1])
               - T(20.0) * b[1, wid])
        a[len,wid] = T(1.0/6.0) * (T(4.0) * (b[len, 1] + b[1, wid] + b[len, wid-1] + b[len-1, wid])
               + T(1.0) * (b[1, 1] + b[1, wid-1] + b[len-1, wid-1] + b[len-1, 1])
               - T(20.0) * b[len, wid])

    end
    return nothing
end

function kernel_laplace_working(a,b, T=Float32)
    len = size(a,1)
    wid = size(a,2)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1) * blockDim().y + threadIdx().y
    for i = 1:len, j = 1:wid
        if i >= 1 && i <= len && j >= 1 && j <= wid
            a[i,j] = T(1.0/6.0) * (T(4.0) * (b[i+1, j] + b[i-1, j] + b[i, j+1] + b[i, j-1])
                                 + T(1.0) * (b[i+1, j+1] + b[i-1, j+1] + b[i-1, j-1] + b[i+1, j-1])
                                 - T(20.0) * b[i, j])
        end                            
    end
    return nothing
end

blocks = (4,4)
threads = (8,8)

# @cuda blocks=blocks threads=threads  kernel_initial_psi(psi_int, N, k)
# @cuda blocks=blocks threads=threads  kernel_laplace(a,b)


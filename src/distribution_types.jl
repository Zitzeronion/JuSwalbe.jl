abstract type Distributionfunction end

mutable struct DistributionD1Q3{T} <: Distributionfunction
    f0::T
    f1::T
    f2::T
end

mutable struct DistributionD2Q9{T} <: Distributionfunction
    f0::T
    f1::T
    f2::T
    f3::T
    f4::T
    f5::T
    f6::T
    f7::T
    f8::T
end

#=
    A Material class describes the cross-sections and material properties
=#
type Material

    # The number of energy groups for the cross-sections
    _num_groups::Int64

    # The total cross-section for each energy group
    _sigma_t::Array{Float64}
    
    # The scattering cross-section for each energy group
    _sigma_s::Array{Float64,2}
    
    # The fission cross-section for each energy group
    _sigma_f::Array{Float64}
    
    # The neutron production for each energy group (the fission cross-section
    # times the average number of neutrons produced per fission)
    _nu_sigma_f::Array{Float64}

    # The propability of a neutron being born from fission in each energy group
    _chi::Array{Float64}

    # The constructor for the Material
    Material(num_groups::Int64) = new(num_groups,
                                      Array(Float64,num_groups),
                                      Array(Float64,(num_groups),num_groups),
                                      Array(Float64,num_groups),
                                      Array(Float64,num_groups),
                                      Array(Float64,num_groups))

end

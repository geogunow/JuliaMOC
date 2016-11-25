#=
    A Segment describes a short characteristic in terms of its length and the
    ID of the FSR it traverses
=#
type Segment

    # The ID of the FSR the Segment traverses
    _id::Int64
    
    # The length of the Segment
    _length::Float64
end

#=
    A Track describes a long characterstic in terms of its starting cooridnates
    and its associated segments (short characteristics)
=#
type Track

    # The starting coordinates of the Track
    _coords::Tuple{Float64, Float64}

    # An array of the segments associated with the Track
    _segments::Array{Segment}

    # The ID of the Track
    _uid::Int64
    
    # The constructor for the Track
    Track() = new((0,0), Array(Segment,0), 0)

end

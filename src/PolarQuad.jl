#=
    Returns the TY quadrature set for a given number of polar angles in
    terms of the sine of the abscissa (polar angles) and the associated weights
=#
function getPolarQuadrature(num_polar_2::Int64)

    if (num_polar_2 > 3 || num_polar_2 < 1)
        println("ERROR: Number of polar angles not supported")
    end

    # Allocate arrays
    sin_thetas = Array(Float64, num_polar_2)
    weights = Array(Float64, num_polar_2)

    # Fill in the sine of the polar angles and the weights for the given number
    # of polar angles
    if num_polar_2 == 1
        sin_thetas[1] = 0.798184
        weights[1] = 1.0 / 2
    elseif num_polar_2 == 2
        sin_thetas[1] = 0.363900
        sin_thetas[2] = 0.899900
        weights[1] = 0.212854 / 2
        weights[2] = 0.787146 / 2
    elseif num_polar_2 == 3
        sin_thetas[1] = 0.166648
        sin_thetas[2] = 0.537707
        sin_thetas[3] = 0.932954
        weights[1] = 0.046233 / 2
        weights[2] = 0.283619 / 2
        weights[3] = 0.670148 / 2
    else
        println("ERROR: Only 1, 2, and 3 polar angles are supported in [0, π]")
    end

    return sin_thetas, weights
end

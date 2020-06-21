using ReachabilityAnalysis, Plots

const ne3x = [-1.0, 0.0, 0.0]
const ne3y = [0.0, -1.0, 0.0]
const ne3z = [0.0, 0.0, -1.0]

const positive_orthant = HPolyhedron([HalfSpace(ne3x, 0.0),   # x >= 0
                                    HalfSpace(ne3y, 0.0),   # y >= 0
                                    HalfSpace(ne3z, 0.0)])  # z >= 0

@taylorize function prod_dest_1!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    return du
end

@taylorize function prod_dest_2!(du, u, params, t)
    x, y, z, a = u[1], u[2], u[3], u[4]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    du[4] = zero(x)
    return du
end


@inline function prod_dest_property(solz)
    X = project(solz(100.0), vars=(1, 2, 3))

    # check that all variables are nonnegative
    nonnegative = X ⊆ positive_orthant

    # compute the volume of the last reach-set
    H = overapproximate(X, Hyperrectangle)

    # check that that 10.0 belongs to the minkowski sum of the reach-sets projected in each coordinate
    B = convert(IntervalBox, H) # get the product-of-intervals representation

    return nonnegative && (10 ∈ sum(B)), volume(H)
end


function production_destruction(; case=1)

    if case == 1
        X0=(9.5 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01)
        prob = @ivp(x'= prod_dest_1!(x), dim:3, x(0) ∈ X0)
    elseif case == 2
        X0=(9.98 .. 9.98) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.296 .. 0.304)
        prob = @ivp(x'= prod_dest_2!(x), dim:4, x(0) ∈ X0)
    elseif case == 3
        X0 = (9.7 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.298 .. 0.302)
        prob = @ivp(x'= prod_dest_2!(x), dim:4, x(0) ∈ X0)
    else
        error("the case = $case is not implemented")
    end

    return prob
end

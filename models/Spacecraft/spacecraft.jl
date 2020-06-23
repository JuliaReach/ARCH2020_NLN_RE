const μ = 3.986e14 * 60^2
const r = 42164.0e3
const r² = r^2
const mc = 500.0
const n² = μ / r^3
const n = sqrt(n²)

const two_n = 2*n
const μ_r² = μ/r²

# columns correspond to x, y, vx, vy, t
const K₁ = [-28.8287 0.1005 -1449.9754 0.0046 0.0;
            -0.087 -33.2562 0.00462 -1451.5013 0.0]
const K₂ = [-288.0288 0.1312 -9614.9898 0.0 0.0;
            -0.1312 -288.0 0.0 -9614.9883   0.0]

const K₁mc = K₁/mc
const K₂mc = K₂/mc

function mymul!(v, A, x)
    @inbounds for ind = 1:length(v)
        v[ind] = zero(x[1])
        for jind = 1:length(x)
            v[ind] += A[ind, jind] * x[jind]
        end
    end
    return nothing
end

# dynamics in the 'approaching' mode
@taylorize function spacecraft_approaching!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    uxy = Vector{typeof(x)}(undef, 2)
    mymul!(uxy, K₁mc, u)

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    # t' = 1
    du[5] = one(x)

    return du
end

# dynamics in the 'rendezvous attempt' mode
@taylorize function spacecraft_attempt!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    uxy = Vector{typeof(x)}(undef, 2)
    mymul!(uxy, K₂mc, u)

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    # t' = 1
    du[5] = one(x)

    return du
end

# dynamics in the 'aborting' mode
@taylorize function spacecraft_aborting!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x)
    du[3] = (n²*x + two_n*vy) + (μ_r² - μ_rc³*rx)

    # vy' = n²y - 2n*vx - μ/(rc^3)y
    du[4] = (n²*y - two_n*vx) - μ_rc³*y

    # t' = 1
    du[5] = one(x)

    return du
end

function spacecraft(; X0 = Hyperrectangle([-900., -400., 0., 0., 0.],
                                          [25., 25., 0., 0., 0.]),
                         init=[(1, X0)],
                         abort_time=(120.0, 150.0))

    # variables
    x = 1   # x position
    y = 2   # y position
    vx = 3  # x velocity
    vy = 4  # y velocity
    t = 5   # time
    n = 4 + 1  # number of variables
    t_abort_lower, t_abort_upper = abort_time[1], abort_time[2]

    automaton = LightAutomaton(3)

    # mode 1 "approaching"
    invariant = HalfSpace(sparsevec([x], [1.], n), -100.) # x <= -100
    approaching = @system(x' = spacecraft_approaching!(x), dim:5, x ∈ invariant)

    # mode 2 ("rendezvous attempt")
    invariant = HalfSpace(sparsevec([x], [-1.], n), 100.) # x >= -100
    attempt = @system(x' = spacecraft_attempt!(x), dim:5, x ∈ invariant)

    # mode 3 "aborting"
    invariant = Universe(n)
    aborting = @system(x' = spacecraft_aborting!(x), dim: 5, x ∈ invariant)

    # transition "approach" -> "attempt"
    add_transition!(automaton, 1, 2, 1)
    guard = HalfSpace(sparsevec([x], [-1.], n), 100.) # x >= -100
    t1 = @map(x -> x, dim: n, x ∈ guard)


    # transition "approach" -> "abort"
    add_transition!(automaton, 1, 3, 2)
    guard_time = HPolyhedron([HalfSpace(sparsevec([t], [-1.], n), -t_abort_lower),  # t >= t_abort_lower
                              HalfSpace(sparsevec([t], [1.], n), t_abort_upper)])   # t <= t_abort_upper
    t2 = @map(x -> x, dim: n, x ∈ guard_time)

    # transition "attempt" -> "abort"
    add_transition!(automaton, 2, 3, 3)
    t3 = @map(x -> x, dim: n, x ∈ guard_time)

    H = HybridSystem(automaton=automaton,
                     modes=[approaching, attempt, aborting],
                     resetmaps=[t1, t2, t3])

    return InitialValueProblem(H, init)
end

# variables
const x, y, vx, vy, t = 1:5
const numvars = 5
const tan30 = tand(30)

LineOfSightCone = HPolyhedron([HalfSpace(sparsevec([x], [-1.], numvars), 100.),   # x >= -100
                    HalfSpace(sparsevec([x, y], [tan30, -1.], numvars), 0.),  # y >= x tan(30°)
                    HalfSpace(sparsevec([x, y], [tan30, 1.], numvars), 0.),   # -y >= x tan(30°)
                   ])

Target = Hyperrectangle(zeros(2), [0.2, 0.2]);

function line_of_sight(sol)
    all_idx = findall(x -> x == 2, location.(sol))  # attempt
    for idx in all_idx
        verif = all(set(R) ⊆ LineOfSightCone for R in sol[idx][2:end])
        !verif && return false
    end
    return true
end

function absolute_velocity(R::RA.AbstractLazyReachSet)
    vx = 3
    vy = 4
    vx2 = set(overapproximate(project(R, vars=(vx,)), Interval)).dat
    vy2 = set(overapproximate(project(R, vars=(vy,)), Interval)).dat
    sqrt(vx2 + vy2)
end

function absolute_velocity(R::RA.TaylorModelReachSet)
    Z = overapprximate(R, Zonotope)
    absolute_velocity(Z)
end

function velocity_constraint(sol)
    all_idx = findall(x -> x == 2, location.(sol)) # attempt
    for idx in all_idx
        # maximum velocity measured in m/min
        verif = all(absolute_velocity(R) < 0.055 * 60. for R in sol[idx])
        !verif && return false
    end
    return true
end

function target_avoidance(sol)
    all_idx = findall(x -> x == 3, location.(sol)) # aborting
    for idx in all_idx
        verif = all(is_intersection_empty(set(Projection(R, [x, y])), Target) for R in sol[idx])
        !verif && return false
    end
    return true
end

export seg_mean
export seg_extract
export seg_select

"""
    seg_mean(seg)

Calculate mean of a segment (e.g. spectrogram).

# Arguments

  - `seg::AbstractArray`

# Returns

  - `sm::Vector{Float64}`: averaged segment
"""
function seg_mean(seg::AbstractArray)::Vector{Float64}

    _chk3d(seg)
    sm = reshape(mean(mean(seg; dims = 1), dims = 2), size(seg, 3))

    return sm

end

"""
    seg2_mean(seg1, seg2)

Calculate mean of two segments (e.g. spectrograms).

# Arguments

  - `seg1::AbstractArray`
  - `seg2::AbstractArray`

# Returns

Named tuple containing:

  - `seg1::Vector{Float64}`: averaged segment 1
  - `seg2::Vector{Float64}`: averaged segment 2
"""
function seg_mean(seg1::AbstractArray, seg2::AbstractArray)::@NamedTuple{seg1::Vector{Float64}, seg2::Vector{Float64}}

    seg1 = seg_mean(seg1)
    seg2 = seg_mean(seg2)

    return (seg1 = seg1, seg2 = seg2)

end

"""
    seg_extract(m, rc; <keyword arguments>)

Extract segment from a matrix.

# Arguments

  - `m::AbstractMatrix`
  - `rc::NTuple{4, Int64}`: upper-left corner row and column, bottom-right corner row and column
  - `c::Bool=false`: if true, use circular segment; for circular segment the segment is always returned as vector
  - `v::Bool=false`: if true, return as vector (matrix m by rows over columns)

# Returns

  - `seg::Union{AbstractMatrix, AbstractVector}`
"""
function seg_extract(
        m::AbstractMatrix, rc::NTuple{4, Int64}; v::Bool = false, c::Bool = false
    )::Union{AbstractMatrix, AbstractVector}

    r1 = rc[1]
    c1 = rc[2]
    r2 = rc[3]
    c2 = rc[4]

    @assert r1 > 0 "r1 must be > 0."
    @assert r2 > 0 "r2 must be > 0."
    @assert c1 > 0 "c1 must be > 0."
    @assert c2 > 0 "c2 must be > 0."
    @assert r1 <= size(m, 1) "r1 must be ≤ $(size(m, 1))."
    @assert c1 <= size(m, 2) "r2 must be ≤ $(size(m, 2))."
    @assert r2 <= size(m, 1) "c1 must be ≤ $(size(m, 1))."
    @assert c2 <= size(m, 2) "c2 must be ≤ $(size(m, 2))."

    if !c
        seg = !v ? m[r1:r2, c1:c2] : vec(m[r1:r2, c1:c2])
    else
        seg = zeros(Bool, size(m))
        seg_radius = distance((r1, c1), (r2, c2))
        for idx_r in axes(m, 1), idx_c in axes(m, 2)
            if distance((r1, c1), (idx_r, idx_c)) <= seg_radius
                seg[idx_r, idx_c] = true
            end
        end

        seg = m[seg .== true]
    end

    return seg

end

"""
    seg_extract(m; <keyword arguments>)

Extract segment from a matrix using thresholding.

# Arguments

  - `m::AbstractMatrix`
  - `threshold::Union{Real, Tuple{Real, Real}}=0`: threshold
  - `threshold_type::Symbol=:neq`: rule for thresholding:
      + `:eq`: return equal to threshold
      + `:neq`: return not equal to threshold
      + `:geq`: return ≥ to threshold
      + `:leq`: return ≤ to threshold
      + `:g`: return > to threshold
      + `:l`: return < to threshold
      + `:in`: draw region is values are in the threshold values, including threshold boundaries
      + `:bin`: draw region is values are between the threshold values, excluding threshold boundaries

# Returns

Named tuple containing:

  - `idx::Vector{CartesianIndex{2}}`: Cartesian coordinates of matrix elements
  - `bm::Matrix{Bool}`: map of the segment
"""
function seg_extract(
        m::AbstractMatrix; threshold::Union{Real, Tuple{Real, Real}} = 0, threshold_type::Symbol = :neq
    )::@NamedTuple{idx::Vector{CartesianIndex{2}}, bm::Matrix{Bool}}

    _check_var(threshold_type, [:eq, :neq, :geq, :leq, :g, :l, :in, :bin], "threshold_type")

    if threshold_type in [:eq, :neq, :geq, :leq, :g, :l]
        @assert length(threshold) == 1 "threshold must contain a single value."
    else
        @assert length(threshold) == 2 "threshold must contain two values."
        _check_tuple(threshold, extrema(m), "threshold")
    end

    if threshold_type === :eq
        idx = findall(x -> x == threshold, m)
    elseif threshold_type === :neq
        idx = findall(x -> x != threshold, m)
    elseif threshold_type === :geq
        idx = findall(x -> x >= threshold, m)
    elseif threshold_type === :leq
        idx = findall(x -> x <= threshold, m)
    elseif threshold_type === :g
        idx = findall(x -> x > threshold, m)
    elseif threshold_type === :l
        idx = findall(x -> x < threshold, m)
    elseif threshold_type === :in
        idx = findall(x -> (x >= threshold[1] && x <= threshold[2]), m)
    elseif threshold_type === :bin
        idx = findall(x -> (x > threshold[1] && x < threshold[2]), m)
    end

    bm = zeros(Bool, size(m))
    bm[idx] .= true

    return (idx = idx, bm = bm)

end

export seg_select

"""
    seg_select(m; <keyword arguments>)

Interactive selection of a matrix area.

# Arguments

  - `m::AbstractMatrix`
  - `shape::Symbol=:r`: selection shape:
      + `:r`: rectangular
      + `:p`: point
      + `:c`: circular
  - `extract::Bool=false`: if true, return values of the matrix
  - `v::Bool=false`: if true, return as vector (matrix m by rows over columns), always true if `shape=:c`

# Returns

  - `seg::Union{Nothing, <:Real, Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64, Int64}, Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}}}`: extracted segment or its coordinates
"""
function seg_select(
        m::AbstractMatrix; shape::Symbol = :r, extract::Bool = false, v::Bool = false
    )::Union{
        Nothing,
        <:Real,
        Tuple{Int64, Int64},
        Tuple{Int64, Int64, Int64, Int64},
        Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}},
    }

    _check_var(shape, [:r, :p, :c], "shape")

    size_x = size(m, 2)
    size_y = size(m, 1)

    p = GLMakie.Figure(size = (size_x, size_y))
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = "",
        ylabel = "",
        title = "",
        aspect = DataAspect(),
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    hidedecorations!(ax)
    hm = GLMakie.heatmap!(m[end:-1:1, :]', colormap = :darktest)

    poins = nothing
    if shape in [:p, :r]
        points = Observable(Point2i[])
    else
        points = Observable(Point2i(div(size_x, 2), div(size_y, 2)))
    end
    radius = Observable(1)

    if shape === :p

        GLMakie.scatter!(ax, points, marker = :rect, markersize = 10, color = :red)

        on(events(ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                pos = round.(Int64, mouseposition(ax))
                if length(points[]) == 0
                    push!(points[], pos)
                    notify(points)
                elseif length(points[]) == 1
                    pop!(points[])
                    push!(points[], pos)
                    notify(points)
                end
            end
            if event.button == Mouse.right && event.action == Mouse.press
                if length(points[]) > 0
                    pop!(points[])
                    notify(points)
                end
            end
        end

    elseif shape === :c

        GLMakie.arc!(ax, points, radius, -pi, pi, linewidth = 5, color = :red)

        on(events(p).scroll, priority = 1) do (dx, dy)
            if dy == 1.0
                radius[] <= size_x && (radius[] += dy)
            elseif dy == -1.0
                radius[] > 1 && (radius[] += dy)
            end
            notify(radius)
        end

        on(events(ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                pos = round.(Int64, mouseposition(ax))
                points[] = Point2i(pos)
                notify(points)
            end
            if event.button == Mouse.right && event.action == Mouse.press
                points[] = Point2i(0, 0)
                radius[] = 1
                notify(points)
            end
        end

    elseif shape === :r

        GLMakie.lines!(ax, points, linewidth = 5, color = :red)
        on(events(ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                pos = round.(Int64, mouseposition(ax))
                if length(points[]) == 0
                    push!(points[], pos)
                    push!(points[], pos)
                    push!(points[], pos)
                    push!(points[], pos)
                    notify(points)
                elseif length(points[]) == 1
                    push!(points[], Point2i(points[][end][1], pos[2]))
                    push!(points[], pos)
                    push!(points[], Point2i(pos[1], points[][end - 2][2]))
                    push!(points[], Point2i(points[][end - 3]))
                    notify(points)
                elseif length(points[]) == 4
                    pop!(points[])
                    pop!(points[])
                    pop!(points[])
                    push!(points[], Point2i(points[][end][1], pos[2]))
                    push!(points[], pos)
                    push!(points[], Point2i(pos[1], points[][end - 2][2]))
                    push!(points[], Point2i(points[][end - 3]))
                    notify(points)
                elseif length(points[]) == 5
                    pop!(points[])
                    pop!(points[])
                    pop!(points[])
                    pop!(points[])
                    push!(points[], Point2i(points[][end][1], pos[2]))
                    push!(points[], pos)
                    push!(points[], Point2i(pos[1], points[][end - 2][2]))
                    push!(points[], Point2i(points[][end - 3]))
                    notify(points)
                end
            end
            if event.button == Mouse.right && event.action == Mouse.press
                for _ in length(points[]):-1:1
                    pop!(points[])
                end
                notify(points)
            end
        end
    end

    on(events(p).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.enter
                close(display(p))
            end
        end
    end

    wait(display(p))

    if length(points[]) == 0
        return nothing
    elseif length(points[]) == 1
        c = points[][1][1]
        r = size_y - points[][1][2]
        c < 1 && (c = 1)
        r < 1 && (r = 1)
        c > size_x && (c = size_x)
        r > size_y && (r = size_y)
    elseif length(points[]) == 2
        if radius[] == 0
            return nothing
        else
            c1 = points[][1]
            r1 = size_y - points[][2]
            c1 < 1 && (c1 = 1)
            r1 < 1 && (r1 = 1)
            c1 > size_x && (c1 = size_x)
            r1 > size_y && (r1 = size_y)

            c2 = c1 + radius[]
            r2 = r1 + radius[]
            c2 < 1 && (c2 = 1)
            r2 < 1 && (r2 = 1)
            c2 > size_x && (c2 = size_x)
            r2 > size_y && (r2 = size_y)
        end
    elseif length(points[]) == 5
        c1 = points[][1][1]
        r1 = size_y - points[][1][2]
        c1 < 1 && (c1 = 1)
        r1 < 1 && (r1 = 1)
        c1 > size_x && (c1 = size_x)
        r1 > size_y && (r1 = size_y)

        c2 = points[][3][1]
        r2 = size_y - points[][3][2]
        c2 < 1 && (c2 = 1)
        r2 < 1 && (r2 = 1)
        c2 > size_x && (c2 = size_x)
        r2 > size_y && (r2 = size_y)
    end

    if shape in [:r, :c]
        r1 > r2 && ((r1, r2) = _swap(r1, r2))
        c1 > c2 && ((c1, c2) = _swap(c1, c2))
        c = shape == :c
        return !extract ? (r1, c1, r2, c2) : seg_extract(m, (r1, c1, r2, c2), v = v, c = c)
    else
        return !extract ? (r, c) : m[r, c]
    end

end

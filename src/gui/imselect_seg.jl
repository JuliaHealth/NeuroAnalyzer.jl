export imselect_seg

"""
    imselect_seg(m; <keyword arguments>)

Interactive selection of matrix area.

# Arguments

- `m::AbstractMatrix`
- `shape::Symbol=:r`: selection shape:
    - `:r`: rectangular
    - `:p`: point
- `extract::Bool=false`: if true, return values of the matrix
- `v::Bool=false`: if true, return as vector (matrix m by rows over columns)

# Returns

- `r1::Int64`: upper-left corner (row)
- `r2::Int64`: bottom-right corner (row)
- `c1::Int64`: upper-left corner (column)
- `c2::Int64`: bottom-right corner (column)

or

- `seg::Union{Nothing, <:Real, Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64, Int64}, Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}}}`: extracted segment or its coordinates
"""
function imselect_seg(m::AbstractMatrix; shape::Symbol=:r, extract::Bool=false, v::Bool=false)::Union{Nothing, <:Real, Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64, Int64}, Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}}}

    _check_var(shape, [:r, :p], "shape")

    size_x = size(m, 2)
    size_y = size(m, 1)

    p = GLMakie.Figure(size=(size_x, size_y))
    ax = GLMakie.Axis(p[1, 1],
                      xlabel="",
                      ylabel="",
                      title="",
                      aspect=DataAspect(),
                      xticksvisible=false,
                      yticksvisible=false,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0),
                      xzoomlock=true,
                      yzoomlock=true,
                      xpanlock=true,
                      ypanlock=true,
                      xrectzoom=false,
                      yrectzoom=false)
    hidedecorations!(ax)
    hm = GLMakie.heatmap!(m[end:-1:1, :]',
                          colormap=:darktest)

    points = Observable(Point2i[])

    if shape === :p

        GLMakie.scatter!(ax,
                         points,
                         marker=:rect,
                         markersize=10,
                         color=:red)

        on(events(p).mousebutton) do event
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

    elseif shape === :r

        GLMakie.lines!(ax,
                       points,
                       linewidth=5,
                       color=:red)
        on(events(p).mousebutton) do event
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

    if shape === :r
        r1 > r2 && ((r1, r2) = _swap(r1, r2))
        c1 > c2 && ((c1, c2) = _swap(c1, c2))
        c = shape == :c
        return !extract ? (r1, c1, r2, c2) : seg_extract(m, (r1, c1, r2, c2), v=v, c=c)
    else
        return !extract ? (r, c) : m[r, c]
    end

end

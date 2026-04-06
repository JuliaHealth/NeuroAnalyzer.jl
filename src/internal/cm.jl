"""Print confusion matrix."""
function _cm(cm::Matrix{Int64})::Nothing
    cm_len = length.(string.(cm)) .+ 1
    ml = maximum(cm_len)
    vbar = repeat("─", ml)
    println("              $(lpad("", ml, " ")) group")
    println("                $(lpad("0", div(ml, 2) + 1, " "))   $(lpad("1", ml, " "))")
    print("              ┌─")
    print(vbar)
    print("─┬─")
    print(vbar)
    println("─┐")
    println("            0 │ $(lpad(cm[1], ml, " ")) │ $(lpad(cm[3], ml, " ")) │")
    print(" prediction   ├─")
    print(vbar)
    print("─┼─")
    print(vbar)
    println("─┤")
    println("            1 │ $(lpad(cm[2], ml, " ")) │ $(lpad(cm[4], ml, " ")) │")
    print("              └─")
    print(vbar)
    print("─┴─")
    print(vbar)
    println("─┘")
    return nothing
end

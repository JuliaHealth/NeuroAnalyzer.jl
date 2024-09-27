function _cm(cm::Matrix{Int64})::Nothing
    println("""
                     group
                    0     1
                 ┌─────┬─────┐
               0 │ $(lpad(cm[1], 3, " ")) │ $(lpad(cm[3], 3, " ")) │
    prediction   ├─────┼─────┤
               1 │ $(lpad(cm[2], 3, " ")) │ $(lpad(cm[4], 3, " ")) │
                 └─────┴─────┘
             """)
    return nothing
end
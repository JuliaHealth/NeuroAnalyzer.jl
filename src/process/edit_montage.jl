export edit_montage

"""
    edit_montage(file_name; <keyword arguments>)

Edit montage file in the OS editor.

# Arguments

- `file_name::String`: name of the file to load

# Returns

Nothing
"""
function edit_montage(file_name::String)::Nothing

    @assert isfile(file_name) "File $file_name cannot be loaded."

    if Sys.iswindows()
        run(`notepad.exe $file_name`)
    elseif Sys.islinux()
        editor = ENV["EDITOR"]
        @assert editor != "" "Default editor not set in the system \$EDITOR variable."
        run(`$editor $file_name`)
    elseif Sys.isapple()
        run(`open -e $file_name`)
    else
        error("Couldn't open $file_name")
    end

    return nothing

end
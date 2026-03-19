export edit_montage

"""
    edit_montage(file_name; <keyword arguments>)

Open a montage file in the operating-system default text editor.

- **Windows**: opens with `notepad.exe`.
- **Linux**: opens with the editor specified in the `\$EDITOR` environment variable.
- **macOS**: opens with TextEdit via `open -e`.

# Arguments

- `file_name::String`: path to the montage file to edit; must exist

# Returns

- `Nothing`

# Throws

- `ArgumentError`: if `file_name` does not exist or the `\$EDITOR` variable is not set on Linux
- `ErrorException`: if the operating system is not Windows, Linux, or macOS
"""
function edit_montage(file_name::String)::Nothing

    isfile(file_name) || throw(ArgumentError("File '$file_name' not found."))

    if Sys.iswindows()
        run(`notepad.exe $file_name`)
    elseif Sys.islinux()
        editor = get(ENV, "EDITOR", "")
        isempty(editor) && throw(ArgumentError("No editor set; please define the \$EDITOR environment variable."))
        run(`$editor $file_name`)
    elseif Sys.isapple()
        run(`open -e $file_name`)
    else
        throw(ErrorException("Unsupported operating system; cannot open '$file_name'."))
    end

    return nothing

end

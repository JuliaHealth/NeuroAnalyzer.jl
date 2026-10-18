export test_run

"""
    test_run(file_name)

Run a test described in a JSON file: a full-screen window shows the test objects in order.

The file holds an array of objects, each with a `"type"`:

- `"delay"`: wait `"duration"` milliseconds
- `"image"`: show the image `"file"`

Press `ESC` to close the window. Work in progress.

# Arguments

- `file_name::String`: JSON file with the test structure

# Returns

- `Nothing`

# Example test file

```json
[
    {"type": "image", "file": "/path/to/fixation.png"},
    {"type": "delay", "duration": 500},
    {"type": "image", "file": "/path/to/stimulus.png"},
    {"type": "delay", "duration": 1000}
]
```
"""
function test_run(file_name::String)::Nothing
    _wip()

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    # read and validate the test structure
    test_structure = try
        JSON.parsefile(file_name)
    catch err
        throw(ArgumentError("File $file_name cannot be parsed: $err"))
    end
    test_structure isa AbstractVector ||
        throw(ArgumentError("$file_name must contain an array of test objects."))
    for (idx, object) in enumerate(test_structure)
        object isa AbstractDict && haskey(object, "type") ||
            throw(ArgumentError("Test object $idx must be an object with a \"type\"."))
        if object["type"] == "delay"
            d = object["duration"]
            # duration may be given as a number or a numeric string
            d = d isa Real ? Float64(d) : parse(Float64, string(d))
            d >= 0 || throw(ArgumentError("Test object $idx: duration must be ≥ 0."))
            object["duration"] = d
        elseif object["type"] == "image"
            isfile(object["file"]) ||
                throw(ArgumentError("Test object $idx: file $(object["file"]) not found."))
        else
            throw(ArgumentError("Test object $idx: unknown type \"$(object["type"])\"."))
        end
    end

    css = """
    window {
        background-color: black;
    }
    """
    cssprov = GtkCssProvider(css)

    win = GtkWindow("", 0, 0)
    Gtk4.decorated(win, false)
    Gtk4.cursor(win, "none")
    # attach the stylesheet to the window's display
    push!(Gtk4.display(win), cssprov)

    eck = GtkEventControllerKey(win)
    signal_connect(eck, "key-pressed") do controller, keyval, keycode, state
        # ESC
        keyval == 0xff1b && close(win)
        return nothing
    end

    Gtk4.show(win)
    Gtk4.fullscreen(win)

    @async begin
        try
            for object in test_structure
                if object["type"] == "delay"
                    sleep(object["duration"] / 1000)
                elseif object["type"] == "image"
                    win.child = GtkPicture(object["file"])
                end
            end
            close(win)
        catch e
            @error "Test run failed!" exception = (e, catch_backtrace())
        end
    end

    if !isinteractive()
        Gtk4.GLib.waitforsignal(win, :close_request)
    end

    return nothing
end

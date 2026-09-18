export test_run

# file_name = homedir() * "/Documents/Code/test.json"

function test_run(file_name::String)::Nothing
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    # read test structure
    test_structure = JSON.parsefile(file_name; dicttype = Dict, inttype = Int64, use_mmap = true)

    # parse test structure
    # test for all required properties
    # test properties data types and ranges
    # test if external objects are available

    # preload external objects

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
        if keyval == 0xff1b
            close(win)
        end
    end

    show(win)
    fullscreen(win)

    @async begin
        try
            for structure_idx = eachindex(test_structure)
                object = test_structure[structure_idx]
                object_type = object[:type]
                if object_type == "delay"
                    delay_duration = parse(Float64, object[:duration])
                    sleep(delay_duration / 1000)
                    continue
                end
                if object_type == "image"
                    img = object[:file]
                    win.child = GtkPicture(img)
                    continue
                end
            end
            close(win)
        catch e
            @error "Timer task failed!" exception=(e, catch_backtrace())
        end
    end

    if !isinteractive()
        Gtk4.GLib.waitforsignal(win, :close_request)
    end

    return nothing
end
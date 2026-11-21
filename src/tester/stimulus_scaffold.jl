"""
stimulus_scaffold.jl

Minimal scaffold for presenting visual stimuli with frame-locked timing
and keyboard response capture, using GLMakie (native OpenGL window).

This is a starting point, not a framework — extend the pieces you need:
  - stimulus generation (images, gratings, shapes)
  - audio (pair with PortAudio.jl for sound onset)
  - trial/condition structure, logging to disk
  - staircase / adaptive procedures

Run with:
    julia --project stimulus_scaffold.jl
"""

using GLMakie
using Dates

# ---------------------------------------------------------------------------
# 1. Window setup
# ---------------------------------------------------------------------------

"""
Query the primary monitor's resolution via GLFW (accessible through
GLMakie's own GLFW dependency, so no extra package needed).
"""
function primary_monitor_resolution()
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    vidmode = GLMakie.GLFW.GetVideoMode(monitor)
    return (Int(vidmode.width), Int(vidmode.height))
end

"""
Create a fullscreen GLMakie window and return (figure, screen, axis).

NOTE: GLMakie's `fullscreen=true` kwarg is unreliable across versions
(known open issue — often silently ignored), and a borderless *windowed*
fullscreen (decorated=false + monitor-sized window) is still treated as
an ordinary window by some compositors (e.g. KDE/KWin), so panels/docks
can stay stacked on top of it. The robust fix is to request *true*
exclusive fullscreen from GLFW directly via `SetWindowMonitor`, which
changes the window's mode at the compositor level rather than just its
size/decorations — this is what actually goes above panels.

Set `fullscreen=false` while developing so you keep window decorations
and a manageable size.
"""
function make_window(; fullscreen::Bool=true, bg=:gray)
    GLMakie.activate!(decorated=!fullscreen, focus_on_show=true)

    size_px = fullscreen ? primary_monitor_resolution() : (1200, 800)

    fig = Figure(; size=size_px, backgroundcolor=bg, figure_padding=0)
    ax = Axis(fig[1, 1],
        backgroundcolor=bg,
        aspect=DataAspect(),
    )
    hidedecorations!(ax)   # no axis ticks/labels — this is a stimulus canvas
    hidespines!(ax)

    screen = display(fig)

    if fullscreen
        try
            window = GLMakie.to_native(screen)
            monitor = GLMakie.GLFW.GetPrimaryMonitor()
            vidmode = GLMakie.GLFW.GetVideoMode(monitor)
            # Handing GLFW a monitor here switches the window into real
            # exclusive fullscreen (not just "big and borderless").
            GLMakie.GLFW.SetWindowMonitor(
                window, monitor, 0, 0,
                vidmode.width, vidmode.height, vidmode.refreshrate,
            )
        catch e
            @warn "Could not switch to true GLFW fullscreen; falling back to borderless window. Panels/docks may still appear on top." exception=e
        end
    else
        try
            window = GLMakie.to_native(screen)
            monitor = GLMakie.GLFW.GetPrimaryMonitor()
            vidmode = GLMakie.GLFW.GetVideoMode(monitor)
            x = (vidmode.width - size_px[1]) ÷ 2
            y = (vidmode.height - size_px[2]) ÷ 2
            GLMakie.GLFW.SetWindowPos(window, x, y)
        catch e
            @warn "Could not center window." exception=e
        end
    end

    return fig, screen, ax
end

"""
Close the experiment window and stop its render loop. Always call this
when you're done — GLMakie doesn't do this automatically when a script
finishes, since normally it's left open for interactive/REPL use.
"""
function close_window!(screen)
    try
        GLMakie.destroy!(screen)
    catch e
        @warn "Error while closing window" exception=e
    end
end

# ---------------------------------------------------------------------------
# 2. Frame-locked presentation
# ---------------------------------------------------------------------------

"""
Block until the next actual screen render, then return the timestamp
at which it completed.

GLMakie doesn't expose a raw vsync/swap-buffer callback the way
Psychtoolbox does, and its internal `render_frame`/`to_native` functions
are private and change signature across versions — don't call them
directly. Instead, GLMakie already runs its own render loop in the
background that redraws whenever an Observable in the scene changes.
`screen.render_tick` is notified once per actual render, so we
subscribe to it, then wait for the *next* tick to fire — that tick
corresponds to the render that includes whatever change triggered it.
This is a practical proxy for onset time, not a hardware-verified
timestamp; see the note on photodiode validation at the end of this
file.
"""
function wait_for_next_render(screen)
    result = Ref{Union{Nothing,Float64}}(nothing)
    listener = on(screen.render_tick) do _
        result[] = time()
    end
    while result[] === nothing
        sleep(0.0005)
        yield()
    end
    off(listener)
    return result[]
end

"""
Present a stimulus by mutating a Makie plot's observable(s), then
wait for the resulting render to actually happen and return that
onset timestamp. Keep `update_fn!` as cheap as possible (ideally just
toggling a pre-loaded Observable) so the next tick reflects your change
promptly.
"""
function show_stimulus!(screen, update_fn!::Function)
    update_fn!()                     # e.g. stim_obs[] = new_image
    return wait_for_next_render(screen)
end

# ---------------------------------------------------------------------------
# 3. Keyboard input with timestamps
# ---------------------------------------------------------------------------

"""
Simple response collector. Registers a keyboard callback on `fig.scene`
and records (key, timestamp) pairs into `responses` until `deadline`
(wall-clock seconds, from `time()`) or until a key in `accepted_keys`
is pressed (whichever comes first). Returns the first accepted
response as (key, rt_seconds) relative to `onset_time`, or `nothing`
on timeout.
"""
function wait_for_response(fig, onset_time::Float64;
                            accepted_keys::Vector{Keyboard.Button}=[Keyboard.left, Keyboard.right],
                            timeout::Float64=2.0)
    result = Ref{Union{Nothing, Tuple{Keyboard.Button,Float64}}}(nothing)
    deadline = onset_time + timeout

    listener = on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press && event.key in accepted_keys
            rt = time() - onset_time
            result[] = (event.key, rt)
        end
    end

    # Block until we get a response or hit the deadline.
    # `yield()` lets GLMakie's event loop keep processing input/render events.
    while result[] === nothing && time() < deadline
        sleep(0.001)
        yield()
    end

    off(listener)   # always clean up the callback
    return result[]
end

# ---------------------------------------------------------------------------
# 4. Example trial loop
# ---------------------------------------------------------------------------

"""
Runs a toy experiment: on each trial, show a colored square (the
"stimulus"), wait for a left/right keypress, log RT and accuracy,
then show a blank inter-trial interval.

Replace `stim_obs` / the drawn primitive with your real stimuli
(images via `image!`, gratings via `heatmap!` of a generated matrix,
text via `text!`, etc.).
"""
function run_experiment(; n_trials::Int=5, fullscreen::Bool=true)
    fig, screen, ax = make_window(; fullscreen=fullscreen)

    # Pre-allocate the stimulus as an Observable so updates are cheap.
    stim_color = Observable(:gray)
    poly!(ax, Rect(-0.5, -0.5, 1.0, 1.0); color=stim_color)
    xlims!(ax, -1, 1); ylims!(ax, -1, 1)

    wait_for_next_render(screen)  # let the initial blank/gray screen actually draw

    log = NamedTuple[]

    for trial in 1:n_trials
        # Randomly choose the "correct" side for this toy trial.
        correct_key = rand([Keyboard.left, Keyboard.right])
        stim_color_value = correct_key == Keyboard.left ? :blue : :red

        # --- Stimulus onset (frame-locked) ---
        onset_time = show_stimulus!(screen, () -> stim_color[] = stim_color_value)

        # --- Wait for response ---
        response = wait_for_response(fig, onset_time; timeout=2.0)

        # --- Clear stimulus (ITI) ---
        show_stimulus!(screen, () -> stim_color[] = :gray)

        if response === nothing
            push!(log, (trial=trial, key=missing, rt=missing, correct=false, timestamp=onset_time))
            println("Trial $trial: no response (timeout)")
        else
            key, rt = response
            correct = key == correct_key
            push!(log, (trial=trial, key=key, rt=rt, correct=correct, timestamp=onset_time))
            println("Trial $trial: key=$key rt=$(round(rt*1000, digits=1))ms correct=$correct")
        end

        sleep(0.5)  # fixed ITI; replace with jittered ITI as needed
    end

    println("\nDone. $(length(log)) trials logged.")
    return log, fig, screen
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    # Set fullscreen=false during development.
    log, fig, screen = run_experiment(; n_trials=5, fullscreen=true)

    # Save results
    open("results_$(Dates.format(now(), "yyyymmdd_HHMMSS")).csv", "w") do io
        println(io, "trial,key,rt,correct,timestamp")
        for row in log
            println(io, "$(row.trial),$(row.key),$(row.rt),$(row.correct),$(row.timestamp)")
        end
    end
end
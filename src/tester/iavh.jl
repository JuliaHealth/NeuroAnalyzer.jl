export iavh

# left and right channel weights of a sound source at (x, y) (normalized head diagram
# coordinates in [-0.5, 0.5]², up = front): power-law attenuation with the distance to each
# simulated ear, the farther ear further attenuated; sources too close to the ears (inside
# the head) are silent
function _avh_weights(x::Real, y::Real)::Tuple{Float64, Float64}
    dist_exponent = 4 # power applied to (1 - distance) for attenuation
    lateral_atten = 0.75 # quieter-ear reduction when sound is lateralized
    ear_l_x = -0.05 # left-ear x position
    ear_r_x = 0.05 # right-ear x position
    min_dist = 0.1 # inside this radius the sound is treated as "in-head"

    dist_l = sqrt((x - ear_l_x)^2 + y^2)
    dist_r = sqrt((x - ear_r_x)^2 + y^2)
    (dist_l < min_dist || dist_r < min_dist) && return (0.0, 0.0)

    d_l = clamp((1.0 - dist_l)^dist_exponent, 0.0, 1.0)
    d_r = clamp((1.0 - dist_r)^dist_exponent, 0.0, 1.0)
    if d_l > d_r
        d_r = round(d_r * lateral_atten; digits = 2)
    elseif d_r > d_l
        d_l = round(d_l * lateral_atten; digits = 2)
    end
    return (round(d_l; digits = 2), round(d_r; digits = 2))
end

# session settings exported by the Save button
function _avh_save(
    file_name::String,
    type::String,
    lang::String,
    gender::String,
    character::String,
    d_l::Real,
    d_r::Real,
    vol::Real,
)::Nothing
    open(file_name, "w") do f
        println(f, "\"AH type\",$type")
        # non-speech types have no language/gender/affect settings
        is_voice = type == "voice"
        println(f, "\"AVH language\",$(is_voice ? lang : "NA")")
        println(f, "\"AVH gender\",$(is_voice ? gender : "NA")")
        println(f, "\"AVH emotional aspect\",$(is_voice ? character : "NA")")
        println(f, "\"distance L\",$d_l")
        println(f, "\"distance R\",$d_r")
        return print(f, "\"volume\",$vol")
    end
    return nothing
end

# window of iavh(); `play(snd, vol, d_l, d_r)` plays a sound (replaced in tests)
function _iavh_build(play::Function)
    vol_step = 0.1 # volume increment / decrement per button press
    vol_min = 0.1 # minimum allowed volume
    vol_max = 1.0 # maximum allowed volume

    # voice file index ranges - 15 files per language×gender set, 5 per valence
    idx_negative = 1:5
    idx_positive = 6:10
    idx_neutral = 11:15

    # =========================================================================
    # sounds
    # =========================================================================
    snd_whisper = wavread(joinpath(res_path, "avh/wav/whisper_2s01fifo.wav"))
    snd_noise = wavread(joinpath(res_path, "avh/wav/noise_2s01fifo.wav"))
    snd_sine = wavread(joinpath(res_path, "avh/wav/sine_8k2s01fifo.wav"))

    # each language × gender set contains 15 WAV files named "<lang>_<g>_NN.wav":
    #   files 01–05 → negative emotional content
    #   files 06–10 → positive emotional content
    #   files 11–15 → neutral emotional content
    _load_voices(prefix::String) = [
        wavread(joinpath(res_path, "avh/wav/$(prefix)_$(lpad(i, 2, '0')).wav"))
        for i = 1:15
    ]
    voices = Dict(
        ("EN", "male")   => _load_voices("en_m"),
        ("EN", "female") => _load_voices("en_w"),
        ("DE", "male")   => _load_voices("de_m"),
        ("DE", "female") => _load_voices("de_w"),
        ("SP", "male")   => _load_voices("sp_m"),
        ("SP", "female") => _load_voices("sp_w"),
        ("PL", "male")   => _load_voices("pl_m"),
        ("PL", "female") => _load_voices("pl_w"),
    )

    # for "voice", a random utterance from the 5-file pool of the language, gender and
    # emotional valence; for other types the fixed non-speech sound
    function select_sound(type::String, lang::String, character::String, gender::String)
        if type == "voice"
            idx_range = character == "negative" ? idx_negative :
                        character == "positive" ? idx_positive : idx_neutral
            return voices[(lang, gender)][rand(idx_range)]
        end
        return type == "whisper" ? snd_whisper : type == "noise" ? snd_noise : snd_sine
    end

    # =========================================================================
    # state
    # =========================================================================
    d_l = Ref(1.0) # left-channel spatial weight (updated by clicks on the head)
    d_r = Ref(1.0) # right-channel spatial weight
    vol = Observable(vol_max) # playback volume
    snd = Ref{Any}(snd_whisper) # current sound
    pending_save = Ref(false)
    status = Observable("Please use headphones for the best results.")
    msg!(s) = (status[] = s; nothing)

    langs = ["EN", "DE", "SP", "PL"]
    types = ["voice", "whisper", "noise", "ringing"]
    genders = ["male", "female"]
    characters = ["negative", "neutral", "positive"]

    # =========================================================================
    # layout: controls | head diagram
    # =========================================================================
    img = rotr90(FileIO.load(joinpath(res_path, "avh/head.png")))
    fig = GLMakie.Figure(; size = (1200, 850))
    ctrl = fig[1, 1] = GLMakie.GridLayout(; valign = :top, tellheight = false)
    GLMakie.colsize!(fig.layout, 1, GLMakie.Fixed(340))

    GLMakie.Label(ctrl[1, 1], "Type"; halign = :right)
    menu_type = GLMakie.Menu(ctrl[1, 2]; options = uppercase.(types), width = 160)
    lb_lang = GLMakie.Label(ctrl[2, 1], "Language"; halign = :right)
    menu_lang = GLMakie.Menu(ctrl[2, 2]; options = langs, width = 160)
    lb_gender = GLMakie.Label(ctrl[3, 1], "Gender"; halign = :right)
    menu_gender = GLMakie.Menu(ctrl[3, 2]; options = uppercase.(genders), width = 160)
    lb_character = GLMakie.Label(ctrl[4, 1], "Emotional aspect"; halign = :right)
    menu_character =
        GLMakie.Menu(ctrl[4, 2]; options = uppercase.(characters), width = 160)

    lb_vol = GLMakie.Label(ctrl[5, 1], lift(v -> "Volume: $v", vol); halign = :right)
    g_vol = ctrl[5, 2] = GLMakie.GridLayout(; halign = :left)
    bt_vol_down = GLMakie.Button(g_vol[1, 1]; label = "−", width = 50)
    bt_vol_up = GLMakie.Button(g_vol[1, 2]; label = "+", width = 50)

    bt_play = GLMakie.Button(ctrl[6, 1:2]; label = "Play", width = 260)
    GLMakie.Label(ctrl[7, 1:2], "Settings file"; halign = :left)
    tb_file = GLMakie.Textbox(
        ctrl[8, 1:2];
        width = 260,
        placeholder = "file path (.csv)",
    )
    bt_save = GLMakie.Button(ctrl[9, 1:2]; label = "Save", width = 260)
    bt_close = GLMakie.Button(ctrl[10, 1:2]; label = "Close", width = 260)
    GLMakie.Label(
        ctrl[11, 1:2],
        status;
        halign = :left,
        word_wrap = true,
        width = 320,
        color = :gray25,
    )
    GLMakie.rowgap!(ctrl, 12)

    ax = GLMakie.Axis(
        fig[1, 2];
        aspect = GLMakie.DataAspect(),
        title = "Click to place the sound source",
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.hidedecorations!(ax)
    GLMakie.hidespines!(ax)
    GLMakie.image!(ax, img)
    GLMakie.xlims!(ax, 0, size(img, 1))
    GLMakie.ylims!(ax, 0, size(img, 2))
    source = Observable(GLMakie.Point2f[])
    GLMakie.scatter!(
        ax,
        source;
        marker = :circle,
        markersize = 20,
        color = :transparent,
        strokecolor = :red,
        strokewidth = 2,
    )

    # =========================================================================
    # behavior
    # =========================================================================
    selections() = (
        types[menu_type.i_selected[]],
        langs[menu_lang.i_selected[]],
        characters[menu_character.i_selected[]],
        genders[menu_gender.i_selected[]],
    )
    new_sound!() = (snd[] = select_sound(selections()...); nothing)
    play_current() = play(snd[], vol[], d_l[], d_r[])

    # speech settings are used only for voices: their labels are grayed out otherwise
    on(menu_type.i_selected) do i
        c = types[i] == "voice" ? :black : :gray70
        for lb in (lb_lang, lb_gender, lb_character)
            lb.color[] = GLMakie.to_color(c)
        end
        return new_sound!()
    end
    for m in (menu_lang, menu_gender, menu_character)
        on(_ -> new_sound!(), m.i_selected)
    end

    # volume changes preview the current sound without choosing another utterance
    on(bt_vol_up.clicks) do _
        vol[] = round(min(vol[] + vol_step, vol_max); digits = 1)
        return play_current()
    end
    on(bt_vol_down.clicks) do _
        vol[] = round(max(vol[] - vol_step, vol_min); digits = 1)
        return play_current()
    end

    # each Play chooses another utterance, as AVH are unpredictable
    on(bt_play.clicks) do _
        new_sound!()
        return play_current()
    end

    # click on the head diagram: place the sound source and play
    function place!(x::Real, y::Real)
        source[] = [GLMakie.Point2f(x, y)]
        d_l[], d_r[] = _avh_weights(x / size(img, 1) - 0.5, y / size(img, 2) - 0.5)
        new_sound!()
        play_current()
        return msg!("Sound source: left $(d_l[]), right $(d_r[]).")
    end
    on(GLMakie.events(ax.scene).mousebutton; priority = 2) do ev
        if ev.button == GLMakie.Mouse.left &&
           ev.action == GLMakie.Mouse.press &&
           GLMakie.is_mouseinside(ax.scene)
            p = GLMakie.mouseposition(ax)
            place!(p[1], p[2])
            return GLMakie.Consume(true)
        end
        return GLMakie.Consume(false)
    end

    # save the session settings; an existing file is overwritten after a second click
    on(_ -> (pending_save[] = false), tb_file.displayed_string)
    on(bt_save.clicks) do _
        f = String(strip(something(tb_file.displayed_string[], "")))
        isempty(f) && return msg!("Enter the file path first.")
        if isfile(f) && !pending_save[]
            pending_save[] = true
            return msg!("$f exists: click Save again to overwrite it.")
        end
        pending_save[] = false
        try
            type, lang, character, gender = selections()
            _avh_save(f, type, lang, gender, character, d_l[], d_r[], vol[])
        catch err
            return msg!("File cannot be saved: $(sprint(showerror, err))")
        end
        return msg!("Settings saved to $f.")
    end

    ctl = Dict{Symbol, Any}(
        :type => menu_type,
        :lang => menu_lang,
        :gender => menu_gender,
        :character => menu_character,
        :vol_up => bt_vol_up,
        :vol_down => bt_vol_down,
        :play => bt_play,
        :file => tb_file,
        :save => bt_save,
        :close => bt_close,
        :status => status,
    )

    return (; fig, ctl, vol, d_l, d_r, snd, source, place!, labels = (lb_lang, lb_gender, lb_character))
end

# play `snd` (as returned by wavread()) with the volume and channel weights; playback runs in
# another task, as wavplay() blocks until the sound ends
function _avh_play(snd, vol::Real, d_l::Real, d_r::Real)::Nothing
    channel_gain = 0.25 # master scalar applied to both channel amplitudes
    y = snd[1] .* [vol * d_l * channel_gain vol * d_r * channel_gain]
    Threads.@spawn try
        wavplay(y, snd[2])
    catch err
        @error "Playback failed" exception = (err, catch_backtrace())
    end
    return nothing
end

"""
    iavh()

**Interactive Auditory (Verbal) Hallucination (iAVH) Simulator** - a clinical and research tool for recreating auditory hallucination experiences in a controlled, configurable environment.

Practitioners can combine spatial positioning, volume, sound type, and (for speech) linguistic and emotional content to explore what a patient might perceive. At the end of a session the complete parameter set can be exported to a CSV file for further analysis or record keeping.

Sound types:

| Type     | Description                                               |
|----------|-----------------------------------------------------------|
| Voice    | Auditory Verbal Hallucinations (AVH) - intelligible speech|
| Whisper  | Unintelligible whispered background voice                 |
| Noise    | Broadband background noise                                |
| Ringing  | Pure 8 000 Hz sine tone (tinnitus simulation)             |

Speech settings (active only when *Voice* is selected)

- **Language** - English (EN), German (DE), Spanish (SP), Polish (PL)
- **Gender** - Male / Female voice
- **Emotional aspect** - Negative / Neutral / Positive content

Each language × gender × valence combination draws randomly from a pool of five pre-recorded utterances on every playback event, reflecting the unpredictable nature of AVH.

Spatial positioning

Click anywhere on the head diagram to place the sound source in space. The left and right channel amplitudes are derived from the Euclidean distance between the clicked point and each simulated ear position using a power-law attenuation model `amplitude = (1 − distance)⁴`, with additional inter-aural asymmetry applied when the source is lateralized (the farther ear receives 75 % of the nearer ear's amplitude). Sounds originating inside the minimum distance threshold receive zero spatialization (equal silence in both ears at that distance).

Volume control

Use the **+** / **−** buttons to adjust playback volume in increments of 0.1 (range: 0.1 – 1.0). Volume changes preview the *current* sound at the new level without re-randomizing the selected utterance.

Export (Save button)

Enter the file path and press **Save** to write a CSV containing: AH type, AVH language, AVH gender, AVH emotional aspect, left channel distance weight, right channel distance weight, and volume. Non-speech types record `NA` for the speech-only fields.

**Headphones are required**. Accurate binaural spatialization requires headphones. The window shows a reminder on startup.

# Arguments

None.

# Returns

- `Nothing`
"""
function iavh()::Nothing
    w = _iavh_build(_avh_play)
    screen = GLMakie.Screen(; title = "NeuroTester: iavh()")
    on(_ -> close(screen), w.ctl[:close].clicks)
    display(screen, w.fig)
    wait(screen)
    return nothing
end

export iavh

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

Writes a CSV containing: AH type, AVH language, AVH gender, AVH emotional aspect, left channel distance weight, right channel distance weight, and volume. Non-speech types record `NA` for the speech-only fields.

!!! note "Headphones required"

Accurate binaural spatialization requires headphones. The application shows a reminder dialog on startup.

# Arguments

None.

# Returns

- `Nothing`
"""
function iavh()::Nothing

    # =========================================================================
    # configuration - named locals used in place of magic numbers throughout
    # =========================================================================
    vol_step      = 0.1   # volume increment / decrement per button press
    vol_min       = 0.1   # minimum allowed volume
    vol_max       = 1.0   # maximum allowed volume
    dist_exponent = 4     # power applied to (1 - distance) for attenuation
    channel_gain  = 0.25  # master scalar applied to both channel amplitudes
    lateral_atten = 0.75  # quieter-ear reduction when sound is lateralized
    ear_l_x       = -0.05 # left-ear x position in normalized [-0.5, 0.5] space
    ear_r_x       = 0.05 # right-ear x position in normalized [-0.5, 0.5] space
    min_dist      = 0.1   # inside this radius the sound is treated as "in-head"
    canvas_size   = 800   # canvas width and height in pixels

    # voice file index ranges - 15 files per language×gender set, 5 per valence
    idx_negative = 1:5
    idx_positive = 6:10
    idx_neutral  = 11:15

    # =========================================================================
    # state - mutable values updated by UI interactions
    # =========================================================================
    d_l = 1.0     # left-channel spatial weight  (updated by canvas clicks)
    d_r = 1.0     # right-channel spatial weight (updated by canvas clicks)
    vol = vol_max # current playback volume

    # =========================================================================
    # load non-speech sound assets
    # =========================================================================
    snd_whisper = wavread(joinpath(res_path, "avh/wav/whisper_2s01fifo.wav"))
    snd_noise   = wavread(joinpath(res_path, "avh/wav/noise_2s01fifo.wav"))
    snd_sine    = wavread(joinpath(res_path, "avh/wav/sine_8k2s01fifo.wav"))

    # =========================================================================
    # load speech samples
    # each language × gender set contains 15 WAV files named "<lang>_<g>_NN.wav":
    #   files 01–05 → negative emotional content
    #   files 06–10 → positive emotional content
    #   files 11–15 → neutral emotional content
    # =========================================================================
    function _load_voices(prefix::String)
        return [
            wavread(joinpath(res_path, "avh/wav/$(prefix)_$(lpad(i, 2, '0')).wav"))
            for i in 1:15
        ]
    end

    # keyed by (language_code, gender_string) for O(1) lookup in signal handlers
    voices = Dict{
        Tuple{String, String},
        Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}},
    }(
        ("EN", "male")   => _load_voices("en_m"),
        ("EN", "female") => _load_voices("en_w"),
        ("DE", "male")   => _load_voices("de_m"),
        ("DE", "female") => _load_voices("de_w"),
        ("SP", "male")   => _load_voices("sp_m"),
        ("SP", "female") => _load_voices("sp_w"),
        ("PL", "male")   => _load_voices("pl_m"),
        ("PL", "female") => _load_voices("pl_w"),
    )

    # currently active raw sound (updated whenever settings change or Play is pressed)
    snd = deepcopy(snd_whisper)

    img = read_from_png(joinpath(res_path, "avh/head.png"))

    # =========================================================================
    # helper: pick a sound clip from the appropriate pool
    #
    # for "voice" type, selects randomly from the 5-file pool that matches
    # the given language, gender, and emotional valence
    # for all other types the fixed non-speech asset is returned
    # =========================================================================
    function _select_sound(type::String, lang::String, character::String, gender::String)
        if type == "voice"
            idx_range = if character == "negative"
                idx_negative
            elseif character == "positive"
                idx_positive
            else  # neutral
                idx_neutral
            end
            return voices[(lang, gender)][rand(idx_range)]
        elseif type == "whisper"
            return snd_whisper
        elseif type == "noise"
            return snd_noise
        else  # ringing
            return snd_sine
        end
    end

    # =========================================================================
    # helper: apply volume and spatial weights to a sound, then play it
    #
    # creates a temporary copy so the original raw samples are never mutated
    # channel amplitudes = raw_sample × vol × (channel_weight × channel_gain)
    # =========================================================================
    function _play(snd_raw)
        snd_out = deepcopy(snd_raw)
        snd_out[1][:, 1] = snd_raw[1][:, 1] .* (vol * (d_l * channel_gain))
        snd_out[1][:, 2] = snd_raw[1][:, 2] .* (vol * (d_r * channel_gain))
        return wavplay(snd_out[1], snd_out[2])
    end

    # =========================================================================
    # helper: read current selections from the four combo boxes
    # returns (type, lang, character, gender) as plain lowercase strings
    # =========================================================================
    function _get_ui_state(combo_type, combo_lang, combo_character, combo_gender,
        types, langs, characters, genders)
        type      = types[Int64(combo_type.active) + 1]
        lang      = langs[Int64(combo_lang.active) + 1]
        character = characters[Int64(combo_character.active) + 1]
        gender    = genders[Int64(combo_gender.active) + 1]
        return type, lang, character, gender
    end

    # =========================================================================
    # GTK application
    # =========================================================================
    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroTester: iavh()")
        Gtk4.default_size(win, 1100, 820)

        # canvas that renders the head diagram and the clicked sound-source marker
        can                = GtkCanvas()
        can.content_width  = canvas_size
        can.content_height = canvas_size

        # outer grid: options panel (column 1) | canvas (column 2)
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 10
        g.row_spacing = 10
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        # inner grid for all control widgets
        g_opts = GtkGrid()
        g_opts.column_homogeneous = false
        g_opts.column_spacing = 10
        g_opts.row_spacing = 10
        g_opts.margin_start = 5
        g_opts.margin_end = 5
        g_opts.margin_top = 5
        g_opts.margin_bottom = 5

        # --- language selector ---
        lab_lang = GtkLabel("Language");
        lab_lang.halign = 2
        langs = ["EN", "DE", "SP", "PL"]
        combo_lang = GtkComboBoxText()
        foreach(l -> push!(combo_lang, l), langs)
        combo_lang.active       = 0
        combo_lang.tooltip_text = "Language of AVH speech samples"

        # --- sound type selector ---
        lab_type = GtkLabel("Type");
        lab_type.halign = 2
        types = ["voice", "whisper", "noise", "ringing"]
        combo_type = GtkComboBoxText()
        foreach(t -> push!(combo_type, uppercase(t)), types)
        combo_type.active       = 0
        combo_type.tooltip_text = "Hallucination sound type"

        # --- voice gender selector (speech only) ---
        lab_gender = GtkLabel("Gender");
        lab_gender.halign = 2
        genders = ["male", "female"]
        combo_gender = GtkComboBoxText()
        foreach(g -> push!(combo_gender, uppercase(g)), genders)
        combo_gender.active       = 0
        combo_gender.tooltip_text = "Voice gender (speech only)"

        # --- emotional aspect selector (speech only) ---
        lab_character = GtkLabel("Emotional aspect");
        lab_character.halign = 2
        characters = ["negative", "neutral", "positive"]
        combo_character = GtkComboBoxText()
        foreach(c -> push!(combo_character, uppercase(c)), characters)
        combo_character.active       = 0
        combo_character.tooltip_text = "Emotional valence of AVH content (speech only)"

        # --- volume controls ---
        lab_vol_up = GtkLabel("Volume");
        lab_vol_up.halign = 2
        lab_vol_down = GtkLabel("Volume");
        lab_vol_down.halign = 2
        bt_vol_up = GtkButton("+")
        bt_vol_down = GtkButton("-")
        # initialise sensitivity: at max vol the up button is disabled
        bt_vol_up.sensitive      = (vol < vol_max) ? 1 : 0
        bt_vol_down.sensitive    = (vol > vol_min) ? 1 : 0
        bt_vol_up.tooltip_text   = "Increase volume (step $(vol_step))"
        bt_vol_down.tooltip_text = "Decrease volume (step $(vol_step))"

        # --- action buttons ---
        bt_play               = GtkButton("Play");
        bt_play.tooltip_text  = "Play the sound with current settings"
        bt_save               = GtkButton("Save");
        bt_save.tooltip_text  = "Export session settings to CSV"
        bt_close              = GtkButton("Close");
        bt_close.tooltip_text = "Close this window"

        # --- populate options grid ---
        g_opts[1, 1]    = lab_lang;
        g_opts[2, 1]    = combo_lang
        g_opts[1, 2]    = lab_type;
        g_opts[2, 2]    = combo_type
        g_opts[1, 3]    = lab_gender;
        g_opts[2, 3]    = combo_gender
        g_opts[1, 4]    = lab_character;
        g_opts[2, 4]    = combo_character
        g_opts[1, 5]    = lab_vol_up;
        g_opts[2, 5]    = bt_vol_up
        g_opts[1, 6]    = lab_vol_down;
        g_opts[2, 6]    = bt_vol_down
        g_opts[1:2, 7]  = GtkLabel("")
        g_opts[1:2, 8]  = bt_play
        g_opts[1:2, 9]  = GtkLabel("")
        g_opts[1:2, 10] = bt_save
        g_opts[1:2, 11] = GtkLabel("")
        g_opts[1:2, 12] = bt_close

        vbox = GtkBox(:v)
        push!(vbox, g_opts)

        g[1, 1] = vbox
        g[2, 1] = can
        push!(win, g)
        Gtk4.show(win)

        info_dialog("Please use headphones for the best results.", win) do
            return nothing
        end

        # --- canvas: initial draw (head diagram, no marker) ---
        @guarded draw(can) do widget
            ctx = getgc(widget)
            h, w = Cairo.height(can), Cairo.width(can)
            Cairo.rectangle(ctx, 0, 0, w, h)
            Cairo.set_source_rgb(ctx, 1, 1, 1)
            Cairo.fill(ctx)
            Cairo.set_source_surface(ctx, img, 1, 1)
            return Cairo.paint(ctx)
        end

        # -----------------------------------------------------------------------
        # shared handler for any speech-setting combo-box change.
        # updates the active sound without triggering playback; the user must
        # press Play or click the canvas to hear the result
        # -----------------------------------------------------------------------
        function _on_speech_setting_changed()
            type, lang, character, gender = _get_ui_state(
                combo_type, combo_lang, combo_character, combo_gender,
                types, langs, characters, genders,
            )
            return snd = _select_sound(type, lang, character, gender)
        end

        # --- signal: sound type changed ---
        # enable/disable speech-only controls depending on whether "voice" is active
        signal_connect(combo_type, "changed") do widget
            type                      = types[Int64(combo_type.active) + 1]
            is_voice                  = (type == "voice")
            combo_lang.sensitive      = is_voice ? 1 : 0
            combo_character.sensitive = is_voice ? 1 : 0
            combo_gender.sensitive    = is_voice ? 1 : 0
            return _on_speech_setting_changed()
        end

        # --- signal: language changed (was missing in original) ---
        signal_connect(combo_lang, "changed") do widget
            return _on_speech_setting_changed()
        end

        # --- signal: voice gender changed ---
        signal_connect(combo_gender, "changed") do widget
            return _on_speech_setting_changed()
        end

        # --- signal: emotional aspect changed ---
        signal_connect(combo_character, "changed") do widget
            return _on_speech_setting_changed()
        end

        # --- signal: volume up ---
        # previews the *current* sound at the new volume without re-randomizing
        signal_connect(bt_vol_up, "clicked") do widget
            vol                   = round(min(vol + vol_step, vol_max); digits = 1)
            bt_vol_up.sensitive   = (vol < vol_max) ? 1 : 0
            bt_vol_down.sensitive = 1  # definitely above minimum after an increase
            return _play(snd)
        end

        # --- signal: volume down ---
        signal_connect(bt_vol_down, "clicked") do widget
            vol                   = round(max(vol - vol_step, vol_min); digits = 1)
            bt_vol_down.sensitive = (vol > vol_min) ? 1 : 0
            bt_vol_up.sensitive   = 1  # definitely below maximum after a decrease
            return _play(snd)
        end

        # --- signal: play button ---
        # re-randomizes the utterance on each press to simulate unpredictable AVH
        signal_connect(bt_play, "clicked") do widget
            type, lang, character, gender = _get_ui_state(
                combo_type, combo_lang, combo_character, combo_gender,
                types, langs, characters, genders,
            )
            snd = _select_sound(type, lang, character, gender)
            return _play(snd)
        end

        # -----------------------------------------------------------------------
        # canvas left-click: reposition the sound source and play.
        #
        # coordinate mapping:
        #   canvas pixel (px, py) → normalized (x, y) in [-0.5, 0.5]²
        #   x = px / canvas_size - 0.5
        #   y = 0.5 - py / canvas_size   (y-axis flipped: up = positive)
        #
        # ear positions (normalized): L = (ear_l_x, 0), R = (ear_r_x, 0)
        #
        # attenuation model:
        #   amplitude = clamp((1 - dist)^dist_exponent, 0, 1)
        # where dist is the Euclidean distance from the sound source to the ear
        # sources within min_dist of either ear are treated as "in-head" and
        # both channels are set to zero
        # -----------------------------------------------------------------------
        function _lmb_click(_, _, x_px, y_px)
            return Threads.@spawn begin
                # redraw the head diagram with a red circle at the clicked position
                @guarded draw(can) do widget
                    ctx = getgc(widget)
                    h, w = Cairo.height(can), Cairo.width(can)
                    Cairo.rectangle(ctx, 0, 0, w, h)
                    Cairo.set_source_rgb(ctx, 1, 1, 1)
                    Cairo.fill(ctx)
                    Cairo.set_source_surface(ctx, img, 1, 1)
                    Cairo.paint(ctx)
                    Gtk4.arc(ctx, x_px, y_px, 10, 0, 2π)
                    Gtk4.set_source_rgb(ctx, 1, 0, 0)
                    Gtk4.stroke(ctx)
                    return Gtk4.reveal(widget)
                end

                # map pixel coordinates to normalized [-0.5, 0.5] space
                x_norm = round(x_px / canvas_size - 0.5, digits = 2)
                y_norm = round(0.5 - y_px / canvas_size, digits = 2)

                # euclidean distances from the sound source to each simulated ear
                dist_l = sqrt((x_norm - ear_l_x)^2 + y_norm^2)
                dist_r = sqrt((x_norm - ear_r_x)^2 + y_norm^2)

                if dist_l < min_dist || dist_r < min_dist
                    # source is inside (or too close to) the head - silence both channels
                    d_l = 0.0
                    d_r = 0.0
                else
                    # power-law attenuation, clamped to [0, 1]
                    d_l = clamp((1.0 - dist_l)^dist_exponent, 0.0, 1.0)
                    d_r = clamp((1.0 - dist_r)^dist_exponent, 0.0, 1.0)

                    # attenuate the farther ear to enhance perceived lateralization
                    if d_l > d_r
                        d_r = round(d_r * lateral_atten; digits = 2)
                    elseif d_r > d_l
                        d_l = round(d_l * lateral_atten; digits = 2)
                    end

                    d_l = round(d_l; digits = 2)
                    d_r = round(d_r; digits = 2)
                end

                # re-randomize and play with the updated spatial weights
                type, lang, character, gender = _get_ui_state(
                    combo_type, combo_lang, combo_character, combo_gender,
                    types, langs, characters, genders,
                )
                snd = _select_sound(type, lang, character, gender)
                _play(snd)
            end
        end

        ggc_l = GtkGestureClick()
        ggc_l.button = 1
        push!(can, ggc_l)
        signal_connect(_lmb_click, ggc_l, "pressed")

        # --- signal: save button - export session parameters to CSV ---
        signal_connect(bt_save, "clicked") do widget
            save_dialog("Save as...", win, ["*.csv", ".txt"]) do file_name
                isempty(file_name) && return

                type, lang, character, gender = _get_ui_state(
                    combo_type, combo_lang, combo_character, combo_gender,
                    types, langs, characters, genders,
                )

                try
                    open(file_name, "w") do f
                        println(f, "\"AH type\",$type")
                        if type == "voice"
                            println(f, "\"AVH language\",$lang")
                            println(f, "\"AVH gender\",$gender")
                            println(f, "\"AVH emotional aspect\",$character")
                        else
                            # non-speech types have no language/gender/affect settings
                            println(f, "\"AVH language\",NA")
                            println(f, "\"AVH gender\",NA")
                            println(f, "\"AVH emotional aspect\",NA")
                        end
                        println(f, "\"distance L\",$d_l")
                        println(f, "\"distance R\",$d_r")
                        return print(f, "\"volume\",$vol")
                    end
                catch
                    warn_dialog(_nill, "File cannot be saved!", win)
                end
            end
        end

        # --- signal: Close button ---
        return signal_connect(bt_close, "clicked") do widget
            return close(win)
        end
    end # _activate

    app = GtkApplication("org.neuroanalyzer.iavh")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing
end

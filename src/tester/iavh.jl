export iavh

"""
    iavh()

Interactive tool for subjective recreation of auditory hallucinations. Four languages are available (English, German, Spanish and Polish). Four types of auditory hallucinations are available: speech (Auditory Verbal Hallucinations), whispers, noise and ringing (8000 Hz).

Patient may locate the sounds in space to recreate how distant or close they are. Also, volume of the auditory experiences may be modified. In case of AVH, emotional content (negative, neutral or positive) and voice gender (male or female) may be chosen.

When done, patient's settings may be easily exported to a CSV file for further analysis.

# Arguments

Nothing

# Returns

Nothing
"""
function iavh()::Nothing

    info_dialog("Please use headphones for the best results.")

    d_l = 1
    d_r = 1
    vol = 1.0
    snd_whisper = wavread(joinpath(res_path, "avh/wav/whisper_2s01fifo.wav"))
    snd_noise = wavread(joinpath(res_path, "avh/wav/noise_2s01fifo.wav"))
    snd_sine = wavread(joinpath(res_path, "avh/wav/sine_8k2s01fifo.wav"))
    voices_en_m = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_en_m, wavread(joinpath(res_path, "avh/wav/en_m_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_en_w = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_en_w, wavread(joinpath(res_path, "avh/wav/en_w_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_de_m = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_de_m, wavread(joinpath(res_path, "avh/wav/de_m_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_de_w = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_de_w, wavread(joinpath(res_path, "avh/wav/de_w_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_sp_m = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_sp_m, wavread(joinpath(res_path, "avh/wav/sp_m_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_sp_w = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_sp_w, wavread(joinpath(res_path, "avh/wav/sp_w_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_pl_m = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_pl_m, wavread(joinpath(res_path, "avh/wav/pl_m_$(lpad(string(idx), 2, "0")).wav")))
    end
    voices_pl_w = Vector{Tuple{Matrix{Float64}, Float32, UInt16, Vector{WAVChunk}}}()
    for idx in 1:15
        push!(voices_pl_w, wavread(joinpath(res_path, "avh/wav/pl_w_$(lpad(string(idx), 2, "0")).wav")))
    end

    snd = deepcopy(snd_whisper)
    snd_tmp = deepcopy(snd)

    img = read_from_png(joinpath(res_path, "avh/head.png"))

    win = GtkWindow("NeuroTester: iavh", 1100, 820)
    set_gtk_property!(win, :border_width, 10)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    can = GtkCanvas(800, 800)

    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_spacing, 10)
    set_gtk_property!(g_opts, :row_spacing, 10)

    lab_lang = GtkLabel("Language")
    set_gtk_property!(lab_lang, :halign, 2)
    langs = ["EN", "DE", "SP", "PL"]
    combo_lang = GtkComboBoxText()
    for idx in langs
        push!(combo_lang, idx)
    end
    set_gtk_property!(combo_lang, :active, 0)
    set_gtk_property!(combo_lang, :tooltip_text, "AVH language")

    lab_character = GtkLabel("Emotional aspect")
    set_gtk_property!(lab_character, :halign, 2)
    characters = ["negative", "neutral", "positive"]
    combo_character = GtkComboBoxText()
    for idx in uppercase.(String.(characters))
        push!(combo_character, idx)
    end
    set_gtk_property!(combo_character, :active, 0)
    set_gtk_property!(combo_character, :tooltip_text, "AVH content emotional aspect")

    lab_type = GtkLabel("Type")
    set_gtk_property!(lab_type, :halign, 2)
    types = ["voice", "whisper", "noise", "ringing"]
    combo_type = GtkComboBoxText()
    for idx in uppercase.(String.(types))
        push!(combo_type, idx)
    end
    set_gtk_property!(combo_type, :active, 0)
    set_gtk_property!(combo_type, :tooltip_text, "VH sounds type")

    lab_gender = GtkLabel("Gender")
    set_gtk_property!(lab_gender, :halign, 2)
    genders = ["male", "female"]
    combo_gender = GtkComboBoxText()
    for idx in uppercase.(String.(genders))
        push!(combo_gender, idx)
    end
    set_gtk_property!(combo_gender, :active, 0)
    set_gtk_property!(combo_gender, :tooltip_text, "AVH voice gender characteristic")

    lab_vol_up = GtkLabel("Volume")
    set_gtk_property!(lab_vol_up, :halign, 2)
    lab_vol_down = GtkLabel("Volume")
    set_gtk_property!(lab_vol_down, :halign, 2)
    bt_vol_up = GtkButton("+")
    vol == 1 && set_gtk_property!(bt_vol_up, :sensitive, 0)
    set_gtk_property!(bt_vol_up, :tooltip_text, "Increase volume")
    bt_vol_down = GtkButton("-")
    set_gtk_property!(bt_vol_down, :tooltip_text, "Decrease volume")

    bt_play = GtkButton("Play")
    set_gtk_property!(bt_play, :tooltip_text, "Play the sound with current settings")

    bt_save = GtkButton("Save")
    set_gtk_property!(bt_save, :tooltip_text, "Save results")

    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    g_opts[1, 1] = lab_lang
    g_opts[2, 1] = combo_lang
    g_opts[1, 2] = lab_type
    g_opts[2, 2] = combo_type
    g_opts[1, 3] = lab_gender
    g_opts[2, 3] = combo_gender
    g_opts[1, 4] = lab_character
    g_opts[2, 4] = combo_character
    g_opts[1, 5] = lab_vol_up
    g_opts[2, 5] = bt_vol_up
    g_opts[1, 6] = lab_vol_down
    g_opts[2, 6] = bt_vol_down
    g_opts[1:2, 7] = ""
    g_opts[1:2, 8] = bt_play
    g_opts[1:2, 9] = ""
    g_opts[1:2, 10] = bt_save
    g_opts[1:2, 11] = ""
    g_opts[1:2, 12] = bt_close
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2, 1] = can
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(widget)
        h = Cairo.height(can)
        w = Cairo.width(can)
        Cairo.rectangle(ctx, 0, 0, w, h)
        Cairo.set_source_rgb(ctx, 1, 1, 1)
        Cairo.fill(ctx)
        Cairo.set_source_surface(ctx, img, 1, 1)
        Cairo.paint(ctx)
    end

    signal_connect(combo_type, "changed") do widget
        type = types[get_gtk_property(combo_type, :active, Int64) + 1]
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if type == "voice"
            if lang == "EN"
                if character == "negative"
                    if gender == "male"
                        snd = voices_en_m[rand(1:5, 1)[]]
                    else
                        snd = voices_en_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_en_m[rand(6:10, 1)[]]
                    else
                        snd = voices_en_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_en_m[rand(11:15, 1)[]]
                    else
                        snd = voices_en_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "DE"
                if character == "negative"
                    if gender == "male"
                        snd = voices_de_m[rand(1:5, 1)[]]
                    else
                        snd = voices_de_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_de_m[rand(6:10, 1)[]]
                    else
                        snd = voices_de_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_de_m[rand(11:15, 1)[]]
                    else
                        snd = voices_de_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "SP"
                if character == "negative"
                    if gender == "male"
                        snd = voices_sp_m[rand(1:5, 1)[]]
                    else
                        snd = voices_sp_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_sp_m[rand(6:10, 1)[]]
                    else
                        snd = voices_sp_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_sp_m[rand(11:15, 1)[]]
                    else
                        snd = voices_sp_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "PL"
                if character == "negative"
                    if gender == "male"
                        snd = voices_pl_m[rand(1:5, 1)[]]
                    else
                        snd = voices_pl_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_pl_m[rand(6:10, 1)[]]
                    else
                        snd = voices_pl_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_pl_m[rand(11:15, 1)[]]
                    else
                        snd = voices_pl_w[rand(11:15, 1)[]]
                    end
                end
            end
            set_gtk_property!(combo_lang, :sensitive, 1)
            set_gtk_property!(combo_character, :sensitive, 1)
            set_gtk_property!(combo_gender, :sensitive, 1)
        elseif type == "whisper"
            snd = snd_whisper
            set_gtk_property!(combo_lang, :sensitive, 0)
            set_gtk_property!(combo_character, :sensitive, 0)
            set_gtk_property!(combo_gender, :sensitive, 0)
        elseif type == "noise"
            snd = snd_noise
            set_gtk_property!(combo_lang, :sensitive, 0)
            set_gtk_property!(combo_character, :sensitive, 0)
            set_gtk_property!(combo_gender, :sensitive, 0)
        elseif type == "ringing"
            snd = snd_sine
            set_gtk_property!(combo_lang, :sensitive, 0)
            set_gtk_property!(combo_character, :sensitive, 0)
            set_gtk_property!(combo_gender, :sensitive, 0)
        end
        snd_tmp = deepcopy(snd)
    end

    signal_connect(combo_gender, "changed") do widget
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if lang == "EN"
            if character == "negative"
                if gender == "male"
                    snd = voices_en_m[rand(1:5, 1)[]]
                else
                    snd = voices_en_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_en_m[rand(6:10, 1)[]]
                else
                    snd = voices_en_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_en_m[rand(11:15, 1)[]]
                else
                    snd = voices_en_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "DE"
            if character == "negative"
                if gender == "male"
                    snd = voices_de_m[rand(1:5, 1)[]]
                else
                    snd = voices_de_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_de_m[rand(6:10, 1)[]]
                else
                    snd = voices_de_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_de_m[rand(11:15, 1)[]]
                else
                    snd = voices_de_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "SP"
            if character == "negative"
                if gender == "male"
                    snd = voices_sp_m[rand(1:5, 1)[]]
                else
                    snd = voices_sp_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_sp_m[rand(6:10, 1)[]]
                else
                    snd = voices_sp_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_sp_m[rand(11:15, 1)[]]
                else
                    snd = voices_sp_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "PL"
            if character == "negative"
                if gender == "male"
                    snd = voices_pl_m[rand(1:5, 1)[]]
                else
                    snd = voices_pl_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_pl_m[rand(6:10, 1)[]]
                else
                    snd = voices_pl_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_pl_m[rand(11:15, 1)[]]
                else
                    snd = voices_pl_w[rand(11:15, 1)[]]
                end
            end
        end
        snd_tmp = deepcopy(snd)
    end

    signal_connect(combo_character, "changed") do widget
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if lang == "EN"
            if character == "negative"
                if gender == "male"
                    snd = voices_en_m[rand(1:5, 1)[]]
                else
                    snd = voices_en_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_en_m[rand(6:10, 1)[]]
                else
                    snd = voices_en_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_en_m[rand(11:15, 1)[]]
                else
                    snd = voices_en_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "DE"
            if character == "negative"
                if gender == "male"
                    snd = voices_de_m[rand(1:5, 1)[]]
                else
                    snd = voices_de_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_de_m[rand(6:10, 1)[]]
                else
                    snd = voices_de_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_de_m[rand(11:15, 1)[]]
                else
                    snd = voices_de_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "SP"
            if character == "negative"
                if gender == "male"
                    snd = voices_sp_m[rand(1:5, 1)[]]
                else
                    snd = voices_sp_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_sp_m[rand(6:10, 1)[]]
                else
                    snd = voices_sp_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_sp_m[rand(11:15, 1)[]]
                else
                    snd = voices_sp_w[rand(11:15, 1)[]]
                end
            end
        elseif lang == "PL"
            if character == "negative"
                if gender == "male"
                    snd = voices_pl_m[rand(1:5, 1)[]]
                else
                    snd = voices_pl_w[rand(1:5, 1)[]]
                end
            elseif character == "positive"
                if gender == "male"
                    snd = voices_pl_m[rand(6:10, 1)[]]
                else
                    snd = voices_pl_w[rand(6:10, 1)[]]
                end
            else
                if gender == "male"
                    snd = voices_pl_m[rand(11:15, 1)[]]
                else
                    snd = voices_pl_w[rand(11:15, 1)[]]
                end
            end
        end
        snd_tmp = deepcopy(snd)
    end

    signal_connect(bt_vol_up, "clicked") do widget
        vol < 1.0 && (vol += 0.1)
        vol > 0.1 && set_gtk_property!(bt_vol_down, :sensitive, 1)
        vol == 1.0 && set_gtk_property!(bt_vol_up, :sensitive, 0)
        vol = round(vol, digits=1)

        type = types[get_gtk_property(combo_type, :active, Int64) + 1]
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if type == "voice"
            if lang == "EN"
                if character == "negative"
                    if gender == "male"
                        snd = voices_en_m[rand(1:5, 1)[]]
                    else
                        snd = voices_en_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_en_m[rand(6:10, 1)[]]
                    else
                        snd = voices_en_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_en_m[rand(11:15, 1)[]]
                    else
                        snd = voices_en_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "DE"
                if character == "negative"
                    if gender == "male"
                        snd = voices_de_m[rand(1:5, 1)[]]
                    else
                        snd = voices_de_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_de_m[rand(6:10, 1)[]]
                    else
                        snd = voices_de_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_de_m[rand(11:15, 1)[]]
                    else
                        snd = voices_de_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "SP"
                if character == "negative"
                    if gender == "male"
                        snd = voices_sp_m[rand(1:5, 1)[]]
                    else
                        snd = voices_sp_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_sp_m[rand(6:10, 1)[]]
                    else
                        snd = voices_sp_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_sp_m[rand(11:15, 1)[]]
                    else
                        snd = voices_sp_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "PL"
                if character == "negative"
                    if gender == "male"
                        snd = voices_pl_m[rand(1:5, 1)[]]
                    else
                        snd = voices_pl_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_pl_m[rand(6:10, 1)[]]
                    else
                        snd = voices_pl_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_pl_m[rand(11:15, 1)[]]
                    else
                        snd = voices_pl_w[rand(11:15, 1)[]]
                    end
                end
            end
        elseif type == "whisper"
            snd = snd_whisper
        elseif type == "noise"
            snd = snd_noise
        elseif type == "ringing"
            snd = snd_sine
        end
        snd_tmp = deepcopy(snd)

        snd_tmp[1][:, 1] = snd[1][:, 1] .* (vol * (d_l * 0.25))
        snd_tmp[1][:, 2] = snd[1][:, 2] .* (vol * (d_r * 0.25))
        wavplay(snd_tmp[1], snd_tmp[2])
    end

    signal_connect(bt_vol_down, "clicked") do widget
        vol > 0.1 && (vol -= 0.1)
        vol == 0.1 && set_gtk_property!(bt_vol_down, :sensitive, 0)
        vol < 1.0 && set_gtk_property!(bt_vol_up, :sensitive, 1)
        vol = round(vol, digits=1)

        type = types[get_gtk_property(combo_type, :active, Int64) + 1]
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if type == "voice"
            if lang == "EN"
                if character == "negative"
                    if gender == "male"
                        snd = voices_en_m[rand(1:5, 1)[]]
                    else
                        snd = voices_en_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_en_m[rand(6:10, 1)[]]
                    else
                        snd = voices_en_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_en_m[rand(11:15, 1)[]]
                    else
                        snd = voices_en_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "DE"
                if character == "negative"
                    if gender == "male"
                        snd = voices_de_m[rand(1:5, 1)[]]
                    else
                        snd = voices_de_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_de_m[rand(6:10, 1)[]]
                    else
                        snd = voices_de_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_de_m[rand(11:15, 1)[]]
                    else
                        snd = voices_de_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "SP"
                if character == "negative"
                    if gender == "male"
                        snd = voices_sp_m[rand(1:5, 1)[]]
                    else
                        snd = voices_sp_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_sp_m[rand(6:10, 1)[]]
                    else
                        snd = voices_sp_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_sp_m[rand(11:15, 1)[]]
                    else
                        snd = voices_sp_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "PL"
                if character == "negative"
                    if gender == "male"
                        snd = voices_pl_m[rand(1:5, 1)[]]
                    else
                        snd = voices_pl_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_pl_m[rand(6:10, 1)[]]
                    else
                        snd = voices_pl_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_pl_m[rand(11:15, 1)[]]
                    else
                        snd = voices_pl_w[rand(11:15, 1)[]]
                    end
                end
            end
        elseif type == "whisper"
            snd = snd_whisper
        elseif type == "noise"
            snd = snd_noise
        elseif type == "ringing"
            snd = snd_sine
        end
        snd_tmp = deepcopy(snd)

        snd_tmp[1][:, 1] = snd[1][:, 1] .* (vol * (d_l * 0.25))
        snd_tmp[1][:, 2] = snd[1][:, 2] .* (vol * (d_r * 0.25))
        wavplay(snd_tmp[1], snd_tmp[2])
    end

    signal_connect(bt_play, "clicked") do widget
        type = types[get_gtk_property(combo_type, :active, Int64) + 1]
        lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
        character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
        gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
        if type == "voice"
            if lang == "EN"
                if character == "negative"
                    if gender == "male"
                        snd = voices_en_m[rand(1:5, 1)[]]
                    else
                        snd = voices_en_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_en_m[rand(6:10, 1)[]]
                    else
                        snd = voices_en_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_en_m[rand(11:15, 1)[]]
                    else
                        snd = voices_en_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "DE"
                if character == "negative"
                    if gender == "male"
                        snd = voices_de_m[rand(1:5, 1)[]]
                    else
                        snd = voices_de_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_de_m[rand(6:10, 1)[]]
                    else
                        snd = voices_de_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_de_m[rand(11:15, 1)[]]
                    else
                        snd = voices_de_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "SP"
                if character == "negative"
                    if gender == "male"
                        snd = voices_sp_m[rand(1:5, 1)[]]
                    else
                        snd = voices_sp_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_sp_m[rand(6:10, 1)[]]
                    else
                        snd = voices_sp_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_sp_m[rand(11:15, 1)[]]
                    else
                        snd = voices_sp_w[rand(11:15, 1)[]]
                    end
                end
            elseif lang == "PL"
                if character == "negative"
                    if gender == "male"
                        snd = voices_pl_m[rand(1:5, 1)[]]
                    else
                        snd = voices_pl_w[rand(1:5, 1)[]]
                    end
                elseif character == "positive"
                    if gender == "male"
                        snd = voices_pl_m[rand(6:10, 1)[]]
                    else
                        snd = voices_pl_w[rand(6:10, 1)[]]
                    end
                else
                    if gender == "male"
                        snd = voices_pl_m[rand(11:15, 1)[]]
                    else
                        snd = voices_pl_w[rand(11:15, 1)[]]
                    end
                end
            end
        elseif type == "whisper"
            snd = snd_whisper
        elseif type == "noise"
            snd = snd_noise
        elseif type == "ringing"
            snd = snd_sine
        end
        snd_tmp = deepcopy(snd)

        snd_tmp[1][:, 1] = snd[1][:, 1] .* (vol * (d_l * 0.25))
        snd_tmp[1][:, 2] = snd[1][:, 2] .* (vol * (d_r * 0.25))
        wavplay(snd_tmp[1], snd_tmp[2])
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = event.x
        y_pos = event.y
        Threads.@spawn begin
            @guarded draw(can) do widget
                ctx = getgc(widget)
                h = Cairo.height(can)
                w = Cairo.width(can)
                Cairo.rectangle(ctx, 0, 0, w, h)
                Cairo.set_source_rgb(ctx, 1, 1, 1)
                Cairo.fill(ctx)

                Cairo.set_source_surface(ctx, img, 1, 1)
                Cairo.paint(ctx)

                Gtk.arc(ctx, x_pos, y_pos, 10, 0, 2*pi)
                Gtk.set_source_rgb(ctx, 1, 0, 0)
                Gtk.stroke(ctx)
                Gtk.reveal(widget)
            end

            x_pos = round((x_pos / 800) - 0.5, digits=2)
            y_pos = round((y_pos / 800) - 0.5, digits=2)

            y_pos = -y_pos
            y_pos == -0.0 && (y_pos = 0)
            d_l = sqrt((x_pos - -0.05)^2 + (y_pos^2))
            d_r = sqrt((x_pos - 0.05)^2 + (y_pos^2))
            if d_l < 0.10 || d_r < 0.10
                d_l = 0
                d_r = 0
            end
            d_l = (1 - d_l)^4
            d_r = (1 - d_r)^4

            d_l > d_r && (d_r *= 0.75)
            d_r > d_l && (d_l *= 0.75)

            d_l = round(d_l, digits=2)
            d_r = round(d_r, digits=2)

            type = types[get_gtk_property(combo_type, :active, Int64) + 1]
            lang = langs[get_gtk_property(combo_lang, :active, Int64) + 1]
            character = characters[get_gtk_property(combo_character, :active, Int64) + 1]
            gender = genders[get_gtk_property(combo_gender, :active, Int64) + 1]
            if type == "voice"
                if lang == "EN"
                    if character == "negative"
                        if gender == "male"
                            snd = voices_en_m[rand(1:5, 1)[]]
                        else
                            snd = voices_en_w[rand(1:5, 1)[]]
                        end
                    elseif character == "positive"
                        if gender == "male"
                            snd = voices_en_m[rand(6:10, 1)[]]
                        else
                            snd = voices_en_w[rand(6:10, 1)[]]
                        end
                    else
                        if gender == "male"
                            snd = voices_en_m[rand(11:15, 1)[]]
                        else
                            snd = voices_en_w[rand(11:15, 1)[]]
                        end
                    end
                elseif lang == "DE"
                    if character == "negative"
                        if gender == "male"
                            snd = voices_de_m[rand(1:5, 1)[]]
                        else
                            snd = voices_de_w[rand(1:5, 1)[]]
                        end
                    elseif character == "positive"
                        if gender == "male"
                            snd = voices_de_m[rand(6:10, 1)[]]
                        else
                            snd = voices_de_w[rand(6:10, 1)[]]
                        end
                    else
                        if gender == "male"
                            snd = voices_de_m[rand(11:15, 1)[]]
                        else
                            snd = voices_de_w[rand(11:15, 1)[]]
                        end
                    end
                elseif lang == "SP"
                    if character == "negative"
                        if gender == "male"
                            snd = voices_sp_m[rand(1:5, 1)[]]
                        else
                            snd = voices_sp_w[rand(1:5, 1)[]]
                        end
                    elseif character == "positive"
                        if gender == "male"
                            snd = voices_sp_m[rand(6:10, 1)[]]
                        else
                            snd = voices_sp_w[rand(6:10, 1)[]]
                        end
                    else
                        if gender == "male"
                            snd = voices_sp_m[rand(11:15, 1)[]]
                        else
                            snd = voices_sp_w[rand(11:15, 1)[]]
                        end
                    end
                elseif lang == "PL"
                    if character == "negative"
                        if gender == "male"
                            snd = voices_pl_m[rand(1:5, 1)[]]
                        else
                            snd = voices_pl_w[rand(1:5, 1)[]]
                        end
                    elseif character == "positive"
                        if gender == "male"
                            snd = voices_pl_m[rand(6:10, 1)[]]
                        else
                            snd = voices_pl_w[rand(6:10, 1)[]]
                        end
                    else
                        if gender == "male"
                            snd = voices_pl_m[rand(11:15, 1)[]]
                        else
                            snd = voices_pl_w[rand(11:15, 1)[]]
                        end
                    end
                end
            elseif type == "whisper"
                snd = snd_whisper
            elseif type == "noise"
                snd = snd_noise
            elseif type == "ringing"
                snd = snd_sine
            end
            snd_tmp = deepcopy(snd)

            snd_tmp[1][:, 1] = snd[1][:, 1] .* (vol * (d_l * 0.25))
            snd_tmp[1][:, 2] = snd[1][:, 2] .* (vol * (d_r * 0.25))
            wavplay(snd_tmp[1], snd_tmp[2])
        end
    end

    signal_connect(bt_save, "clicked") do widget
        fname = save_dialog("Save as...", GtkNullContainer(), (GtkFileFilter("*.csv", name="All supported formats"), "*.txt"))
        if fname != ""
            f = open(fname, "w")
            println(f, "\"AH type\",$(types[get_gtk_property(combo_type, :active, Int64) + 1])")
            if get_gtk_property(combo_type, :active, Int64) == 0
                println(f, "\"AVH language\",$(langs[get_gtk_property(combo_lang, :active, Int64) + 1])")
                println(f, "\"AVH gender\",$(genders[get_gtk_property(combo_character, :active, Int64) + 1])")
                println(f, "\"AVH emotional aspect\",$(characters[get_gtk_property(combo_gender, :active, Int64) + 1])")
            else
                println(f, "\"AVH language\",NA")
                println(f, "\"AVH gender\",NA")
                println(f, "\"AVH emotional aspect\",NA")
            end
            println(f, "\"distance L\",$d_l")
            println(f, "\"distance R\",$d_r")
            print(f, "\"volume\",$vol")
            close(f)
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
        return nothing
    end

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end

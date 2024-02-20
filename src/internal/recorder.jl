function _beep()
    beep, fs = wavread(joinpath(res_path, "beep.wav"))
    wavplay(beep, fs)
end
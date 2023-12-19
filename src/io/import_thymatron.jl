export import_thymatron

"""
    import_thymatron(file_name)

Load scanned images of EEG printed by Thymatron ECT equipment. Usually there are three images: baseline, during seizure activity and after seizure activity. If more then one image is provided, images are added as consecutive channels in the same object.

Image properties:

- 100 DPI
- 490 x 100 px
- height 2.5 cm
- width 12.5 cm

# Arguments

- `file_name::String`: name of the file to load
- `dpi::Int64=100`: DPI of the scanned images

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

Wysokiński A. EEG_ADC: Digitizer and Analyzer of Electroconvulsive Therapy Paper Electroencephalogram Recordings. JECT 2022; 4: 255-256
"""
function import_thymatron(file_name::Union{String, Vector{String}})

    data_tmp = Vector{Vector{Float64}}()
    sampling_rate = 100 # = DPI

    typeof(file_name) == String && (file_name = [file_name])

    for file_idx in 1:length(file_name)

        @assert isfile(file_name[file_idx]) "File $(file_name[file_idx]) does not exist!"

        # load data
        img = Gray.(FileIO.load(file_name[file_idx]))

        size(img, 1) > size(img, 2) && (img = rotl90(img))

        size(img, 1) > 100 && (img = img[1:100, :])
        size(img, 2) > 490 && (img = img[:, 490])

        dimx = size(img, 1)
        dimy = size(img, 2)

        # adjust histogram
        # in some cases modify the `nbins` parameter for better results 
        # alg = Equalization(nbins = 3)
        # img = adjust_histogram(img, alg)

        # denoise image with median filter
        img = mapwindow(median, img, (3, 3))

        # unsharp mask
        # img = unsharp_mask(img, radius=5, amount=2)

        # edge detection
        # img = imedge(img, KernelFactors.sobel)[3]
        img = imedge(img, KernelFactors.sobel)[3]

        # binarize
        img_bin = binarize(img, Otsu())
        # img_bin = binarize(img, MinimumError())
        img_bin = dilate(img_bin)
        # img_bin = Bool.(img_bin)
        img_bin = Int.(img_bin)
        # fill in discontinuities
        for idx in 2:dimx
            c = img_bin[:, idx]
            sum(c) == 0 && (img_bin[:, idx] = img_bin[:, idx - 1])
        end
        # thinning
        for idx in 1:dimy
            c = img_bin[:, idx]
            l = sum(c) ÷ 2
            l_idx = findfirst(isequal(1), c)
            img_bin[:, idx] .= 0
            l_idx !== nothing && (img_bin[l_idx + l, idx] = 1)
        end

        # interpolate
        # TODO: input argument for user-defined px_uv and px_s
        t = 1:dimy
        signal = zeros(Int64, dimy)
        for idx in 1:dimy
            signal[idx] = dimx - findfirst(isequal(1), img_bin[:, idx])
        end
        s = round.(Int64, CubicSplineInterpolation(t, signal))
        img_bin = zeros(Int64, size(img_bin))
        for idx in 1:dimy
            img_bin[(dimx - s[idx]), idx] = 1
        end

        # digitize: convert to EEG signal
        px_cm = 0.03937008 * sampling_rate * 10      # 1 dpi = 0.03937008 px/mm => px/mm => px/cm
        px_uv = 200 / px_cm                          # 1 cm = 200 uV => convert to pixels
        px_s = 0.00025 * px_cm                       # 2.5 cm = 1 s => convert to pixels
        eeg_signal = zeros(dimy)
        eeg_time = zeros(dimy)
        for idx in 1:dimy
            eeg_signal[idx] = px_uv * ((dimx ÷ 2) - findfirst(isequal(1), img_bin[:, idx]))
            eeg_time[idx] = idx * px_s
        end
    
        eeg_signal = round.(eeg_signal, digits=3)

        push!(data_tmp, eeg_signal)

    end

    data = zeros(length(data_tmp), length(data_tmp[1]), 1)

    for idx in 1:length(file_name)
        data[idx, :, 1] = data_tmp[idx]
    end

    ch_n = size(data, 1)
    clabels = repeat(["ch"], ch_n)
    for idx in 1:ch_n
        clabels[idx] *= string(idx)
    end
    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type="eeg",
                              file_name=file_name[1],
                              file_size_mb=0,
                              file_type="Thymatron",
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes="",
                              channel_type=repeat(["eeg"], ch_n),
                              reference="physical",
                              clabels=clabels,
                              units=repeat(["μV"], ch_n),
                              transducers=repeat([""], ch_n),
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=ones(ch_n))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)
    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj
    
end

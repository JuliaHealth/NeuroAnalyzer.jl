#=
[eegplotter.jl]
Version =  0.022
Author = "William Herrera"
Copyright = "Copyright 2018 William Herrera"
Created = "12 Jan 2018"
Purpose = "EEG file routines viewer example"
=#


using EDFPlus
using DSP
using Plots
import FileIO
pyplot()
using PyPlot


linspace(start, stop, len) = LinRange{Float64}(start, stop, len)

function averagereference(edfh, channels=edfh.mapped_signals)
    data = EDFPlus.signaldata(edfh)
    recs = edfh.datarecords
    siglen = maximum(x->x.smp_per_record, edfh.signalparam)
    avgref = zeros(siglen*recs)
    spans = map(chan->EDFPlus.signalindices(edfh, chan), channels)
    chancount = length(channels)
    for recnum in 1:size(data)[1], channum in 1:chancount
        span = spans[channum]
        addlen = span[2] - span[1]
        addstart = (recnum-1)*siglen+1
        avgref[addstart:addstart+addlen] = data[recnum, span[1]:span[2]][:] ./ chancount
    end
    avgref
end


function eegpages(edfh; channels=collect(3:8), secsperpage=15.0)
    epages = Array{Array{Array{Float64,1},1},1}([])
    fs = samplerate(edfh, channels[1])
    starts = linspace(0.0, edfh.file_duration, div(edfh.file_duration,secsperpage))
    for t1 in starts[1:end-1]
        t2 = t1 + secsperpage
        epage = multichanneltimesegment(edfh, channels, t1, t2, true)
        for i in 1:size(epage)[1]
            epage[i] = highpassfilter(epage[i], fs, 1.5)
            epage[i] = lowpassfilter(epage[i], fs, 24.0)
        end
        push!(epages, epage)
    end
    avgref = averagereference(edfh, collect(1:8))
    spectpage = spectrogram(avgref)
    spectsegment = Plots.contourf(spectpage.time, spectpage.freq, log.(spectpage.power),
                        xaxis=false, yaxis=false, colorbar=false)
    imgname = "tmp.png"
    Plots.savefig(spectsegment, imgname)
    timepoints = linspace(0.0, secsperpage, length(epages[1][1]))
    channelnames = [trim(edfh.signalparam[chan].label) for chan in channels]
    epages, timepoints, channelnames, imgname, edfh.path
end


function plotpage(timepoints, page, ylabels, imgname, title)
    nchannels = size(page)[1]
    plt = Plots.plot(timepoints, page, layout=(nchannels+1,1),
                           xticks=collect(timepoints[1]:1:timepoints[end]),
                           yticks=false, legend=false, linetype=:line)
    for i in 1:nchannels
        Plots.plot!(yaxis=true, xaxis=false, subplot=i, tight_layout=true, ylabel = ylabels[i])
    end
    plt[nchannels][:xaxis][:showaxis] = true
    Plots.plot!(title = title, subplot=1, size=dims)
    img = FileIO.load(imgname)
    plot!(img, aspect_ratio="auto",subplot=nchannels+1)
    plt
end


function vieweeg(filename)
    currentpage = 1
    maxpage = 1
    channelcount = 8
    needupdate = false
    function onclick(event)
        println(event)
        println(event[:x], "  ", event[:y], "   ", event[:canvas][:get_width_height]())
        x = event[:x]
        y = event[:y]
        (width, height) = event[:canvas][:get_width_height]()
        global dims = (width, height)
        yfrac = div(height, channelcount+2)
        ybottom = div(yfrac, 2)
        ytop = yfrac + ybottom
        if y < ytop && y > ybottom && x < width && x > 0
            tmp = Int(round(x * maxpage / width))
            oldpage = currentpage
            currentpage = (tmp < 1) ? 1 : (tmp > maxpage ? maxpage : tmp)
            if oldpage != currentpage
                needupdate = true
            end
        end
    end
    function onkeypress(event)
        println(event)
        println(event[:key])
        (width, height) = event[:canvas][:get_width_height]()
        global dims = (width, height)
        tmp = currentpage
        if event[:key] in ("up", "right", "pageup", "+")
            tmp = currentpage + 1
        elseif event[:key] in ("down", "left", "pagedown", "-")
            tmp = currentpage - 1
        else
            return
        end
        oldpage = currentpage
        currentpage = (tmp < 1) ? 1 : (tmp > maxpage ? maxpage : tmp)
        if oldpage != currentpage
            needupdate = true
        end
    end
    edfh = loadfile(filename)
    epages, timepoints, channelnames, imgname, title = eegpages(edfh)
    maxpage = length(epages)
    plt = plotpage(timepoints, epages[2], channelnames, imgname, title)
    display(plt)
    fig = PyPlot.gcf()
    fig[:canvas][:mpl_connect]("button_press_event", onclick)
    fig[:canvas][:mpl_connect]("key_press_event", onkeypress)
    println("entering loop...")
    while true
        yield()
        if needupdate
            needupdate = false
            newpage = epages[currentpage]
            xdelta = (currentpage - 1) * timepoints[end]
            mins, secs = Int.(floor.(divrem(xdelta, 60)))
            newtitle = replace(title, r".+\/([^\/]+)$" => s"\1")
            newtitle *= " (Page $currentpage of $maxpage: at $mins minutes, $secs seconds)"
            display(plotpage(timepoints .+ xdelta, newpage, channelnames, imgname, newtitle))
            println("updated")
            yield()
        end
        sleep(0.1)
    end
end


global dims = (1200, 400)
const filename = "../test/EDFPlusTestFile.edf"
vieweeg(filename)

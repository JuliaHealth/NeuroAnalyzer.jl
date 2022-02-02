using EDFPlus
using Test


edfh = loadfile("EDFPlusTestFile.edf")

@test edfh.gender == "Female"
@test edfh.annotationchannel == 30
@test edfh.filetype == FileStatus(1)

fs = samplerate(edfh,1)
f1 = EDFPlus.recordslice(edfh, 21, 22)[1,:]
f2 = highpassfilter(reshape(f1, length(f1)), fs)
f3 = lowpassfilter(f2, fs)
f4 = notchfilter(f3, fs)
@test round(f4[end-3], digits=2) == 6803.28
@test EDFPlus.recordindexat(edfh, edfh.file_duration - 0.05) == edfh.datarecords
@test EDFPlus.epochmarkers(edfh, 10.525)[7] == (632, 10)

eegpages = epoch_iterator(edfh, 12, channels=[7,8,9])
annots = annotation_epoch_iterator(edfh, 12)

for (pgnum, page) in enumerate(eegpages)
    if pgnum == 1
        @test round(page[3][end], digits=3) == -24.121
    elseif pgnum == 30
        @test round(page[1][100], digits=3) == -87.012
    end
end

for (i, ann) in enumerate(annots)
    if i == 9
        @test ann[1][1].onset == 95.9
    end
end

closefile!(edfh)
@test edfh.filetype == EDFPlus.CLOSED


bdfh = loadfile("samplefrombiosemicom.bdf")

@test "$(bdfh.filetype)" == "BDF" && bdfh.filetype == EDFPlus.BDF
@test trim(bdfh.signalparam[4].label) == "A4"
@test bdfh.startdate_year == 2001
@test bdfh.signalparam[end].physdimension == "Boolean"
@test EDFPlus.recordslice(bdfh, 4, 14)[1,end-3] == 1835263

@test EDFPlus.version() == 0.59
sig = physicalchanneldata(bdfh, 1)
@test round(sig[100], digits=3) == 523818.031
@test EDFPlus.recordindexat(bdfh, bdfh.file_duration - 0.3) == bdfh.datarecords

true

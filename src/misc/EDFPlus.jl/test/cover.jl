using EDFPlus
using Test


edfh = loadfile("EDFPlusTestFile.edf")

EDFPlus.version() != nothing
epoch_iterator(edfh, 1) != nothing
annotation_epoch_iterator(edfh, 1) != nothing
channeltimesegment(edfh, 1, 0.5, 1.5, true) != nothing
multichanneltimesegment(edfh, [1], 0.25, 1.25, false) != nothing
sg = physicalchanneldata(edfh, 2) != nothing
(fs = samplerate(edfh, 2)) != nothing



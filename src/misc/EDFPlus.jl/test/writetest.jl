using EDFPlus
using Test

edfh = loadfile("EDFPlusTestFile.edf")
sz = size(edfh.EDFsignals)
@test sz == (20010,601)

eann = Annotation()
@test eann.onset == 0.0
ann = Annotation(61.04, "5.25", "They said schöner")
addannotation!(edfh, ann.onset, ann.duration, ann.annotation)
EDFPlus.latintoascii("")
edfh.gender = "Male"
edfh.signalparam[1].label = "Signal"
@test EDFPlus.readBiosemiStatus(edfh) isa Dict

newedfh = writefile!(edfh, "NEWedfplustestfile.edf")
@test size(newedfh.EDFsignals) == sz
@test EDFPlus.bytesperdatapoint(newedfh) == 2
closefile!(newedfh)

bdfh = loadfile("samplefrombiosemicom.bdf")
bsz = size(bdfh.BDFsignals)
@test bsz == (60, 34816)

newbdfh = writefile!(bdfh, "NEWsamplefrombiosemicom.bdf")
@test size(newbdfh.BDFsignals) == bsz
closefile!(newbdfh)
closefile!(bdfh)

newbdfh = writefile!(edfh, "NEWbdfplusfromedfplus.bdf", sigformat=EDFPlus.bdfplus)
@test size(newbdfh.BDFsignals) == sz
@test EDFPlus.bytesperdatapoint(newbdfh) == 3

newedfh2 = writefile!(newbdfh, "NEWedfplusfrombdfplus.edf", sigformat=EDFPlus.edfplus)
@test size(newedfh2.EDFsignals) == sz
@test EDFPlus.bytesperdatapoint(newedfh2) == 2

true

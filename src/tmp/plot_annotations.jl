plot((x, 0))
x = edf.eeg_annotations[!, :onset]
e = edf.eeg_annotations[!, :event]
using Plots
x = [1, 2, 3, 4]
e = ["AA", "BB", "CC", "DD"]

plot(xlims=(0, x[10]), xticks=x, ylims=(0, 1), yticks=[0, 1], tick_direction=:none, grid=false)
for idx in 1:10
    p = Plots.annotate!(x[idx], 0, Plots.text(" $(e[idx])", pointsize=3, halign=:left, valign=:center, rotation=90))
end
plot!()
export tes_model

# initial version, simplified model -- just superficial spread of the electric field
# next version: model spread of electric field at the cortical surface -- reduce charge for skull resistance

function tes_model(; anode::String, cathode::String, anode_curr::Real=2.0, cathode_curr::Real=-2.0)

    locs = import_locs(joinpath(NeuroAnalyzer.PATH, "locs", "standard-10-10-cap47.ced"))

    labels = locs[!, :labels]
    anode_ch = findfirst(isequal(anode), labels)
    cathode_ch = findfirst(isequal(cathode), labe)

    if cart == false
        loc_x = zeros(nrow(locs))
        loc_y = zeros(nrow(locs))
        for idx in 1:nrow(locs)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end
    loc_x = NeuroAnalyzer._s2v(loc_x)
    loc_y = NeuroAnalyzer._s2v(loc_y)
    anode_pos = (loc_x[anode_ch], loc_y[anode_ch])
    cathode_pos = (loc_x[cathode_ch], loc_y[cathode_ch])

    # kₑ = 8.988 * 10^9 # for air / vacuum
    kₑ = 1

    # 1st column: distance from anode
    # 2nd column: distance from cathode
    r = zeros(nrow(locs), 2)
    for idx in 1:nrow(locs)
        r[idx, 1] = euclidean((locs[!, :loc_x][idx], locs[!, :loc_y][idx]), anode_pos)
        r[idx, 2] = euclidean((locs[!, :loc_x][idx], locs[!, :loc_y][idx]), cathode_pos)
    end
    # replace distance for anode and cathode with 1 (should be 0, but it would result in division by 0 error)
    r[anode_ch, :] = [1, 1]
    r[cathode_ch, :] = [1, 1]
    r .+= 1

    E = zeros(nrow(locs))
    for idx in 1:nrow(locs)
        E[idx] = ((kₑ * anode_curr) / r[idx, 1]) + ((kₑ * cathode_curr) / r[idx, 2])
    end
    E[anode_ch] = anode_curr
    E[cathode_ch] = cathode_curr

    obj_tmp = ecog = create(data_type="eeg")
    obj_tmp.data = reshape(repeat(E, 1, 2), nrow(locs), 2, 1)
    obj_tmp.header.recording[:labels] = locs[!, :labels]
    obj_tmp.header.recording[:channel_type] = repeat(["EEG"], nrow(locs))
    obj_tmp.locs = locs
    create_time!(obj_tmp, fs=1)
    obj_tmp.time_pts
    plot_topo(obj_tmp, seg=(0, 1), title="", nmethod=:none)
    
    # Plots.plot(E, xticks=(1:nrow(locs), locs[!, :labels]), xtickfontsize=3)

end
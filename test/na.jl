using NeuroAnalyzer
using Test
const load_preference = NeuroAnalyzer.Preferences.load_preference

# capture what f() prints to stdout
function capture(f)
    mktemp() do path, io
        redirect_stdout(io) do
            f()
        end
        flush(io)
        read(path, String)
    end
end

@info "Test: messages"
v0, c0 = NeuroAnalyzer.verbose, NeuroAnalyzer.colors
NeuroAnalyzer.verbose = true
NeuroAnalyzer.colors = false
@test capture(() -> NeuroAnalyzer._info("abc")) == "[ Info: abc\n"
@test capture(() -> NeuroAnalyzer._warn("abc")) == "[ Warning: abc\n"
@test capture(() -> NeuroAnalyzer._error("abc")) == "[ Error: abc\n"
@test capture(() -> NeuroAnalyzer._deprecated("f")) == "[ Warning: Function f() is deprecated.\n"
@test capture(() -> NeuroAnalyzer._deprecated("f", "g")) == "[ Warning: Function f() is deprecated, please use g() instead.\n"
@test occursin("Work In Progress", capture(() -> NeuroAnalyzer._wip()))
@test capture(() -> NeuroAnalyzer._msg("X", NeuroAnalyzer.Crayons.Box.RED_FG, "y")) == "[ X: y\n"
NeuroAnalyzer.colors = true
@test occursin("abc", capture(() -> NeuroAnalyzer._info("abc")))
# verbose = false: nothing is printed
NeuroAnalyzer.verbose = false
@test capture(() -> NeuroAnalyzer._info("abc")) == ""
NeuroAnalyzer.verbose, NeuroAnalyzer.colors = v0, c0
# logging on/off
NeuroAnalyzer._log_off()
@test NeuroAnalyzer.Logging.min_enabled_level(NeuroAnalyzer.Logging.current_logger()) > NeuroAnalyzer.Logging.Warn ||
      (@test_logs @warn("hidden"); true)
NeuroAnalyzer._log_on()

@info "Test: version() / info()"
@test version() == "$(NeuroAnalyzer.VER.major).$(NeuroAnalyzer.VER.minor).$(NeuroAnalyzer.VER.patch)"
out = capture(info)
@test occursin("NeuroAnalyzer: $(NeuroAnalyzer.VER)", out)
@test occursin("Imported packages:", out)
@test occursin("DSP", out)
@test occursin("Threads:", capture(NeuroAnalyzer._print_settings))
@test capture(NeuroAnalyzer._print_thread_warnings) isa String
@test occursin("FFTW", capture(NeuroAnalyzer._print_package_versions))

@info "Test: set_*()"
# settings are stored as preferences: set and restore
pb0, vb0, eb0, cl0 = NeuroAnalyzer.progress_bar, NeuroAnalyzer.verbose, NeuroAnalyzer.exclude_bads, NeuroAnalyzer.colors
set_progress_bar(!pb0)
@test NeuroAnalyzer.progress_bar == !pb0
@test load_preference(NeuroAnalyzer, "progress_bar") == !pb0
set_verbose(!vb0)
@test NeuroAnalyzer.verbose == !vb0
set_exclude_bads(!eb0)
@test NeuroAnalyzer.exclude_bads == !eb0
set_colors(!cl0)
@test NeuroAnalyzer.colors == !cl0
@test load_preference(NeuroAnalyzer, "colors") == !cl0
set_prefs(; progress_bar = pb0, verbose = vb0, exclude_bads = eb0, colors = cl0)
@test (NeuroAnalyzer.progress_bar, NeuroAnalyzer.verbose, NeuroAnalyzer.exclude_bads, NeuroAnalyzer.colors) == (pb0, vb0, eb0, cl0)
@test load_preference(NeuroAnalyzer, "verbose") == vb0

@info "Test: plugins"
# a temporary plugins folder (the user's plugins are not touched)
pp0 = NeuroAnalyzer.plugins_path
tmp = mktempdir()
NeuroAnalyzer.plugins_path = tmp
mkpath(joinpath(tmp, "naplugin", "src"))
write(joinpath(tmp, "naplugin", "src", "naplugin.jl"), "na_test_plugin_fn() = 42\n")
NeuroAnalyzer.plugins_reload()
@test NeuroAnalyzer.na_test_plugin_fn() == 42
@test occursin("1. naplugin", capture(plugins_list))
# install from a .tar.gz archive
arch = joinpath(mktempdir(), "naplugin2.tar.gz")
src2 = mktempdir()
mkpath(joinpath(src2, "naplugin2", "src"))
write(joinpath(src2, "naplugin2", "src", "naplugin2.jl"), "na_test_plugin_fn2() = 43\n")
run(Cmd(`tar -czf $arch naplugin2`; dir = src2))
plugins_install(arch)
@test isdir(joinpath(tmp, "naplugin2"))
@test NeuroAnalyzer.na_test_plugin_fn2() == 43
@test_throws ArgumentError NeuroAnalyzer._install_from_archive(joinpath(tmp, "naplugin2.rar"))
write(joinpath(tmp, "x.rar"), "x")
@test_throws ArgumentError NeuroAnalyzer._install_from_archive(joinpath(tmp, "x.rar"))
# install from a git repository (a local one: no network)
repo = joinpath(mktempdir(), "naplugin3")
mkpath(joinpath(repo, "src"))
write(joinpath(repo, "src", "naplugin3.jl"), "na_test_plugin_fn3() = 44\n")
run(`$(NeuroAnalyzer.git()) -C $repo init -q`)
run(`$(NeuroAnalyzer.git()) -C $repo -c user.name=t -c user.email=t@t add .`)
run(`$(NeuroAnalyzer.git()) -C $repo -c user.name=t -c user.email=t@t commit -qm init`)
NeuroAnalyzer._install_from_remote(repo)
@test isfile(joinpath(tmp, "naplugin3", "src", "naplugin3.jl"))
@test_throws ErrorException NeuroAnalyzer._install_from_remote(joinpath(tmp, "no_such_repo"))
# a git plugin is updated (git pull)
@test NeuroAnalyzer._update_plugin(joinpath(tmp, "naplugin3")) === nothing
# update: plugins that are not git repositories are skipped
@test plugins_update() === nothing
@test NeuroAnalyzer._update_plugin(joinpath(tmp, "naplugin")) === nothing
@test_throws ArgumentError plugins_update("no_such_plugin")
# remove
plugins_remove("naplugin2")
@test !isdir(joinpath(tmp, "naplugin2"))
@test_throws ArgumentError plugins_remove("no_such_plugin")
NeuroAnalyzer.plugins_path = joinpath(tmp, "no_such_folder")
@test_throws ArgumentError plugins_list()
@test_throws ArgumentError plugins_reload()
NeuroAnalyzer.plugins_path = pp0
rm(tmp; recursive = true)

@info "Test: recorder / tester object builders"
e = NeuroAnalyzer._build_eda_object([1.0, 2.0, 3.0, 4.0])
@test datatype(e) == "eda"
@test e.data[1, :, 1] == [1.0, 2.0, 3.0, 4.0]
@test labels(e) == ["eda1"] && e.header.recording[:unit] == ["μS"]
@test sr(e) == NeuroAnalyzer._EDA_FS
@test_throws ArgumentError NeuroAnalyzer._build_eda_object(Float64[])
o = NeuroAnalyzer._build_tpt_object([1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 10.0], [11.0, 12.0])
@test datatype(o) == "tpt"
@test o.data[:, :, 1] == [1.0 2.0; 3.0 4.0; 5.0 6.0; 7.0 8.0; 9.0 10.0; 11.0 12.0]
@test labels(o) == ["pos_x", "pos_y", "pos_z", "acc_x", "acc_y", "acc_z"]
@test sr(o) == NeuroAnalyzer._TPT_FS
@test NeuroAnalyzer._ftt_check(5, 2, 1, -1, "") === nothing
@test_throws ArgumentError NeuroAnalyzer._ftt_check(0, 2, 1, -1, "")
@test_throws ArgumentError NeuroAnalyzer._ftt_check(5, 0, 1, -1, "")
@test_throws ArgumentError NeuroAnalyzer._ftt_check(5, 2, -1, -1, "")
@test_throws ArgumentError NeuroAnalyzer._ftt_check(5, 2, 1, -1, "/dev/ttyUSB0")

true

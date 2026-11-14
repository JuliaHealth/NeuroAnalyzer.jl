# Rewrite Artifacts.toml of NeuroAnalyzer.jl from the local artifact archives.
#
# For every artifact the archive neuroanalyzer-<name>[-dev].tar.gz is read from the artifacts
# directory and its entry is replaced: git-tree-sha1, download URL and sha256.
#
# usage: julia update_artifacts.jl <dev|stable> [--keep-tree]
#   dev          use the "-dev" archives
#   stable       use the stable archives
#   --keep-tree  fail if the contents of an artifact (git-tree-sha1) would change; used for
#                stable releases, which must ship the artifacts that devel was tested with
#
# environment:
#   NA_DIR            NeuroAnalyzer.jl checkout (default: parent of this script directory)
#   NA_ARTIFACTS_DIR  directory with the archives (default: ../NeuroAnalyzer-artifacts next to
#                     the NeuroAnalyzer.jl checkout)

using Pkg.Artifacts: bind_artifact!
using Pkg.TOML
using SHA
using Tar

const BASE_URL = "https://neuroanalyzer.org/artifacts"
# Artifacts.toml entry => archive base name (= top-level directory inside the archive)
const ARTIFACTS = ["NeuroAnalyzer_resources" => "neuroanalyzer-resources", "NeuroAnalyzer_test-files" => "neuroanalyzer-test-files"]

function main(args)
    channel = isempty(args) ? "" : args[1]
    channel in ("dev", "stable") || error("usage: julia update_artifacts.jl <dev|stable> [--keep-tree]")
    keep_tree = "--keep-tree" in args[2:end]

    # NA_DIR: set when the script is run from a copy (release scripts)
    na_dir = get(ENV, "NA_DIR", dirname(@__DIR__))
    toml_file = joinpath(na_dir, "Artifacts.toml")
    art_dir = get(ENV, "NA_ARTIFACTS_DIR", joinpath(dirname(na_dir), "NeuroAnalyzer-artifacts"))
    isfile(toml_file) || error("$toml_file not found")
    isdir(art_dir) || error("artifacts directory $art_dir not found (set NA_ARTIFACTS_DIR)")

    old = TOML.parsefile(toml_file)
    suffix = channel == "dev" ? "-dev" : ""

    # check all archives before Artifacts.toml is modified
    for (_, base) in ARTIFACTS
        archive = joinpath(art_dir, "$base$suffix.tar.gz")
        isfile(archive) || error("archive not found: $archive")
    end

    for (name, base) in ARTIFACTS
        archive = joinpath(art_dir, "$base$suffix.tar.gz")
        url = "$BASE_URL/$base$suffix.tar.gz"

        @info "$name: $archive"
        sha = bytes2hex(open(sha256, archive))

        # the tree hash takes long for large archives: reuse it if the archive is unchanged
        old_entry = get(old, name, Dict{String, Any}())
        old_tree = get(old_entry, "git-tree-sha1", "")
        old_sha = [get(d, "sha256", "") for d in get(old_entry, "download", [])]
        tree = if !isempty(old_tree) && sha in old_sha
            old_tree
        else
            # git-tree-sha1 is the git hash of the unpacked tree (empty directories are ignored)
            open(`gzip -dc $archive`) do io
                Tar.tree_hash(io; skip_empty = true)
            end
        end

        if !isempty(old_tree) && tree != old_tree
            keep_tree && error(
                "$name: contents differ from those referenced in Artifacts.toml ($old_tree → $tree). Publish the artifacts on devel first (scripts/publish_artifacts.sh).",
            )
            @info "$name: contents changed ($old_tree → $tree)"
        end

        bind_artifact!(toml_file, name, Base.SHA1(tree); download_info = [(url, sha)], force = true)
        println("  git-tree-sha1 = $tree")
        println("  url           = $url")
        println("  sha256        = $sha")
    end

    return nothing
end

main(ARGS)

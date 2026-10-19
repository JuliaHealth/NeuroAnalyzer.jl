using Pkg

@info "Activating LocalRegistry package"
Pkg.add("LocalRegistry")
using LocalRegistry

reg_path = expanduser(ARGS[1])
@assert isdir(reg_path) "The LocalRegistry folder does not exist or was not provided."
@info "Updating LocalRegistry"
Pkg.activate(".")
register(registry = reg_path, push = false)
@info "Done"

exit()

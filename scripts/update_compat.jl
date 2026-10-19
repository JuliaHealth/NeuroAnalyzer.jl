using Pkg

@info "Activating PackageCompatUI package"
Pkg.add("PackageCompatUI")
using PackageCompatUI

Pkg.activate(".")
compat_ui()

exit()

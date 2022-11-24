


![](assets/neuroanalyzer.png)


<a id='NeuroAnalyzer.jl-documentation'></a>

<a id='NeuroAnalyzer.jl-documentation-1'></a>

# NeuroAnalyzer.jl documentation


This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).


<a id='NeuroAnalyzer'></a>

<a id='NeuroAnalyzer-1'></a>

## NeuroAnalyzer

<a id='NeuroAnalyzer.na_info-Tuple{}' href='#NeuroAnalyzer.na_info-Tuple{}'>#</a>
**`NeuroAnalyzer.na_info`** &mdash; *Method*.



```julia
na_info()
```

Show NeuroAnalyzer and imported packages versions.

<a id='NeuroAnalyzer.na_plugins_reload-Tuple{}' href='#NeuroAnalyzer.na_plugins_reload-Tuple{}'>#</a>
**`NeuroAnalyzer.na_plugins_reload`** &mdash; *Method*.



```julia
na_plugins_reload()
```

Reload NeuroAnalyzer plugins.

<a id='NeuroAnalyzer.na_plugins_list-Tuple{}' href='#NeuroAnalyzer.na_plugins_list-Tuple{}'>#</a>
**`NeuroAnalyzer.na_plugins_list`** &mdash; *Method*.



```julia
na_plugins_list()
```

List NeuroAnalyzer plugins.

<a id='NeuroAnalyzer.na_plugins_remove-Tuple{String}' href='#NeuroAnalyzer.na_plugins_remove-Tuple{String}'>#</a>
**`NeuroAnalyzer.na_plugins_remove`** &mdash; *Method*.



```julia
na_plugins_remove(plugin)
```

Remove NeuroAnalyzer `plugin`.

**Arguments**

  * `plugin::String`: plugin name

<a id='NeuroAnalyzer.na_plugins_install-Tuple{String}' href='#NeuroAnalyzer.na_plugins_install-Tuple{String}'>#</a>
**`NeuroAnalyzer.na_plugins_install`** &mdash; *Method*.



```julia
na_plugins_install(plugin)
```

Install NeuroAnalyzer `plugin`.

**Arguments**

  * `plugin::String`: plugin Git repository URL

<a id='NeuroAnalyzer.na_plugins_update' href='#NeuroAnalyzer.na_plugins_update'>#</a>
**`NeuroAnalyzer.na_plugins_update`** &mdash; *Function*.



```julia
na_plugins_update(plugin)
```

Install NeuroAnalyzer `plugin`.

**Arguments**

  * `plugin::String`: plugin to update; if empty, update all

<a id='NeuroAnalyzer.na_set_use_cuda-Tuple{Bool}' href='#NeuroAnalyzer.na_set_use_cuda-Tuple{Bool}'>#</a>
**`NeuroAnalyzer.na_set_use_cuda`** &mdash; *Method*.



```julia
na_set_use_cuda(use_cuda)
```

Change `use_cuda` preference.

**Arguments**

  * `use_cuda::Bool`: value

<a id='NeuroAnalyzer.na_set_progress_bar-Tuple{Bool}' href='#NeuroAnalyzer.na_set_progress_bar-Tuple{Bool}'>#</a>
**`NeuroAnalyzer.na_set_progress_bar`** &mdash; *Method*.



```julia
na_set_progress_bar(progress_bar)
```

Change `progress_bar` preference.

**Arguments**

  * `progress_bar::Bool`: value

<a id='NeuroAnalyzer.na_set_plugins_path-Tuple{String}' href='#NeuroAnalyzer.na_set_plugins_path-Tuple{String}'>#</a>
**`NeuroAnalyzer.na_set_plugins_path`** &mdash; *Method*.



```julia
na_set_plugins_path(plugins_path)
```

Change `plugins_path` preference.

**Arguments**

  * `plugins_path::String`: value

<a id='NeuroAnalyzer.na_set_prefs-Tuple{}' href='#NeuroAnalyzer.na_set_prefs-Tuple{}'>#</a>
**`NeuroAnalyzer.na_set_prefs`** &mdash; *Method*.



```julia
na_set_prefs(use_cuda, plugins_path, progress_bar, verbose)
```

Save NeuroAnalyzer preferences.

**Arguments**

  * `use_cuda::Bool`
  * `plugins_path::String`
  * `progress_bar::Bool`
  * `verbose::Bool`

<a id='NeuroAnalyzer.na_set_verbose-Tuple{Bool}' href='#NeuroAnalyzer.na_set_verbose-Tuple{Bool}'>#</a>
**`NeuroAnalyzer.na_set_verbose`** &mdash; *Method*.



```julia
na_set_verbose(verbose)
```

Change `verbose` preference.

**Arguments**

  * `verbose::Bool`: value

<a id='NeuroAnalyzer.na_version-Tuple{}' href='#NeuroAnalyzer.na_version-Tuple{}'>#</a>
**`NeuroAnalyzer.na_version`** &mdash; *Method*.



```julia
na_version()
```

Convert NeuroAnalyzer version to string.

**Returns**

  * `na_ver::String`


<a id='Low-level-functions'></a>

<a id='Low-level-functions-1'></a>

## Low-level functions

<a id='NeuroAnalyzer.linspace-Tuple{Number, Number, Int64}' href='#NeuroAnalyzer.linspace-Tuple{Number, Number, Int64}'>#</a>
**`NeuroAnalyzer.linspace`** &mdash; *Method*.



```julia
linspace(start, stop, length)
```

Generates a sequence of evenly spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Number`
  * `stop::Number`
  * `n::Int64`: sequence length

**Returns**

  * `range::Vector`

<a id='NeuroAnalyzer.logspace-Tuple{Number, Number, Int64}' href='#NeuroAnalyzer.logspace-Tuple{Number, Number, Int64}'>#</a>
**`NeuroAnalyzer.logspace`** &mdash; *Method*.



```julia
logspace(start, stop, n)
```

Generates a sequence of log10-spaced numbers between `start` and `stop`.

**Arguments**

  * `start::Number`
  * `stop::Number`
  * `n::Int64`: sequence length

**Returns**

  * `range::Vector{<:Number}`

<a id='NeuroAnalyzer.m_pad0-Tuple{Matrix{<:Number}}' href='#NeuroAnalyzer.m_pad0-Tuple{Matrix{<:Number}}'>#</a>
**`NeuroAnalyzer.m_pad0`** &mdash; *Method*.



```julia
m_pad0(m)
```

Pad matrix with zeros to make it square.

**Arguments**

  * `m::Matrix{<:Number}`

**Returns**

  * `m::Matrix{<:Number}`

<a id='NeuroAnalyzer.vsearch-Tuple{Real, AbstractVector}' href='#NeuroAnalyzer.vsearch-Tuple{Real, AbstractVector}'>#</a>
**`NeuroAnalyzer.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; return_distance)
```

Return the positions of the `y` value in the vector `x`.

**Arguments**

  * `y::Real`: value of interest
  * `x::AbstractVector`: vector to search within
  * `acc::Bool=false`: if true, return the difference between `y` and `x[y_idx]`

**Returns**

  * `y_idx::Int64`
  * `y_dist::Real`: the difference between `y` and `x[y_idx]`

<a id='NeuroAnalyzer.vsearch-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.vsearch-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.vsearch`** &mdash; *Method*.



```julia
vsearch(y, x; acc)
```

Return the positions of the `y` vector in the vector `x`.

**Arguments**

  * `y::AbstractVector`: vector of interest
  * `x::AbstractVector`: vector to search within
  * `acc::Bool=false`: if true, return the difference between `y` and `x[y_idx:y_idx + length(y)]`

**Returns**

  * `y_idx::Int64`
  * `y_dist::Real`: the difference between `y` and `x[y_idx:y_idx + length(y)]`

<a id='NeuroAnalyzer.cart2pol-Tuple{Real, Real}' href='#NeuroAnalyzer.cart2pol-Tuple{Real, Real}'>#</a>
**`NeuroAnalyzer.cart2pol`** &mdash; *Method*.



```julia
cart2pol(x, y)
```

Convert cartographic coordinates to polar.

**Arguments**

  * `x::Real`
  * `y::Real`

**Returns**

  * `radius::Float64`
  * `theta::Float64`

<a id='NeuroAnalyzer.pol2cart-Tuple{Real, Real}' href='#NeuroAnalyzer.pol2cart-Tuple{Real, Real}'>#</a>
**`NeuroAnalyzer.pol2cart`** &mdash; *Method*.



```julia
pol2cart(radius, theta)
```

Convert polar coordinates to cartographic.

**Arguments**

  * `radius::Real`: polar radius, the distance from the origin to the point, in degrees
  * `theta::Real`: polar angle

**Returns**

  * `x::Float64`
  * `y::Float64`

<a id='NeuroAnalyzer.sph2cart-Tuple{Real, Real, Real}' href='#NeuroAnalyzer.sph2cart-Tuple{Real, Real, Real}'>#</a>
**`NeuroAnalyzer.sph2cart`** &mdash; *Method*.



```julia
sph2cart(radius, theta, phi)
```

Convert spherical coordinates to cartographic.

**Arguments**

  * `radius::Real`: spherical radius, the distance from the origin to the point
  * `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `phi::Real`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

**Returns**

  * `x::Float64`
  * `y::Float64`
  * `z::Float64`

<a id='NeuroAnalyzer.cart2sph-Tuple{Real, Real, Real}' href='#NeuroAnalyzer.cart2sph-Tuple{Real, Real, Real}'>#</a>
**`NeuroAnalyzer.cart2sph`** &mdash; *Method*.



```julia
cart2sph(x, y, z)
```

Convert spherical coordinates to cartographic.

**Arguments**

  * `x::Real`
  * `y::Real`
  * `z::Real`

**Returns**

  * `radius::Float64`: spherical radius, the distance from the origin to the point
  * `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `phi::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

<a id='NeuroAnalyzer.generate_window-Tuple{Symbol, Int64}' href='#NeuroAnalyzer.generate_window-Tuple{Symbol, Int64}'>#</a>
**`NeuroAnalyzer.generate_window`** &mdash; *Method*.



```julia
generate_window(type, n; even)
```

Return the `n`-point long symmetric window `type`.

**Arguments**

  * `type::Symbol`: window type:

      * `:hann`: Hann
      * `:bh`: Blackman-Harris
      * `:bohman`: Bohman
      * `:flat`: Flat-top window
      * `:bn`: Blackman-Nuttall
      * `:nutall`: Nuttall
      * `:triangle`: symmetric triangle (left half ↑, right half ↓)
      * `:exp`: symmetric exponential (left half ↑, right half ↓)
  * `n::Int64`: window length
  * `even::Bool=false`: if true, make the window of even length (+1 for odd n)

**Returns**

  * `w::Vector{Float64}`:: generated window

<a id='NeuroAnalyzer.fft0' href='#NeuroAnalyzer.fft0'>#</a>
**`NeuroAnalyzer.fft0`** &mdash; *Function*.



```julia
fft0(x, n)
```

Zeros-padded FFT.

**Arguments**

  * `x::AbstractArray`
  * `n::Int64`: number of zeros to add

**Returns**

  * `fft0::Vector{ComplexF64}`

<a id='NeuroAnalyzer.ifft0' href='#NeuroAnalyzer.ifft0'>#</a>
**`NeuroAnalyzer.ifft0`** &mdash; *Function*.



```julia
ifft0(x, n)
```

Zeros-padded IFFT.

**Arguments**

  * `x::AbstractArray`
  * `n::Int64`: number of zeros to add

**Returns**

  * `ifft0::Vector{ComplexF64}`

<a id='NeuroAnalyzer.fft2-Tuple{AbstractArray}' href='#NeuroAnalyzer.fft2-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.fft2`** &mdash; *Method*.



```julia
fft2(x)
```

Zeros-padded FFT, so the length of padded vector is a power of 2.

**Arguments**

  * `x::AbstractArray`

**Returns**

  * `fft2::Vector{ComplexF64}`

<a id='NeuroAnalyzer.ifft2-Tuple{AbstractArray}' href='#NeuroAnalyzer.ifft2-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.ifft2`** &mdash; *Method*.



```julia
ifft2(x)
```

Calculate IFFT for the vector `x` padded with zeros, so the length of padded `x` is a power of 2.

**Arguments**

  * `x::AbstractArray`

**Returns**

  * `ifft0::Vector{ComplexF64}`

<a id='NeuroAnalyzer.nextpow2-Tuple{Int64}' href='#NeuroAnalyzer.nextpow2-Tuple{Int64}'>#</a>
**`NeuroAnalyzer.nextpow2`** &mdash; *Method*.



```julia
nextpow2(x)
```

Return the next power of 2 for given number `x`.

**Argument**

  * `x::Int64`

**Returns**

  * `nextpow::Int64`

<a id='NeuroAnalyzer.vsplit' href='#NeuroAnalyzer.vsplit'>#</a>
**`NeuroAnalyzer.vsplit`** &mdash; *Function*.



```julia
vsplit(x, n)
```

Splits the vector `x` into `n`-long pieces.

**Argument**

  * `x::AbstractVector`
  * `n::Int64`

**Returns**

  * `x::Vector{AbstractVector}`

<a id='NeuroAnalyzer.s_rms-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_rms-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_rms`** &mdash; *Method*.



```julia
s_rms(signal)
```

Calculate Root Mean Square.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * rms::Float64`

<a id='NeuroAnalyzer.generate_sine' href='#NeuroAnalyzer.generate_sine'>#</a>
**`NeuroAnalyzer.generate_sine`** &mdash; *Function*.



```julia
generate_sine(f, t, a, p)
```

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

**Arguments**

  * `f::Real`: frequency
  * `t::Union{AbstractVector, AbstractRange}`: time vector
  * `a::Real`: amplitude
  * `p::Real`: initial phase

**Returns**

  * sine::Vector{Float64}`

<a id='NeuroAnalyzer.s_freqs-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_freqs-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_freqs`** &mdash; *Method*.



```julia
s_freqs(t)
```

Return vector of frequencies and Nyquist frequency for time vector.

**Arguments**

  * `t::AbstractVector, AbstractRange}`

**Returns**

  * `hz::Vector{Float64}`
  * `nyquist_freq::Float64`

<a id='NeuroAnalyzer.s_freqs-Tuple{Vector{Float64}, Int64}' href='#NeuroAnalyzer.s_freqs-Tuple{Vector{Float64}, Int64}'>#</a>
**`NeuroAnalyzer.s_freqs`** &mdash; *Method*.



```julia
s_freqs(signal, fs)
```

Return vector of frequencies and Nyquist frequency for signal.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`

**Returns**

  * `hz::Vector{Float64`
  * `nyquist_freq::Float64`

<a id='NeuroAnalyzer.m_sortperm-Tuple{Matrix}' href='#NeuroAnalyzer.m_sortperm-Tuple{Matrix}'>#</a>
**`NeuroAnalyzer.m_sortperm`** &mdash; *Method*.



```julia
m_sortperm(m; rev=false, dims=1)
```

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_idx::Matrix{Int64}`

<a id='NeuroAnalyzer.m_sort-Tuple{Matrix, Vector{Int64}}' href='#NeuroAnalyzer.m_sort-Tuple{Matrix, Vector{Int64}}'>#</a>
**`NeuroAnalyzer.m_sort`** &mdash; *Method*.



```julia
m_sort(m, m_idx; rev=false, dims=1)
```

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

**Arguments**

  * `m::Matrix`
  * `m_idx::Vector{Int64}`
  * `rev::Bool`
  * `dims::Int64`

**Returns**

  * `m_sorted::Matrix`

<a id='NeuroAnalyzer.pad0-Tuple{AbstractArray, Int64}' href='#NeuroAnalyzer.pad0-Tuple{AbstractArray, Int64}'>#</a>
**`NeuroAnalyzer.pad0`** &mdash; *Method*.



```julia
pad0(x, n)
```

Pad vector / rows of matrix / array with zeros. Works with 1-, 2- and 3-dimensional arrays.

**Arguments**

  * `x::Union{AbstractVector, AbstractArray`
  * `n::Int64`: number of zeros to add

**Returns**

  * `v_pad::Union{AbstractVector, AbstractArray`

<a id='NeuroAnalyzer.pad2-Tuple{AbstractArray}' href='#NeuroAnalyzer.pad2-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.pad2`** &mdash; *Method*.



```julia
pad2(x)
```

Pad vector / rows of matrix / array with zeros to the nearest power of 2 length.

**Arguments**

  * `x::Union{AbstractVector, AbstractArray`

**Returns**

  * `v_pad::Union{AbstractVector, AbstractArray`

<a id='NeuroAnalyzer.hz2rads-Tuple{Real}' href='#NeuroAnalyzer.hz2rads-Tuple{Real}'>#</a>
**`NeuroAnalyzer.hz2rads`** &mdash; *Method*.



```julia
hz2rads(f)
```

Convert frequency `f` in Hz to rad/s.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroAnalyzer.rads2hz-Tuple{Real}' href='#NeuroAnalyzer.rads2hz-Tuple{Real}'>#</a>
**`NeuroAnalyzer.rads2hz`** &mdash; *Method*.



```julia
rads2hz(f)
```

Convert frequency `f` in rad/s to Hz.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroAnalyzer.cmax-Tuple{Vector{ComplexF64}}' href='#NeuroAnalyzer.cmax-Tuple{Vector{ComplexF64}}'>#</a>
**`NeuroAnalyzer.cmax`** &mdash; *Method*.



```julia
cmax(x)
```

Return maximum value of the complex vector`x`.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmax::ComplexF64`

<a id='NeuroAnalyzer.cmin-Tuple{Vector{ComplexF64}}' href='#NeuroAnalyzer.cmin-Tuple{Vector{ComplexF64}}'>#</a>
**`NeuroAnalyzer.cmin`** &mdash; *Method*.



```julia
cmin(x)
```

Return minimum value of the complex vector`x`.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmin::ComplexF64`

<a id='NeuroAnalyzer.generate_sinc' href='#NeuroAnalyzer.generate_sinc'>#</a>
**`NeuroAnalyzer.generate_sinc`** &mdash; *Function*.



```julia
generate_sinc(t, f, peak, norm)
```

Generate sinc function.

**Arguments**

  * `t::AbstractRange=-2:0.01:2`: time
  * `f::Real=10.0`: frequency
  * `peak::Real=0`: sinc peak time
  * `norm::Bool=true`: generate normalized function

**Returns**

  * `sinc::Vector{Float64}`

<a id='NeuroAnalyzer.generate_morlet' href='#NeuroAnalyzer.generate_morlet'>#</a>
**`NeuroAnalyzer.generate_morlet`** &mdash; *Function*.



```julia
generate_morlet(fs, f, t; ncyc, complex)
```

Generate Morlet wavelet.

**Arguments**

  * `fs::Int64`: sampling rate
  * `f::Real`: frequency
  * `t::Real=1`: length = -t:1/fs:t
  * `ncyc::Int64=5`: number of cycles
  * `complex::Bool=false`: generate complex Morlet

**Returns**

  * `morlet::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroAnalyzer.generate_gaussian' href='#NeuroAnalyzer.generate_gaussian'>#</a>
**`NeuroAnalyzer.generate_gaussian`** &mdash; *Function*.



```julia
generate_gaussian(fs, f, t; ncyc, a)
```

Generate Gaussian wave.

**Arguments**

  * `fs::Int64`: sampling rate
  * `f::Real`: frequency
  * `t::Real=1`: length = -t:1/fs:t
  * `ncyc::Int64`: : number of cycles, width, SD of the Gaussian
  * `a::Real=1`: peak amp

**Returns**

  * `gaussian::Vector{Float64}`

<a id='NeuroAnalyzer.tuple_order' href='#NeuroAnalyzer.tuple_order'>#</a>
**`NeuroAnalyzer.tuple_order`** &mdash; *Function*.



```julia
tuple_order(t, rev)
```

Order tuple elements in ascending or descending (rev=true) order.

**Arguments**

  * `t::Tuple{Real, Real}`
  * `rev::Bool=false`

**Returns**

  * `t::Tuple{Real, Real}`

<a id='NeuroAnalyzer.s2_rmse-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_rmse-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_rmse`** &mdash; *Method*.



```julia
s2_rmse(signal1, signal2)
```

Calculate RMSE between two signals.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`

**Returns**

  * `r::Float64`

<a id='NeuroAnalyzer.m_norm-Tuple{AbstractArray}' href='#NeuroAnalyzer.m_norm-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.m_norm`** &mdash; *Method*.



```julia
m_norm(m)
```

Normalize matrix.

**Arguments**

  * `m::AbstractArray`

**Returns**

  * `m_norm::AbstractArray`

<a id='NeuroAnalyzer.s_cov-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_cov-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_cov`** &mdash; *Method*.



s_cov(signal; norm=true)

Calculate covariance of `signal * signal'`.

**Arguments**

  * `signal::AbstractVector`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='NeuroAnalyzer.s2_cov-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_cov-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_cov`** &mdash; *Method*.



s2_cov(signal1, signal2; norm=true)

Calculate covariance between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`
  * `norm::Bool`: normalize covariance

**Returns**

  * `cov_mat::Matrix{Float64}`

<a id='NeuroAnalyzer.s_dft-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_dft-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_dft`** &mdash; *Method*.



```julia
s_dft(signal; fs, pad)
```

Return FFT and DFT sample frequencies for a DFT.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `pad::Int64=0`: number of zeros to add

**Returns**

Named tuple containing:

  * `s_fft::Vector{ComplexF64}`
  * `s_sf::Vector{Float64}`

<a id='NeuroAnalyzer.s_msci95-Tuple{Vector{Float64}}' href='#NeuroAnalyzer.s_msci95-Tuple{Vector{Float64}}'>#</a>
**`NeuroAnalyzer.s_msci95`** &mdash; *Method*.



```julia
s_msci95(signal)
```

Calculate mean, std and 95% confidence interval.

**Arguments**

  * `signal::Vector{Float64}`

**Returns**

  * `s_m::Float64`: mean
  * `s_s::Float64`: standard deviation
  * `s_u::Float64`: upper 95% CI
  * `s_l::Float64`: lower 95% CI

<a id='NeuroAnalyzer.s_msci95-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_msci95-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_msci95`** &mdash; *Method*.



```julia
s_msci95(signal; n, method)
```

Calculate mean, std and 95% confidence interval for each the channel.

**Arguments**

  * `signal::AbstractArray`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal method (:normal) or `n`-times boostrapping (:boot)

**Returns**

  * `s_m::Vector{Float64}`: mean
  * `s_s::Vector{Float64}`: standard deviation
  * `s_u::Vector{Float64}`: upper 95% CI
  * `s_l::Vector{Float64}`: lower 95% CI

<a id='NeuroAnalyzer.s2_mean-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroAnalyzer.s2_mean-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroAnalyzer.s2_mean`** &mdash; *Method*.



```julia
s2_mean(signal1, signal2)
```

Calculate mean and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::Vector{Float64}`
  * `signal2:Vector{Float64}`

**Returns**

  * `s_m::Float64`: mean
  * `s_s::Float64`: standard deviation
  * `s_u::Float64`: upper 95% CI
  * `s_l::Float64`: lower 95% CI

<a id='NeuroAnalyzer.s2_difference-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.s2_difference-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.s2_difference`** &mdash; *Method*.



```julia
s2_difference(signal1, signal2; n, method)
```

Calculate mean difference and 95% confidence interval for 2 signals.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:absdiff`: maximum difference (`:absdiff`), integrated area of the squared difference (`:diff2int`)

**Returns**

Named tuple containing:

  * `s_stat::Vector{Float64}`
  * `s_stat_single::Float64`
  * `p::Float64`

<a id='NeuroAnalyzer.s_acov-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_acov-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_acov`** &mdash; *Method*.



s_acov(signal; lag, demean, norm)

Calculate autocovariance.

**Arguments**

  * `signal::AbstractVector`
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean `signal` prior to calculations
  * `norm::Bool=false`: normalize autocovariance

**Returns**

  * `acov::Vector{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroAnalyzer.s2_xcov-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_xcov-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_xcov`** &mdash; *Method*.



s2_xcov(signal1, signal2; lag=1, demean=false, norm=false)

Calculate cross-covariance between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean signal prior to analysis
  * `norm::Bool`: normalize cross-covariance

**Returns**

  * `ccov::Vector{Float64}`
  * `lags::Vector{Int64}`

<a id='NeuroAnalyzer.s_spectrum-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_spectrum-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_spectrum`** &mdash; *Method*.



```julia
s_spectrum(signal; pad)
```

Calculate FFT, amplitudes, powers and phases.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64=0`: number of zeros to add
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `s_fft::Vector{ComplexF64}`
  * `s_amp::Vector{Float64}`
  * `s_pow::Vector{Float64}`
  * `s_pha::Vector{Float64}`

<a id='NeuroAnalyzer.s_total_power-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_total_power-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_total_power`** &mdash; *Method*.



```julia
s_total_power(signal; fs, mt)
```

Calculate total power.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

  * `stp::Float64`: signal total power

<a id='NeuroAnalyzer.s_band_power-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_band_power-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_band_power`** &mdash; *Method*.



```julia
s_band_power(signal; fs, f, mt)
```

Calculate power between `f[1]` and `f[2]`.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

  * `sbp::Float64`: signal band power

<a id='NeuroAnalyzer.s_taper-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_taper-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_taper`** &mdash; *Method*.



```julia
s_taper(signal; taper)
```

Taper the signal.

**Arguments**

  * `signal::AbstractVector`
  * `taper::Union{AbstractVector, Vector{ComplexF64}}`

**Returns**

  * `s_tapered::Vector{Union{Float64, ComplexF64}}`

<a id='NeuroAnalyzer.s_detrend-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_detrend-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_detrend`** &mdash; *Method*.



```julia
s_detrend(signal; type, offset, order, span, fs)
```

Perform piecewise detrending.

**Arguments**

  * `signal::AbstractVector`
  * `type::Symbol`:

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` is subtracted
      * `:loess`: fit and subtract loess approximation
      * `:hp`: use HP filter
  * `offset::Real=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`
  * `fs::Int64=0`: sampling frequency

**Returns**

  * `s_det::Vector{Float64}`

<a id='NeuroAnalyzer.s_demean-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_demean-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_demean`** &mdash; *Method*.



```julia
s_demean(signal)
```

Remove mean value (DC offset).

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `s_demeaned::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_zscore-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_zscore-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_zscore`** &mdash; *Method*.



```julia
s_normalize_zscore(signal)
```

Normalize (by z-score).

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_minmax-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_minmax-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_minmax`** &mdash; *Method*.



```julia
s_normalize_minmax(signal)
```

Normalize in [-1, +1].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::AbstractArray`

<a id='NeuroAnalyzer.s_normalize_max-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_max-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_max`** &mdash; *Method*.



```julia
s_normalize_max(signal)
```

Normalize in [0, +1].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::AbstractArray`

<a id='NeuroAnalyzer.s_normalize_log-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_log-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_log`** &mdash; *Method*.



```julia
s_normalize_log(signal)
```

Normalize using log-transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::AbstractArray`

<a id='NeuroAnalyzer.s_add_noise-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_add_noise-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_add_noise`** &mdash; *Method*.



```julia
s_add_noise(signal)
```

Adds random noise.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_noisy::AbstractArray`

<a id='NeuroAnalyzer.s_resample-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_resample-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_resample`** &mdash; *Method*.



```julia
s_resample(signal; t, new_sr)
```

Resample to `new_sr` sampling frequency.

**Arguments**

  * `signal::AbstractVector`
  * `t::AbstractRange`: time
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_resampled::Vector{Float64}`
  * `t_resampled::AbstractRange`

<a id='NeuroAnalyzer.s_resample-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_resample-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_resample`** &mdash; *Method*.



```julia
s_resample(signal; t, new_sr)
```

Resamples all channels and time vector `t` to `new_sr` sampling frequency.

**Arguments**

  * `signal::AbstractArray`
  * `t::AbstractRange`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_downsampled::Array{Float64, 3}`
  * `t_downsampled::AbstractRange`

<a id='NeuroAnalyzer.s_derivative-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_derivative-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_derivative`** &mdash; *Method*.



```julia
s_derivative(signal)
```

Return derivative of the same length.

**Arguments**

  * `signal::AbstractVector`

<a id='NeuroAnalyzer.s_tconv-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_tconv-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_tconv`** &mdash; *Method*.



```julia
s_tconv(signal; kernel)
```

Performs convolution in the time domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractVector`
  * `kernel::Union{AbstractVector, Vector{ComplexF64}}`

**Returns**

  * `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`

<a id='NeuroAnalyzer.s_filter-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_filter-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_filter`** &mdash; *Method*.



```julia
s_filter(signal; <keyword arguments>)
```

Apply filtering.

**Arguments**

  * `signal::AbstractVector`
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order`
      * `:conv`: convolution
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, `:twopass`), for causal filter use `:onepass`
  * `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
  * `window::Union{AbstractVector, Nothing} - window, required for FIR filter, weighting window for`:mavg`and`:mmed`

**Returns**

  * `s_filtered::Vector{Float64}`

<a id='NeuroAnalyzer.s_psd-Tuple{Vector{Float64}}' href='#NeuroAnalyzer.s_psd-Tuple{Vector{Float64}}'>#</a>
**`NeuroAnalyzer.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm, mt)
```

Calculate power spectrum density.

**Arguments**

  * `signal::Vector{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroAnalyzer.s_psd-Tuple{Matrix{Float64}}' href='#NeuroAnalyzer.s_psd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm, mt)
```

Calculate power spectrum density.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `psd_pow::Matrix{Float64}`
  * `psd_frq::Matrix{Float64}`

<a id='NeuroAnalyzer.s_psd-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_psd-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_psd`** &mdash; *Method*.



```julia
s_psd(signal; fs, norm, mt)
```

Calculate power spectrum density.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `psd_pow::Array{Float64, 3}`
  * `psd_frq::Array{Float64, 3}`

<a id='NeuroAnalyzer.s_stationarity_hilbert-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_stationarity_hilbert-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_stationarity_hilbert`** &mdash; *Method*.



```julia
s_stationarity_hilbert(signal)
```

Calculate phase stationarity using Hilbert transformation.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `phase_stationarity::Vector{Float64}`

<a id='NeuroAnalyzer.s_stationarity_mean-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_stationarity_mean-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_stationarity_mean`** &mdash; *Method*.



```julia
s_stationarity_mean(signal)
```

Calculate mean stationarity.

**Arguments**

  * `signal::AbstractVector`
  * `window::Int64`: time window in samples

**Returns**

  * `mean_stationarity::Vector{Float64}`

<a id='NeuroAnalyzer.s_stationarity_var-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_stationarity_var-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_stationarity_var`** &mdash; *Method*.



```julia
s_stationarity_var(signal)
```

Calculate variance stationarity.

**Arguments**

  * `signal::AbstractVector`
  * `window::Int64`: time window in samples

**Returns**

  * `var_stationarity::Vector{Float64}`

<a id='NeuroAnalyzer.s_trim-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_trim-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_trim`** &mdash; *Method*.



```julia
s_trim(signal; segment)
```

Remove segment from the signal.

**Arguments**

  * `signal::AbstractVector`
  * `segment::Tuple{Int64, Int64}`: segment (from, to) in samples

**Returns**

  * `s_trimmed::Vector{Float64}`

<a id='NeuroAnalyzer.s_trim-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_trim-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_trim`** &mdash; *Method*.



```julia
s_trim(signal; segment)
```

Remove segment from the signal.

**Arguments**

  * `signal::AbstractArray`
  * `segment::Tuple{Int64, Int64}`: segment (from, to) in samples

**Returns**

  * `s_trimmed::Array{Float64}`

<a id='NeuroAnalyzer.s2_mi-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_mi-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_mi`** &mdash; *Method*.



```julia
s2_mi(signal1, signal2)
```

Calculate mutual information between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`

**Returns**

  * `mi::Float64`

<a id='NeuroAnalyzer.s_entropy-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_entropy-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_entropy`** &mdash; *Method*.



```julia
s_entropy(signal)
```

Calculate entropy.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `ent::Float64`
  * `sent::Float64`: Shanon entropy
  * `leent::Float64`: log energy entropy

<a id='NeuroAnalyzer.s_negentropy-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_negentropy-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_negentropy`** &mdash; *Method*.



```julia
s_negentropy(signal)
```

Calculate negentropy.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `ent::Float64`

<a id='NeuroAnalyzer.s_average-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_average-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_average`** &mdash; *Method*.



```julia
s_average(signal)
```

Average all channels.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_averaged::AbstractArray`

<a id='NeuroAnalyzer.s2_average-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.s2_average-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.s2_average`** &mdash; *Method*.



```julia
s2_average(signal1, signal2)
```

Averages `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`

**Returns**

  * `s_averaged::Vector{Float64}`

<a id='NeuroAnalyzer.s2_tcoherence-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_tcoherence-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_tcoherence`** &mdash; *Method*.



```julia
s2_tcoherence(signal1, signal2; pad)
```

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`
  * `pad::Int64=0`: number of zeros to add

**Returns**

Named tuple containing:

  * `c::Vector{Float64}`: coherence
  * `msc::Vector{Float64}`: magnitude-squares coherence
  * `ic::Vector{Float64}`: imaginary part of coherence

<a id='NeuroAnalyzer.s_pca-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_pca-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_pca`** &mdash; *Method*.



```julia
s_pca(signal, n)
```

Calculate `n` first PCs.

**Arguments**

  * `signal::AbstractArray`
  * `n::Int64`: number of PCs

**Returns**

Named tuple containing:

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
  * `pc_m::PCA{Float64}`: PC mean
  * `pca::MultivariateStats.PCA{Float64}`: PC model

<a id='NeuroAnalyzer.s_pca_reconstruct-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_pca_reconstruct-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_pca_reconstruct`** &mdash; *Method*.



```julia
s_pca_reconstruct(signal, pc, pca)
```

Reconstructs `signal` using PCA components.

**Arguments**

  * `signal::AbstractArray`
  * `pc::AbstractArray:`: IC(1)..IC(n) × epoch
  * `pca::MultivariateStats.PCA{Float64}:`: IC(1)..IC(n) × epoch

**Returns**

  * `s_reconstructed::Array{Float64, 3}`

<a id='NeuroAnalyzer.s_fconv-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_fconv-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_fconv`** &mdash; *Method*.



```julia
s_fconv(signal; kernel, norm)
```

Perform convolution in the frequency domain between `signal` and `kernel`.

**Arguments**

  * `signal::AbstractArray`
  * `kernel::Union{AbstractVector, Vector{ComplexF64}}`
  * `pad::Int64=0`: number of zeros to add
  * `norm::Bool=false`: normalize kernel

**Returns**

  * `s_conv::Vector{ComplexF64}`

<a id='NeuroAnalyzer.s_ica-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_ica-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_ica`** &mdash; *Method*.



```julia
s_ica(signal, n, tol, iter, f)
```

Calculate `n` first ICs for `signal`.

**Arguments**

  * `signal::AbstractArray`
  * `n::Int64`: number of PCs
  * `tol::Float64=1.0e-6`: tolerance for ICA
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor (:tanh or :gaus)

**Returns**

Named tuple containing:

  * `ica::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
  * `ica_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch

<a id='NeuroAnalyzer.s_ica_reconstruct-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_ica_reconstruct-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_ica_reconstruct`** &mdash; *Method*.



```julia
s_ica_reconstruct(signal, ica, ica_mw, ic)
```

Reconstructs `signal` using removal of `ic` ICA components.

**Arguments**

  * `signal::AbstractArray`
  * `ica::AbstractArray:`: IC(1)..IC(n) × epoch
  * `ica_mw::AbstractArray:`: IC(1)..IC(n) × epoch
  * `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

**Returns**

  * `s_reconstructed::Array{Float64, 3}`

<a id='NeuroAnalyzer.s_spectrogram-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_spectrogram-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_spectrogram`** &mdash; *Method*.



```julia
s_spectrogram(signal; fs, norm, mt, st, demean)
```

Calculate spectrogram.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=true`: normalize powers to dB
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `st::Bool=false`: if true use short time Fourier transform
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `s_pow::Matrix{Float64}`: powers
  * `s_frq::Vector{Float64}`: frequencies
  * `s_t::Vector{Float64}`: time

<a id='NeuroAnalyzer.s_detect_channel_flat-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_detect_channel_flat-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_detect_channel_flat`** &mdash; *Method*.



```julia
s_detect_channel_flat(signal; w, tol)
```

Detect flat channel.

**Arguments**

  * `signal::AbstractVector`
  * `w::Int64=4`: window width in samples (signal is averaged within `w`-width window)
  * `tol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance

**Returns**

Named tuple containing:

  * `r::Float64`: flat to non-flat segments ratio

<a id='NeuroAnalyzer.s_snr-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_snr-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_snr`** &mdash; *Method*.



```julia
s_snr(signal)
```

Calculate SNR.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `snr::Float64`: SNR

**Source**

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278

<a id='NeuroAnalyzer.s_findpeaks-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_findpeaks-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_findpeaks`** &mdash; *Method*.



```julia
s_findpeaks(signal; d)
```

Find peaks.

**Arguments**

  * `signal::AbstractVector`
  * `d::Int64=32`: distance between peeks in samples

**Returns**

  * `p_idx::Vector{Int64}`

<a id='NeuroAnalyzer.s_wdenoise-Union{Tuple{AbstractVector}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.s_wdenoise-Union{Tuple{AbstractVector}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.s_wdenoise`** &mdash; *Method*.



```julia
s_wdenoise(signal; wt)
```

Perform wavelet denoising.

**Arguments**

  * `signal::AbstractVector`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`

**Returns**

  * `signal_denoised::Vector{Float64}`

<a id='NeuroAnalyzer.s2_ispc-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_ispc-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_ispc`** &mdash; *Method*.



```julia
s2_ispc(signal1, signal2)
```

Calculate ISPC (Inter-Site-Phase Clustering) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`

**Returns**

Named tuple containing:

  * `ispc::Float64`: ISPC value
  * `ispc_angle::Float64`: ISPC angle
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroAnalyzer.s_itpc-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_itpc-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_itpc`** &mdash; *Method*.



```julia
s_itpc(signal; t)
```

Calculate ITPC (Inter-Trial-Phase Clustering) over epochs/trials at time `t`.

**Arguments**

  * `signal::AbstractArray`: one channel over epochs/trials
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc::Float64`: ITPC value
  * `itpcz::Float64`: Rayleigh's ITPC Z value
  * `itpc_angle::Float64`: ITPC angle
  * `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs

<a id='NeuroAnalyzer.s2_pli-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_pli-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_pli`** &mdash; *Method*.



```julia
s2_pli(signal1, signal2)
```

Calculate PLI (Phase-Lag Index) between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`

**Returns**

Named tuple containing:

  * `pli::Float64`: PLI value
  * `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
  * `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase

<a id='NeuroAnalyzer.s2_ged-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.s2_ged-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.s2_ged`** &mdash; *Method*.



```julia
s2_ged(signal1, signal2)
```

Perform generalized eigendecomposition.

**Arguments**

  * `signal1::AbstractArray`: signal to be analyzed
  * `signal2::AbstractArray`: original signal

**Returns**

Named tuple containing:

  * `sged::Matrix{Float64}`
  * `ress::Vector{Float64}`
  * `ress_normalized::Vector{Float64}`: RESS normalized to -1..1

<a id='NeuroAnalyzer.s_frqinst-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_frqinst-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_frqinst`** &mdash; *Method*.



```julia
s_frqinst(signal; fs)
```

Calculate instantaneous frequency.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`

**Returns**

  * `frqinst::Vector{Float64}`

<a id='NeuroAnalyzer.s_hspectrum-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_hspectrum-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_hspectrum`** &mdash; *Method*.



```julia
s_hspectrum(signal; pad=0)
```

Calculate amplitudes, powers and phases using Hilbert transform.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize do dB

**Returns**

Named tuple containing:

  * `h::Vector(ComplexF64}`: Hilbert components
  * `h_amp::Vector{Float64}`
  * `h_pow::Vector{Float64}`
  * `h_pha::Vector{Float64}`

<a id='NeuroAnalyzer.t2f-Tuple{Real}' href='#NeuroAnalyzer.t2f-Tuple{Real}'>#</a>
**`NeuroAnalyzer.t2f`** &mdash; *Method*.



```julia
t2f(t)
```

Convert cycle length in ms `t` to frequency.

**Arguments**

  * `t::Real`: cycle length in ms

**Returns**

  * `f::Float64`: frequency in Hz

<a id='NeuroAnalyzer.f2t-Tuple{Real}' href='#NeuroAnalyzer.f2t-Tuple{Real}'>#</a>
**`NeuroAnalyzer.f2t`** &mdash; *Method*.



```julia
f2t(f)
```

Convert frequency `f` to cycle length in ms.

**Arguments**

  * `f::Real`: frequency in Hz

**Returns**

  * `f::Float64`: cycle length in ms

<a id='NeuroAnalyzer.s_wspectrogram-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_wspectrogram-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_wspectrogram`** &mdash; *Method*.



```julia
s_wspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc, demean)
```

Calculate spectrogram using wavelet convolution.

**Arguments**

  * `signal::AbstractVector`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq*n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq*n) for frq === :lin
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `w_conv::Matrix(ComplexF64}`: convoluted signal
  * `w_powers::Matrix{Float64}`
  * `w_phases::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroAnalyzer.s_fftdenoise-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_fftdenoise-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_fftdenoise`** &mdash; *Method*.



s_fftdenoise(signal; pad, threshold) 

Perform FFT denoising.

**Arguments**

  * `signal::AbstractVector`
  * `pad::Int64=0`: number of zeros to add
  * `threshold::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

**Returns**

Named tuple containing:

  * `s_denoised::Vector{Float64}`
  * `frq_idx::BitVector`: index of components zeroed

<a id='NeuroAnalyzer.s_gfilter-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_gfilter-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_gfilter`** &mdash; *Method*.



```julia
s_gfilter(signal, fs, f, gw)
```

Filter using Gaussian in the frequency domain.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `pad::Int64=0`: number of zeros to add
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

**Returns**

Named tuple containing:

  * `s_f::Vector{Float64}`

<a id='NeuroAnalyzer.s_ghspectrogram-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_ghspectrogram-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_ghspectrogram`** &mdash; *Method*.



```julia
s_ghspectrogram(signal; fs, norm, frq_lim, frq_n, frq, fs, demean)
```

Calculate spectrogram using Gaussian and Hilbert transform.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `gw::Real=5`: Gaussian width in Hz
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `s_pow::Matrix{Float64}`: powers
  * `s_pow::Matrix{Float64}`: phases
  * `frq_list::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.s_tkeo-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_tkeo-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_tkeo`** &mdash; *Method*.



```julia
s_tkeo(signal)
```

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `s_new::Vector{Float64}`

<a id='NeuroAnalyzer.s_mwpsd-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_mwpsd-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_mwpsd`** &mdash; *Method*.



```julia
s_mwpsd(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `signal::AbstractVector`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64=0`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq*n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq*n) for frq === :lin

**Returns**

Named tuple containing:

  * `w_powers::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroAnalyzer.s_mwpsd-Tuple{Matrix{Float64}}' href='#NeuroAnalyzer.s_mwpsd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.s_mwpsd`** &mdash; *Method*.



```julia
s_mwpsd(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `signal::Array{Float64, 2}`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
  * `frq_n::Int64=10`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `fs::Int64`: sampling rate
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq*n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq*n) for frq === :lin

**Returns**

Named tuple containing:

  * `w_powers::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroAnalyzer.a2_cmp-Tuple{Array{<:Real, 3}, Array{<:Real, 3}}' href='#NeuroAnalyzer.a2_cmp-Tuple{Array{<:Real, 3}, Array{<:Real, 3}}'>#</a>
**`NeuroAnalyzer.a2_cmp`** &mdash; *Method*.



```julia
a2_cmp(a1, a2; p, perm_n)
```

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using permutation based statistic.

**Arguments**

  * `a1::Array{<:Real, 3}`: first array
  * `a2::Array{<:Real, 3}`: second array
  * `p::Float64=0.05`: p-value
  * `perm_n::Int64=1000`: number of permutations

**Returns**

Named tuple containing:

  * `zmap::Array{Float64, 3}`: array of Z-values
  * `zmap_b::Array{Float64, 3}`: binarized mask of statistically significant positions

<a id='NeuroAnalyzer.s_fcoherence-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_fcoherence-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_fcoherence`** &mdash; *Method*.



```julia
s_fcoherence(signal; fs, frq_lim)
```

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels.

**Arguments**

  * `signal::AbstractArray`
  * `fs::Int64`: sampling rate
  * `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

**Returns**

Named tuple containing:

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.s2_fcoherence-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.s2_fcoherence-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.s2_fcoherence`** &mdash; *Method*.



```julia
s2_fcoherence(signal1, signal2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)
```

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractArray`
  * `signal2::AbstractArray`
  * `fs::Int64`
  * `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

**Returns**

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.a2_l1-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.a2_l1-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.a2_l1`** &mdash; *Method*.



```julia
a2_l1(a1, a2)
```

Compare two arrays (e.g. two spectrograms), using L1 (Manhattan) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l1::Float64`

<a id='NeuroAnalyzer.a2_l2-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.a2_l2-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.a2_l2`** &mdash; *Method*.



```julia
a2_l2(a1, a2)
```

Compare two arrays (e.g. two spectrograms), using L2 (Euclidean) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l2::Float64`

<a id='NeuroAnalyzer.s_cums-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_cums-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_cums`** &mdash; *Method*.



```julia
s_cums(signal)
```

Calculate cumulative sum.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `signal_cs::Vector{Float64}`

<a id='NeuroAnalyzer.s_cums-Tuple{Array{<:Real, 3}}' href='#NeuroAnalyzer.s_cums-Tuple{Array{<:Real, 3}}'>#</a>
**`NeuroAnalyzer.s_cums`** &mdash; *Method*.



```julia
s_cums(signal)
```

Calculate cumulative sum.

**Arguments**

  * `signal::Array{<:Real, 3}`

**Returns**

  * `signal_cs::Array{Float64, 3}`

<a id='NeuroAnalyzer.s_gfp-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_gfp-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_gfp`** &mdash; *Method*.



```julia
s_gfp(signal)
```

Calculate GFP (Global Field Power).

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `gfp::Float64`

<a id='NeuroAnalyzer.s_gfp_norm-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_gfp_norm-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_gfp_norm`** &mdash; *Method*.



```julia
s_gfp_norm(signal)
```

Calculate signal normalized for GFP (Global Field Power).

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `gfp_norm::Float64`

<a id='NeuroAnalyzer.s2_diss-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_diss-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_diss`** &mdash; *Method*.



```julia
s2_diss(signal1, signal2)
```

Calculate DISS (global dissimilarity) and spatial correlation between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`

**Returns**

Named tuple containing:

  * `diss::Float64`: global dissimilarity
  * `c::Float64`: spatial correlation

<a id='NeuroAnalyzer.generate_morlet_fwhm' href='#NeuroAnalyzer.generate_morlet_fwhm'>#</a>
**`NeuroAnalyzer.generate_morlet_fwhm`** &mdash; *Function*.



```julia
generate_morlet_fwhm(fs, f, t; h)
```

Generate Morlet wavelet using Mike X Cohen formula.

**Arguments**

  * `fs::Int64`: sampling rate
  * `f::Real`: frequency
  * `t::Real=1`: length = -t:1/fs:t
  * `h::Float64=0.25`: full width at half-maximum in seconds (FWHM)

**Returns**

  * `morlet::Vector{ComplexF64}`

**Source**

Cohen MX. A better way to define and describe Morlet wavelets for time-frequency analysis. NeuroImage. 2019 Oct;199:81–6. 

<a id='NeuroAnalyzer.f_nearest-Tuple{Matrix{Tuple{Float64, Float64}}, Tuple{Float64, Float64}}' href='#NeuroAnalyzer.f_nearest-Tuple{Matrix{Tuple{Float64, Float64}}, Tuple{Float64, Float64}}'>#</a>
**`NeuroAnalyzer.f_nearest`** &mdash; *Method*.



```julia
f_nearest(m, pos)
```

Find nearest position tuple `pos` in matrix of positions `m`.

**Arguments**

  * `m::Matrix{Tuple{Float64, Float64}}`
  * `p::Tuple{Float64, Float64}`

**Returns**

  * `pos::Tuple{Int64, Int64}`: row and column in m

<a id='NeuroAnalyzer.s_band_mpower-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_band_mpower-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_band_mpower`** &mdash; *Method*.



```julia
s_band_mpower(signal; fs, f)
```

Calculate mean and maximum band power and its frequency.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds

**Returns**

Named tuple containing:

  * `mbp::Float64`: mean band power [dB]
  * `maxfrq::Float64`: frequency of maximum band power [Hz]
  * `maxbp::Float64`: power at maximum band frequency [dB]

<a id='NeuroAnalyzer.s_rel_psd-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_rel_psd-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_rel_psd`** &mdash; *Method*.



```julia
s_rel_psd(signal; fs, norm, mt, f)
```

Calculate relative power spectrum density.

**Arguments**

  * `signal::AbstractVector`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

**Returns**

Named tuple containing:

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroAnalyzer.s_rel_psd-Tuple{Matrix{Float64}}' href='#NeuroAnalyzer.s_rel_psd-Tuple{Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.s_rel_psd`** &mdash; *Method*.



```julia
s_rel_psd(signal; fs, norm, mt, f)
```

Calculate relative power spectrum density.

**Arguments**

  * `signal::Matrix{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

**Returns**

Named tuple containing:

  * `psd_pow::Vector{Float64}`
  * `psd_frq::Vector{Float64}`

<a id='NeuroAnalyzer.s_wbp-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_wbp-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_wbp`** &mdash; *Method*.



```julia
s_wbp(signal; pad, frq, fs, ncyc, demean)
```

Perform wavelet bandpass filtering.

**Arguments**

  * `signal::AbstractVector`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

  * `signal_new::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_gauss' href='#NeuroAnalyzer.s_normalize_gauss'>#</a>
**`NeuroAnalyzer.s_normalize_gauss`** &mdash; *Function*.



```julia
s_normalize_gauss(signal, dims)
```

Normalize to Gaussian.

**Arguments**

  * `signal::AbstractArray`
  * `dims::Int64=1`: dimension for cumsum()

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_cbp-Tuple{AbstractVector}' href='#NeuroAnalyzer.s_cbp-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.s_cbp`** &mdash; *Method*.



```julia
s_cbp(signal; pad, frq, fs, demean)
```

Perform convolution bandpass filtering.

**Arguments**

  * `signal::AbstractVector`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `fs::Int64`: sampling rate
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

  * `signal_new::Vector{Float64}`

<a id='NeuroAnalyzer.s_specseg-Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}' href='#NeuroAnalyzer.s_specseg-Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroAnalyzer.s_specseg`** &mdash; *Method*.



```julia
s_specseg(sp, sf, st; t, f)
```

Return spectrogram segment.

**Arguments**

  * `sp::Matrix{Float64}`: spectrogram powers
  * `st::Vector{Float64}`: spectrogram time
  * `sf::Vector{Float64}`: spectrogram frequencies
  * `t::Tuple{Real, Real}`: time bounds
  * `f::Tuple{Real, Real}`: frequency bounds

**Returns**

Named tuple containing:

  * `seg_pow::Matrix{Float64}`: powers
  * `seg_shape::Shape{Real, Int64}`: shape for plotting
  * `t_idx::Tuple{Real, Real}`: time indices
  * `f_idx::Tuple{Real, Real}`: frequency indices

<a id='NeuroAnalyzer.s_specseg-Tuple{AbstractArray, AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s_specseg-Tuple{AbstractArray, AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s_specseg`** &mdash; *Method*.



```julia
s_specseg(sp, sf, st; t, f)
```

Return spectrogram segment.

**Arguments**

  * `sp::AbstractArray`: spectrogram powers
  * `st::AbstractVector`: spectrogram time
  * `sf::AbstractVector`: spectrogram frequencies
  * `t::Tuple{Real, Real}`: time bounds
  * `f::Tuple{Real, Real}`: frequency bounds

**Returns**

Named tuple containing:

  * `seg_pow::Array{Float64, 3}`: segment of powers
  * `seg_shape::Shape{Real, Int64}`: segment coordinates (shape for plotting)
  * `t_idx::Tuple{Real, Real}`: time indices
  * `f_idx::Tuple{Real, Real}`: frequency indices

<a id='NeuroAnalyzer.s_denoise_wien-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_denoise_wien-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_denoise_wien`** &mdash; *Method*.



```julia
s_denoise_wien(signal)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `signal_new::Vector{Float64}`

<a id='NeuroAnalyzer.s2_cps-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_cps-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_cps`** &mdash; *Method*.



```julia
s2_cps(signal1, signal2; fs, norm)
```

Calculate cross power spectrum between `signal1` and `signal2`.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize do dB

**Returns**

Named tuple containing:

  * `cps_pw::Vector{Float64}`: cross power spectrum power
  * `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
  * `cps_fq::Vector{Float64}`: cross power spectrum frequencies

<a id='NeuroAnalyzer.s2_phdiff-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_phdiff-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_phdiff`** &mdash; *Method*.



```julia
s_phdiff(signal1, signal2; pad, h)
```

Calculate phase difference between signals.

**Arguments**

  * `signal1::AbstractVector`
  * `signal2::AbstractVector`
  * `pad::Int64=0`: number of zeros to add
  * `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

**Returns**

Named tuple containing:

  * `ph_diff::Vector{Float64}`: phase differences in radians

<a id='NeuroAnalyzer.s_normalize_log10-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_log10-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_log10`** &mdash; *Method*.



```julia
s_normalize_log10(signal)
```

Normalize using log10-transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_neglog-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_neglog-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_neglog`** &mdash; *Method*.



```julia
s_normalize_neglog(signal)
```

Normalize to using -log-transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_neglog10-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_neglog10-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_neglog10`** &mdash; *Method*.



```julia
s_normalize_neglog10(signal)
```

Normalize using -log10-transformation.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_neg-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_neg-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_neg`** &mdash; *Method*.



```julia
s_normalize_neg(signal)
```

Normalize in [0, -∞].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_pos-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_pos-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_pos`** &mdash; *Method*.



```julia
s_normalize_pos(signal)
```

Normalize in [0, +∞].

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize_perc-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_perc-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_perc`** &mdash; *Method*.



```julia
s_normalize_perc(signal)
```

Normalize in percentages.

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_normalize-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize`** &mdash; *Method*.



```julia
s_normalize(signal; method)
```

Normalize.

**Arguments**

  * `signal::AbstractArray`
  * `method::Symbol`: :zscore, :minmax, :max, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :none

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_phases' href='#NeuroAnalyzer.s_phases'>#</a>
**`NeuroAnalyzer.s_phases`** &mdash; *Function*.



```julia
s_phases(signal; pad)
```

Calculate phases.

**Arguments**

  * `signal::AbstractArray`
  * `pad::Int64=0`: number of zeros to add

**Returns**

  * `phases::Vector{Float64}`

<a id='NeuroAnalyzer.s_cwtspectrogram-Union{Tuple{AbstractVector}, Tuple{T}} where T<:CWT' href='#NeuroAnalyzer.s_cwtspectrogram-Union{Tuple{AbstractVector}, Tuple{T}} where T<:CWT'>#</a>
**`NeuroAnalyzer.s_cwtspectrogram`** &mdash; *Method*.



```julia
s_cwtspectrogram(signal; wt, pad, norm, frq_lim, fs, demean)
```

Calculate spectrogram using continuous wavelet transformation (CWT).

**Arguments**

  * `signal::AbstractVector`
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `demean::Bool`=true: demean signal prior to analysis

**Returns**

Named tuple containing:

  * `h_powers::Matrix{Float64}`
  * `frq_list::Vector{Float64}`

<a id='NeuroAnalyzer.s_dwt-Union{Tuple{AbstractVector}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.s_dwt-Union{Tuple{AbstractVector}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.s_dwt`** &mdash; *Method*.



```julia
s_dwt(signal; wt, type, l)
```

Perform discrete wavelet transformation (DWT).

**Arguments**

  * `signal::AbstractVector`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: Stationary Wavelet Transforms (:sdwt) or Autocorrelation Wavelet Transforms (:acdwt)
  * `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

**Returns**

  * `dwt_c::Array{Float64, 2}`: DWT coefficients cAl, cD1, ..., cDl (by rows)

<a id='NeuroAnalyzer.s_idwt-Union{Tuple{AbstractArray}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.s_idwt-Union{Tuple{AbstractArray}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.s_idwt`** &mdash; *Method*.



```julia
s_idwt(dwt_coefs; wt, type)
```

Perform inverse discrete wavelet transformation (iDWT) of the `dwt_coefs`.

**Arguments**

  * `dwt_coefs::AbstractArray`: DWT coefficients cAl, cD1, ..., cDl (by rows)
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: Stationary Wavelet Transforms (:sdwt) or Autocorrelation Wavelet Transforms (:acdwt)

**Returns**

  * `signal::Vector{Float64}`: reconstructed signal

<a id='NeuroAnalyzer.s_normalize_invroot-Tuple{AbstractArray}' href='#NeuroAnalyzer.s_normalize_invroot-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.s_normalize_invroot`** &mdash; *Method*.



```julia
s_normalize_invroot(signal)
```

Normalize in inverse root (1/sqrt(x)).

**Arguments**

  * `signal::AbstractArray`

**Returns**

  * `s_normalized::Vector{Float64}`

<a id='NeuroAnalyzer.s_cwt-Union{Tuple{AbstractVector}, Tuple{T}} where T<:CWT' href='#NeuroAnalyzer.s_cwt-Union{Tuple{AbstractVector}, Tuple{T}} where T<:CWT'>#</a>
**`NeuroAnalyzer.s_cwt`** &mdash; *Method*.



```julia
s_cwt(signal; wt, type, l)
```

Perform continuous wavelet transformation (CWT).

**Arguments**

  * `signal::AbstractVector`
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

  * `cwt_c::Array{Float64, 2}`: CWT coefficients (by rows)

<a id='NeuroAnalyzer.s_icwt-Union{Tuple{AbstractArray}, Tuple{T}} where T<:CWT' href='#NeuroAnalyzer.s_icwt-Union{Tuple{AbstractArray}, Tuple{T}} where T<:CWT'>#</a>
**`NeuroAnalyzer.s_icwt`** &mdash; *Method*.



```julia
s_icwt(dwt_coefs; wt, type)
```

Perform inverse continuous wavelet transformation (iCWT) of the `dwt_coefs`.

**Arguments**

  * `cwt_coefs::AbstractArray`: CWT coefficients (by rows)
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `type::Symbol=df`: inverse style type: NaiveDelta (:nd), PenroseDelta (:pd) or DualFrames (:df)

**Returns**

  * `signal::Vector{Float64}`: reconstructed signal

<a id='NeuroAnalyzer.t2s-Tuple{Real, Int64}' href='#NeuroAnalyzer.t2s-Tuple{Real, Int64}'>#</a>
**`NeuroAnalyzer.t2s`** &mdash; *Method*.



```julia
t2s(t, fs)
```

Convert time to sample number.

**Arguments**

  * `t::Real`: time in s
  * `fs::Int64`: sampling rate

**Returns**

  * `s::Int6464`: sample number

<a id='NeuroAnalyzer.s2t-Tuple{Int64, Int64}' href='#NeuroAnalyzer.s2t-Tuple{Int64, Int64}'>#</a>
**`NeuroAnalyzer.s2t`** &mdash; *Method*.



```julia
s2t(s, fs)
```

Convert sample number to time.

**Arguments**

  * `t::Int64`: sample number
  * `fs::Int64`: sampling rate

**Returns**

  * `t::Float64`: time in s


<a id='statistics'></a>

<a id='statistics-1'></a>

## statistics

<a id='NeuroAnalyzer.hildebrand_rule-Tuple{AbstractVector}' href='#NeuroAnalyzer.hildebrand_rule-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.hildebrand_rule`** &mdash; *Method*.



```julia
hildebrand_rule(x)
```

Calculate Hildebrand rule for vector `x`. If H < 0.2 then the vector `x` is symmetrical.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `h::Float64`

<a id='NeuroAnalyzer.jaccard_similarity-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.jaccard_similarity-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.jaccard_similarity`** &mdash; *Method*.



```julia
jaccard_similarity(x, y)
```

Calculate Jaccard similarity between two vectors `x` and `y`.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `j::Float64`

<a id='NeuroAnalyzer.z_score-Tuple{AbstractVector}' href='#NeuroAnalyzer.z_score-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.z_score`** &mdash; *Method*.



```julia
z_score(x)
```

Calculate Z-scores for each value of the vector `x`.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `z_score::Vector{Float64}`

<a id='NeuroAnalyzer.k_categories-Tuple{Int64}' href='#NeuroAnalyzer.k_categories-Tuple{Int64}'>#</a>
**`NeuroAnalyzer.k_categories`** &mdash; *Method*.



```julia
k_categories(n)
```

Calculate number of categories for a given sample size `n`.

**Arguments**

  * `n::Int64`

**Returns**

Named tuple containing:

  * `k1::Float64`: sqrt(n)
  * `k2::Float64`: 1 + 3.222 * log10(n)

<a id='NeuroAnalyzer.effsize-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.effsize-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.effsize`** &mdash; *Method*.



```julia
effsize(x1, x2)
```

Calculate Cohen's d and Hedges g effect sizes.

**Arguments**

  * `x1::AbstractVector`
  * `x2::AbstractVector`

**Returns**

Named tuple containing:

  * `d::Float64`: Cohen's d
  * `g::Float64`: Hedges g

<a id='NeuroAnalyzer.infcrit-Tuple{T} where T<:StatsModels.TableRegressionModel' href='#NeuroAnalyzer.infcrit-Tuple{T} where T<:StatsModels.TableRegressionModel'>#</a>
**`NeuroAnalyzer.infcrit`** &mdash; *Method*.



```julia
infcrit(m)
```

Calculate Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression model `m`.

**Arguments**

  * `m::StatsModels.TableRegressionModel`

**Returns**

Named tuple containing:

  * `aic::Float64`
  * `bic::Float64`

<a id='NeuroAnalyzer.grubbs-Tuple{AbstractVector}' href='#NeuroAnalyzer.grubbs-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.grubbs`** &mdash; *Method*.



```julia
grubbs(x; alpha, t)
```

Perform Grubbs test for outlier.

**Arguments**

  * `x::AbstractVector`
  * `alpha::Float64=0.95`
  * `t::Int64=0`: test type: -1 test whether the minimum value is an outlier; 0 two-sided test; 1 test whether the maximum value is an outlier

**Returns**

Named tuple containing:

  * `g::Bool`: true: outlier exists, false: there is no outlier

<a id='NeuroAnalyzer.outlier_detect-Tuple{AbstractVector}' href='#NeuroAnalyzer.outlier_detect-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.outlier_detect`** &mdash; *Method*.



```julia
outlier_detect(x; method)
```

Detect outliers.

**Arguments**

  * `x::AbstractVector`
  * `method::Symbol=iqr`: methods: `:iqr` (interquartile range), `:z` (z-score) or `:g` (Grubbs test)

**Returns**

  * `o::Vector{Bool}`: index of outliers

<a id='NeuroAnalyzer.seg_mean-Tuple{AbstractArray}' href='#NeuroAnalyzer.seg_mean-Tuple{AbstractArray}'>#</a>
**`NeuroAnalyzer.seg_mean`** &mdash; *Method*.



```julia
seg_mean(seg)
```

Calculate mean of a segment (e.g. spectrogram).

**Arguments**

  * `seg::AbstractArray`

**Returns**

  * `seg::Vector{Float64}`: averaged segment

<a id='NeuroAnalyzer.seg2_mean-Tuple{AbstractArray, AbstractArray}' href='#NeuroAnalyzer.seg2_mean-Tuple{AbstractArray, AbstractArray}'>#</a>
**`NeuroAnalyzer.seg2_mean`** &mdash; *Method*.



```julia
seg2_mean(seg1, seg2)
```

Calculate mean of two segments (e.g. spectrograms).

**Arguments**

  * `seg1::AbstractArray`
  * `seg2::AbstractArray`

**Returns**

Named tuple containing:

  * `seg1::Vector{Float64}`: averaged segment 1
  * `seg2::Vector{Float64}`: averaged segment 2

<a id='NeuroAnalyzer.binom_prob-Tuple{Float64, Int64, Int64}' href='#NeuroAnalyzer.binom_prob-Tuple{Float64, Int64, Int64}'>#</a>
**`NeuroAnalyzer.binom_prob`** &mdash; *Method*.



```julia
binom_prob(p, r, n)
```

Calculate probability of exactly `r` successes in `n` trials.

**Arguments**

  * `p::Float64`: proportion of successes
  * `r::Int64`: number of successes
  * `n::Int64`: number of trials

**Returns**

  * `binomp::Float64`: probability

<a id='NeuroAnalyzer.binom_stat-Tuple{Float64, Int64}' href='#NeuroAnalyzer.binom_stat-Tuple{Float64, Int64}'>#</a>
**`NeuroAnalyzer.binom_stat`** &mdash; *Method*.



```julia
binom_stat(p, n)
```

Calculate mean and standard deviation for probability `p`.

**Arguments**

  * `p::Float64`: proportion of successes
  * `n::Int64`: number of trials

**Returns**

  * `mean::Float64`
  * `std::Float64`

<a id='NeuroAnalyzer.cvar_mean-Tuple{AbstractVector}' href='#NeuroAnalyzer.cvar_mean-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.cvar_mean`** &mdash; *Method*.



```julia
cvar_mean(x)
```

Calculate coefficient of variation for a mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `cvar::Float64`

<a id='NeuroAnalyzer.cvar_median-Tuple{AbstractVector}' href='#NeuroAnalyzer.cvar_median-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.cvar_median`** &mdash; *Method*.



```julia
cvar_median(x)
```

Calculate coefficient of variation for a median.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `cvar::Float64`

<a id='NeuroAnalyzer.cvar-Tuple{Real, Real}' href='#NeuroAnalyzer.cvar-Tuple{Real, Real}'>#</a>
**`NeuroAnalyzer.cvar`** &mdash; *Method*.



```julia
cvar(se, s)
```

Calculate coefficient of variation for statistic `s`.

**Arguments**

  * `se::Real`: standard error
  * `s::Real`: statistics, e.g. mean value

**Returns**

  * `cvar::Float64`

<a id='NeuroAnalyzer.effsize-Tuple{Float64, Float64}' href='#NeuroAnalyzer.effsize-Tuple{Float64, Float64}'>#</a>
**`NeuroAnalyzer.effsize`** &mdash; *Method*.



```julia
effsize(p1, p2)
```

Calculate effect size for two proportions `p1` and `p2`.

**Arguments**

  * `p1::Float64`: 1st proportion, e.g. 0.7
  * `p2::Float64`: 2nd proportion, e.g. 0.3

**Returns**

  * `e::Float64`

<a id='NeuroAnalyzer.meang-Tuple{AbstractVector}' href='#NeuroAnalyzer.meang-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.meang`** &mdash; *Method*.



```julia
meang(x)
```

Calculate geometric mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.meanh-Tuple{AbstractVector}' href='#NeuroAnalyzer.meanh-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.meanh`** &mdash; *Method*.



```julia
meanh(x)
```

Calculate harmonic mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.meanw-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.meanw-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.meanw`** &mdash; *Method*.



```julia
meanw(x, w)
```

Calculate weighted mean.

**Arguments**

  * `x::AbstractVector`
  * `w::AbstractVector`: weights

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.moe-Tuple{Int64}' href='#NeuroAnalyzer.moe-Tuple{Int64}'>#</a>
**`NeuroAnalyzer.moe`** &mdash; *Method*.



```julia
moe(n)
```

Calculate margin of error for given sample size `n`.

**Arguments**

  * `n::Int64`

**Returns**

  * `moe::Float64`

<a id='NeuroAnalyzer.rng-Tuple{AbstractVector}' href='#NeuroAnalyzer.rng-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.rng`** &mdash; *Method*.



```julia
rng(x)
```

Calculate range.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `r::Float64`

<a id='NeuroAnalyzer.se-Tuple{AbstractVector}' href='#NeuroAnalyzer.se-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.se`** &mdash; *Method*.



```julia
se(x)
```

Calculate standard error.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `se::Float64`

<a id='NeuroAnalyzer.pred_int-Tuple{Int64}' href='#NeuroAnalyzer.pred_int-Tuple{Int64}'>#</a>
**`NeuroAnalyzer.pred_int`** &mdash; *Method*.



```julia
pred_int(n)
```

Calculates the prediction interval (95% CI adjusted for sample size)

**Arguments**

  * `n::Int64`: sample size

**Returns**

  * `pred_int::Float64`

<a id='NeuroAnalyzer.sem_diff-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.sem_diff-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.sem_diff`** &mdash; *Method*.



```julia
sem_diff(x::AbstractVector, y::AbstractVector)
```

Calculate SEM (standard error of the mean) for the difference of two means.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `sem::Float64`

<a id='NeuroAnalyzer.prank-Tuple{AbstractVector}' href='#NeuroAnalyzer.prank-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.prank`** &mdash; *Method*.



```julia
prank(x)
```

Calculate percentile rank.

**Arguments**

  * `x::AbstractVector`: the vector to analyze

**Returns**

  * `prnk::Vector{Float64}`

<a id='NeuroAnalyzer.linreg-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.linreg-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.linreg`** &mdash; *Method*.



```julia
linreg(x, y)
```

Linear regression between `x` and `y`.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Notes**

To predict, use: `new*x = DataFrame(x = [3.5, 7]); predict(lr, new*x)

**Returns**

Named tuple containing:

  * `lr::StatsModels.TableRegressionModel`: model
  * `radj::Flpoat64`: R^2
  * `c::Vector{Float64}`: coefficients
  * `se::Vector{Float64}`: standard error for coefficients
  * `aic::Float64`:: Akaike’s Information Criterion (AIC)
  * `bic::Float64`:: Bayesian Information Criterion (BIC)
  * `lf::Vector{Float64}`: linear fit (plot(x, lf))

<a id='NeuroAnalyzer.s2_cmp-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_cmp-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_cmp`** &mdash; *Method*.



```julia
s2_cmp(seg1, seg2, paired, alpha, type, exact)
```

Compare two vectors; Kruskall-Wallis test is used first, next t-test (paired on non-paired) or non-parametric test (paired: Wilcoxon signed rank, non-paired: Mann-Whitney U test) is applied.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `paired::Bool`
  * `alpha::Float64=0.05`: confidence level
  * `type::Symbol=:auto`: choose test automatically (:auto), parametric (:p) or non-parametric (:np)
  * `exact::Bool=false`: if true, use exact Wilcoxon test

**Returns**

Named tuple containing:

  * `t`: test results
  * `ts::Tuple{Float64, String}`: test statistics
  * `tc::Tuple{Float64, Float64}`: test statistics confidence interval
  * `df::Int64`: degrees of freedom
  * `p::Float64`: p-value

<a id='NeuroAnalyzer.s2_cor-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.s2_cor-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.s2_cor`** &mdash; *Method*.



```julia
s2_cor(seg1, seg2)
```

Calculate correlation between two vectors.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

Named tuple containing:

  * `t::CorrelationTest{Float64}`
  * `r::Float64`: correlation coefficient
  * `rc::Tuple{Float64, Float64}`: correlation coefficient confidence interval
  * `tt::Tuple{Float64, String}`: t-statistics
  * `df::Int64`: degrees of freedom
  * `p::Float64`: p-value


<a id='EEG-I/O'></a>

<a id='EEG-I/O-1'></a>

## EEG I/O

<a id='NeuroAnalyzer.eeg_import-Tuple{String}' href='#NeuroAnalyzer.eeg_import-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import`** &mdash; *Method*.



```julia
eeg_import(file_name; detect_type)
```

Load EEG file and return `NeuroAnalyzer.EEG` object. Supported formats:

  * EDF/EDF+
  * BDF/BDF+
  * BrainVision

This is a meta-function that triggers appropriate `eeg_import_*()` function. File format is detected based on file extension (.edf|.bdf|.vhdr).

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `eeg:EEG`

<a id='NeuroAnalyzer.eeg_import_edf-Tuple{String}' href='#NeuroAnalyzer.eeg_import_edf-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_edf`** &mdash; *Method*.



```julia
eeg_import_edf(file_name; detect_type)
```

Load EDF/EDF+ file and return `NeuroAnalyzer.EEG` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `eeg:EEG`

**Notes**

  * sampling_rate = n.samples / data.record.duration
  * gain = (physical*maximum - physical*minimum) / (digital*maximum - digital*minimum)
  * value = (value - digital*minimum ) * gain + physical*minimum

**Source**

1. Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3.
2. Kemp B, Olivan J. European data format ‘plus’ (EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003;114:1755–61.
3. https://www.edfplus.info/specs/

<a id='NeuroAnalyzer.eeg_import_ced-Tuple{String}' href='#NeuroAnalyzer.eeg_import_ced-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_ced`** &mdash; *Method*.



```julia
eeg_import_ced(file_name)
```

Load electrode positions from CED file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_locs-Tuple{String}' href='#NeuroAnalyzer.eeg_import_locs-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_locs`** &mdash; *Method*.



```julia
eeg_import_locs(file_name)
```

Load electrode positions from LOCS file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_elc-Tuple{String}' href='#NeuroAnalyzer.eeg_import_elc-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_elc`** &mdash; *Method*.



```julia
eeg_import_elc(file_name)
```

Load electrode positions from ELC file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_tsv-Tuple{String}' href='#NeuroAnalyzer.eeg_import_tsv-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_tsv`** &mdash; *Method*.



```julia
eeg_import_tsv(file_name)
```

Load electrode positions from TSV file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_sfp-Tuple{String}' href='#NeuroAnalyzer.eeg_import_sfp-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_sfp`** &mdash; *Method*.



```julia
eeg_import_sfp(file_name)
```

Load electrode positions from SFP file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_csd-Tuple{String}' href='#NeuroAnalyzer.eeg_import_csd-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_csd`** &mdash; *Method*.



```julia
eeg_import_csd(file_name)
```

Load electrode positions from CSD file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_load_electrodes-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_load_electrodes-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_load_electrodes`** &mdash; *Method*.



```julia
eeg_load_electrodes(eeg; file_name)
```

Load electrode positions from `file_name` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Accepted formats:

  * CED
  * LOCS
  * ELC
  * TSV
  * SFP
  * CSD

Electrode locations:

  * `loc_theta`       planar polar angle
  * `loc_radius`      planar polar radius
  * `loc_x`           spherical Cartesian x
  * `loc_y`           spherical Cartesian y
  * `loc_z`           spherical Cartesian z
  * `loc_radius_sph`  spherical radius
  * `loc_theta_sph`   spherical horizontal angle
  * `loc_phi_sph`     spherical azimuth angle

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `file_name::String`

**Returns**

  * `eeg:EEG`

<a id='NeuroAnalyzer.eeg_load_electrodes!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_load_electrodes!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_load_electrodes!`** &mdash; *Method*.



```julia
eeg_load_electrodes!(eeg; file_name)
```

Load electrode positions from `file_name` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Accepted formats:

  * CED
  * LOCS
  * ELC
  * TSV
  * SFP

Electrode locations:

  * `loc_theta`       planar polar angle
  * `loc_radius`      planar polar radius
  * `loc_x`           spherical Cartesian x
  * `loc_y`           spherical Cartesian y
  * `loc_z`           spherical Cartesian z
  * `loc_radius_sph`  spherical radius
  * `loc_theta_sph`   spherical horizontal angle
  * `loc_phi_sph`     spherical azimuth angle

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `file_name::String`

<a id='NeuroAnalyzer.eeg_save-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_save-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_save`** &mdash; *Method*.



```julia
eeg_save(eeg; file_name, overwrite)
```

Save `eeg` to `file_name` file (HDF5-based).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `file_name::String`: file name
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`

<a id='NeuroAnalyzer.eeg_load-Tuple{String}' href='#NeuroAnalyzer.eeg_load-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_load`** &mdash; *Method*.



```julia
eeg_load(file_name)
```

Load `eeg` from `file_name` file (HDF5-based).

**Arguments**

  * `file_name::String`: file name

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_export_csv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_export_csv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_export_csv`** &mdash; *Method*.



```julia
eeg_export_csv(eeg; file_name, header, components, markers, overwrite)
```

Export EEG data as CSV.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `file_name::String`
  * `header::Bool=false`: export header
  * `components::Bool=false`: export components
  * `markers::Bool=false`: export markers
  * `locs::Bool=false`: export locations
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`

<a id='NeuroAnalyzer.eeg_save_electrodes-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_save_electrodes-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_save_electrodes`** &mdash; *Method*.



```julia
eeg_save_electrodes(eeg; file_name, overwrite)
```

Export EEG channel locations data, format is based on `file_name` extension (.ced, .locs or .tsv)

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `file_name::String`
  * `overwrite::Bool=false`

<a id='NeuroAnalyzer.eeg_save_electrodes-Tuple{DataFrame}' href='#NeuroAnalyzer.eeg_save_electrodes-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.eeg_save_electrodes`** &mdash; *Method*.



```julia
eeg_save_electrodes(locs; file_name, overwrite)
```

Export channel locations, format is based on `file_name` extension (.ced, .locs, .tsv)

**Arguments**

  * `locs::DataFrame`
  * `file_name::String`
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`

<a id='NeuroAnalyzer.eeg_add_electrodes-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_electrodes-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_electrodes`** &mdash; *Method*.



```julia
eeg_add_electrodes(eeg; locs)
```

Add electrode positions from `locs`. 

Electrode locations:

  * `channel`         channel number
  * `labels`          channel label
  * `loc_theta`       planar polar angle
  * `loc_radius`      planar polar radius
  * `loc_x`           spherical Cartesian x
  * `loc_y`           spherical Cartesian y
  * `loc_z`           spherical Cartesian z
  * `loc_radius_sph`  spherical radius
  * `loc_theta_sph`   spherical horizontal angle
  * `loc_phi_sph`     spherical azimuth angle

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `locs::DataFrame`

**Returns**

  * `eeg:EEG`

<a id='NeuroAnalyzer.eeg_add_electrodes!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_electrodes!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_electrodes!`** &mdash; *Method*.



```julia
eeg_add_electrodes!(eeg; locs)
```

Load electrode positions from `locs` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Electrode locations:

  * `channel`         channel number
  * `labels`          channel label
  * `loc_theta`       planar polar angle
  * `loc_radius`      planar polar radius
  * `loc_x`           spherical Cartesian x
  * `loc_y`           spherical Cartesian y
  * `loc_z`           spherical Cartesian z
  * `loc_radius_sph`  spherical radius
  * `loc_theta_sph`   spherical horizontal angle
  * `loc_phi_sph`     spherical azimuth angle

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_import_bdf-Tuple{String}' href='#NeuroAnalyzer.eeg_import_bdf-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_bdf`** &mdash; *Method*.



```julia
eeg_import_bdf(file_name; detect_type)
```

Load BDF/BDF+ file and return `NeuroAnalyzer.EEG` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `eeg:EEG`

**Notes**

  * sampling_rate = n.samples / data.record.duration
  * gain = (physical*maximum - physical*minimum) / (digital*maximum - digital*minimum)
  * value = (value - digital*minimum ) * gain + physical*minimum

**Source**

https://www.biosemi.com/faq/file_format.htm

<a id='NeuroAnalyzer.eeg_import_digitrack-Tuple{String}' href='#NeuroAnalyzer.eeg_import_digitrack-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_digitrack`** &mdash; *Method*.



```julia
eeg_import_digitrack(file_name; detect_type)
```

Load Digitrack ASCII file and return `NeuroAnalyzer.EEG` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `eeg:EEG`

**Notes**

<a id='NeuroAnalyzer.eeg_import_bv-Tuple{String}' href='#NeuroAnalyzer.eeg_import_bv-Tuple{String}'>#</a>
**`NeuroAnalyzer.eeg_import_bv`** &mdash; *Method*.



```julia
eeg_import_bv(file_name; detect_type)
```

Load BrainVision BVCDF file and return `NeuroAnalyzer.EEG` object. At least two files are required: .vhdr (header) and .eeg (signal data). If available, markers are loaded from .vmrk file.

**Arguments**

  * `file_name::String`: name of the file to load, should point to .vhdr file.
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `eeg:EEG`


<a id='EEG-edit'></a>

<a id='EEG-edit-1'></a>

## EEG edit

<a id='NeuroAnalyzer.eeg_add_component-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_component-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_component`** &mdash; *Method*.



```julia
eeg_add_component(eeg; c, v)
```

Add component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name
  * `v::Any`: component value

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_component!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_component!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_component!`** &mdash; *Method*.



```julia
eeg_add_component!(eeg; c, v)
```

Add component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name
  * `v::Any`: component value

<a id='NeuroAnalyzer.eeg_list_components-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_list_components-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_list_components`** &mdash; *Method*.



```julia
eeg_list_components(eeg)
```

List component names.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `components::Vector{Symbol}`

<a id='NeuroAnalyzer.eeg_extract_component-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_extract_component-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_extract_component`** &mdash; *Method*.



```julia
eeg_extract_component(eeg, c)
```

Extract component values.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name

**Returns**

  * `component::Any`

<a id='NeuroAnalyzer.eeg_delete_component-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_component-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_component`** &mdash; *Method*.



```julia
eeg_delete_component(eeg; c)
```

Delete component. 

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_component!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_component!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_component!`** &mdash; *Method*.



```julia
eeg_delete_component!(eeg; c)
```

Delete component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name

<a id='NeuroAnalyzer.eeg_reset_components-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reset_components-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reset_components`** &mdash; *Method*.



```julia
eeg_reset_components(eeg)
```

Remove all components.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reset_components!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reset_components!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reset_components!`** &mdash; *Method*.



```julia
eeg_reset_components!(eeg)
```

Remove all components.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_component_idx-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_component_idx-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_component_idx`** &mdash; *Method*.



```julia
eeg_component_idx(eeg, c)
```

Return component index.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name

**Return**

  * `c_idx::Int64`

<a id='NeuroAnalyzer.eeg_component_type-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_component_type-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_component_type`** &mdash; *Method*.



```julia
eeg_component_type(eeg, c)
```

Return component data type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Symbol`: component name

**Return**

  * `c_type::DataType`

<a id='NeuroAnalyzer.eeg_rename_component-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_rename_component-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_rename_component`** &mdash; *Method*.



```julia
eeg_rename_component(eeg, c_old, c_new)
```

Rename component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

**Return**

  * `eeg_new:EEG`

<a id='NeuroAnalyzer.eeg_rename_component!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_rename_component!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_rename_component!`** &mdash; *Method*.



```julia
eeg_rename_component!(eeg, c_old, c_new)
```

Rename component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

<a id='NeuroAnalyzer.eeg_delete_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_channel`** &mdash; *Method*.



```julia
eeg_delete_channel(eeg; channel)
```

Delete EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_channel!`** &mdash; *Method*.



```julia
eeg_delete_channel!(eeg; channel)
```

Delete EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed

<a id='NeuroAnalyzer.eeg_keep_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_channel`** &mdash; *Method*.



```julia
eeg_keep_channel(eeg; channel)
```

Keep EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_keep_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_channel!`** &mdash; *Method*.



```julia
eeg_keep_channel!(eeg; channel)
```

Keep EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep

<a id='NeuroAnalyzer.eeg_get_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_get_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_get_channel`** &mdash; *Method*.



```julia
eeg_get_channel(eeg; channel)
```

Return EEG channel number (if provided by name) or name (if provided by number).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name

**Returns**

  * `channel_idx::Union{Int64, String}`: channel number or name

<a id='NeuroAnalyzer.eeg_rename_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_rename_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_rename_channel`** &mdash; *Method*.



```julia
eeg_rename_channel(eeg; channel, name)
```

Rename EEG channel.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name
  * `name::String`: new name

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_rename_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_rename_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_rename_channel!`** &mdash; *Method*.



```julia
eeg_rename_channel!(eeg; channel, name)
```

Rename EEG channel.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name
  * `name::String`: new name

<a id='NeuroAnalyzer.eeg_extract_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_extract_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_extract_channel`** &mdash; *Method*.



```julia
eeg_extract_channel(eeg; channel)
```

Extract EEG channel data.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name

**Returns**

  * `channel::Vector{Float64}`

<a id='NeuroAnalyzer.eeg_history-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_history-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_history`** &mdash; *Method*.



```julia
eeg_history(eeg)
```

Show EEG processing history.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_labels-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_labels-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_labels`** &mdash; *Method*.



```julia
eeg_labels(eeg)
```

Return EEG channel labels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroAnalyzer.eeg_sr-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_sr-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_sr`** &mdash; *Method*.



```julia
eeg_sr(eeg)
```

Return EEG sampling rate.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `sr::Int64`

<a id='NeuroAnalyzer.eeg_channel_n-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_channel_n-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_channel_n`** &mdash; *Method*.



```julia
eeg_channel_n(eeg; type=:eeg)
```

Return number of EEG channels of `type`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Vector{Symbol}=:all`: channel type: `:all`, `:eeg`, `:meg`, `:ecg`, `:eog`, `:emg`, `:ref`, `:mrk`

**Returns**

  * `channel_n::Int64`

<a id='NeuroAnalyzer.eeg_epoch_n-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epoch_n-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epoch_n`** &mdash; *Method*.



```julia
eeg_epoch_n(eeg)
```

Return number of EEG epochs.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `epoch_n::Int64`

<a id='NeuroAnalyzer.eeg_signal_len-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_signal_len-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_signal_len`** &mdash; *Method*.



```julia
eeg_signal_len(eeg)
```

Return length of EEG signal.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `signal_len::Int64`

<a id='NeuroAnalyzer.eeg_epoch_len-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epoch_len-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epoch_len`** &mdash; *Method*.



```julia
eeg_epoch_len(eeg)
```

Return length of EEG epoch.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `epoch_len::Int64`

<a id='NeuroAnalyzer.eeg_info-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_info-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_info`** &mdash; *Method*.



```julia
eeg_info(eeg)
```

Show EEG info.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_epochs-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epochs-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epochs`** &mdash; *Method*.



```julia
eeg_epochs(eeg; epoch_n=nothing, epoch_len=nothing, average=false)
```

Split EEG into epochs. Return signal that is split either by epoch length or by number of epochs. 

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch_n::Union{Int64, Nothing}=nothing`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples
  * `average::Bool=false`: average all epochs, return one averaged epoch

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_epochs!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epochs!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epochs!`** &mdash; *Method*.



```julia
eeg_epochs!(eeg; epoch_n=nothing, epoch_len=nothing, average=false)
```

Split EEG into epochs. Return signal that is split either by epoch length or by number of epochs. 

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch_n::Union{Int64, Nothing}=nothing`: number of epochs
  * `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples
  * `average::Bool=false`: average all epochs, return one averaged epoch

<a id='NeuroAnalyzer.eeg_extract_epoch-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_extract_epoch-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_extract_epoch`** &mdash; *Method*.



```julia
eeg_extract_epoch(eeg; epoch)
```

Extract EEG epoch.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Int64`: epoch index

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_trim-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_trim-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_trim`** &mdash; *Method*.



```julia
eeg_trim(eeg; segment, remove_epochs)
```

Trim EEG signal by removing parts of the signal.

**Arguments**

  * EEG
  * `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
  * `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

**Returns**

  * EEG

<a id='NeuroAnalyzer.eeg_trim!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_trim!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_trim!`** &mdash; *Method*.



```julia
eeg_trim!(eeg; segment, remove_epochs)
```

Trim EEG signal by removing parts of the signal.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
  * `remove_epochs::Bool=true`: remove epochs containing signal to trim (remove_epochs=true) or remove signal and remove epoching

<a id='NeuroAnalyzer.eeg_edit_header-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_header-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_header`** &mdash; *Method*.



```julia
eeg_edit_header(eeg; field, value)
```

Change value of EEG header.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `field::Symbol`
  * `value::Any`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_edit_header!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_header!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_header!`** &mdash; *Method*.



```julia
eeg_edit_header!(eeg; field, value)
```

Change value of EEG header.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `field::Symbol`
  * `value::Any`

<a id='NeuroAnalyzer.eeg_show_header-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_show_header-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_show_header`** &mdash; *Method*.



```julia
eeg_show_header(eeg)
```

Show keys and values of EEG header.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_epoch-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_epoch-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_epoch`** &mdash; *Method*.



```julia
eeg_delete_epoch(eeg; epoch)
```

Remove EEG epoch(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_epoch!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_epoch!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_epoch!`** &mdash; *Method*.



```julia
eeg_delete_epoch!(eeg; epoch)
```

Remove EEG epoch(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed

<a id='NeuroAnalyzer.eeg_keep_epoch-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_epoch-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_epoch`** &mdash; *Method*.



```julia
eeg_keep_epoch(eeg; epoch)
```

Keep EEG epoch(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_keep_epoch!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_epoch!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_epoch!`** &mdash; *Method*.



```julia
eeg_keep_epoch!(eeg; epoch)
```

Keep EEG epoch(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep

<a id='NeuroAnalyzer.eeg_detect_bad-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_detect_bad-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_detect_bad`** &mdash; *Method*.



```julia
eeg_detect_bad(eeg; method, ch_t)
```

Detect bad EEG channels and epochs.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p, :var]`: detection method:

      * `:flat`: flat channel(s)
      * `:p2p`: peak-to-peak amplitude; good for detecting transient artifacts
      * `:var`: mean signal variance outside of 95%CI and variance inter-quartile outliers
      * `:rmse`: RMSE vs average channel outside of 95%CI
      * `:rmsd`: RMSD
      * `:euclid`: Euclidean distance
  * `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
  * `ftol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
  * `fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
  * `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326)
  * `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad

**Returns**

Named tuple containing:

  * `bad_m::Matrix{Bool}`: matrix of bad channels × epochs
  * `bad_epochs::Vector{Int64}`: list of bad epochs

<a id='NeuroAnalyzer.eeg_add_labels-Tuple{NeuroAnalyzer.EEG, Vector{String}}' href='#NeuroAnalyzer.eeg_add_labels-Tuple{NeuroAnalyzer.EEG, Vector{String}}'>#</a>
**`NeuroAnalyzer.eeg_add_labels`** &mdash; *Method*.



```julia
eeg_add_labels(eeg, labels)
```

Add EEG channel labels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `labels::Vector{String}`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_labels!-Tuple{NeuroAnalyzer.EEG, Vector{String}}' href='#NeuroAnalyzer.eeg_add_labels!-Tuple{NeuroAnalyzer.EEG, Vector{String}}'>#</a>
**`NeuroAnalyzer.eeg_add_labels!`** &mdash; *Method*.



```julia
eeg_add_labels!(eeg::NeuroAnalyzer.EEG, labels::Vector{String})
```

Add EEG channel labels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `labels::Vector{String}`

<a id='NeuroAnalyzer.eeg_edit_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_channel`** &mdash; *Method*.



```julia
eeg_edit_channel(eeg; channel, field, value)
```

Edit EEG channel properties.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Int64`
  * `field::Symbol`
  * `value::Any`

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_edit_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_channel!`** &mdash; *Method*.



```julia
eeg_edit_channel!(eeg; channel, field, value)
```

Edit EEG channel properties.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Int64`
  * `field::Symbol`
  * `value::Any`

<a id='NeuroAnalyzer.eeg_keep_channel_type-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_channel_type-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_channel_type`** &mdash; *Method*.



```julia
eeg_keep_channel_type(eeg; type)
```

Keep EEG channels of `type` type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:eeg`: type of channels to keep

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_keep_channel_type!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_keep_channel_type!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_keep_channel_type!`** &mdash; *Method*.



```julia
eeg_keep_channel_type!(eeg; type)
```

Keep EEG channels of `type` type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:eeg`: type of channels to keep

<a id='NeuroAnalyzer.eeg_view_note-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_view_note-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_view_note`** &mdash; *Method*.



```julia
eeg_view_note(eeg)
```

Return EEG note.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_copy-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_copy-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_copy`** &mdash; *Method*.



```julia
eeg_copy(eeg)
```

Make copy of EEG.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `eeg_copy::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_epochs_time-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epochs_time-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epochs_time`** &mdash; *Method*.



```julia
eeg_epochs_time(eeg; ts)
```

Edit EEG epochs time start.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `ts::Real`: time start in seconds

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_epochs_time!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epochs_time!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epochs_time!`** &mdash; *Method*.



```julia
eeg_epochs_time!(eeg; ts)
```

Edit EEG epochs time start.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `ts::Real`: time start in seconds

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_note-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_note-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_note`** &mdash; *Method*.



```julia
eeg_add_note(eeg; note)
```

Add EEG note.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `note::String`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_note!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_note!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_note!`** &mdash; *Method*.



```julia
eeg_add_note!(eeg; note)
```

Add EEG note.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `note::String`

<a id='NeuroAnalyzer.eeg_delete_note-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_note-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_note`** &mdash; *Method*.



```julia
eeg_delete_note(eeg)
```

Delete EEG note.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_note!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_note!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_note!`** &mdash; *Method*.



```julia
eeg_delete_note!(eeg)
```

Delete EEG note.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_replace_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_replace_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_replace_channel`** &mdash; *Method*.



```julia
eeg_replace_channel(eeg; channel, signal)
```

Replace EEG channel.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name
  * `signal::Array{Float64, 3}

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_replace_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_replace_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_replace_channel!`** &mdash; *Method*.



```julia
eeg_replace_channel!(eeg; channel, signal)
```

Replace EEG channel.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name
  * `signal::Array{Float64, 3}

<a id='NeuroAnalyzer.eeg_interpolate_channel-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_interpolate_channel-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_interpolate_channel`** &mdash; *Method*.



```julia
eeg_interpolate_channel(eeg; channel, epoch, m, q)
```

Interpolate EEG channel(s) using planar interpolation.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `interpolation_factor::Int64=100`: interpolation quality

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_interpolate_channel!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_interpolate_channel!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_interpolate_channel!`** &mdash; *Method*.



```julia
eeg_interpolate_channel!(eeg; channel, epoch, imethod, interpolation_factor)
```

Interpolate EEG channel(s) using planar interpolation.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, Abstractrange}`: channel number(s) to interpolate
  * `epoch::Union{Int64, Vector{Int64}, Abstractrange}`: epoch number(s) within to interpolate
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `interpolation_factor::Int64=100`: interpolation quality

<a id='NeuroAnalyzer.loc_flipy-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipy-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipy`** &mdash; *Method*.



```julia
loc_flipy(locs; planar, spherical)
```

Flip channel locations along y axis.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.loc_flipy!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipy!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipy!`** &mdash; *Method*.



```julia
loc_flipy!(locs; planar, spherical)
```

Flip channel locations along y axis.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.loc_flipx-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipx-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipx`** &mdash; *Method*.



```julia
loc_flipx(locs; planar, spherical)
```

Flip channel locations along x axis.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.loc_flipx!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipx!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipx!`** &mdash; *Method*.



```julia
loc_flipx!(locs; planar, spherical)
```

Flip channel locations along x axis.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.loc_flipz-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipz-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipz`** &mdash; *Method*.



```julia
loc_flipz(locs)
```

Flip channel locations along z axis.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.loc_flipz!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_flipz!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_flipz!`** &mdash; *Method*.



```julia
loc_flipz!(eeg)
```

Flip channel locations along z axis.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_channel_type-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_channel_type-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_channel_type`** &mdash; *Method*.



```julia
eeg_channel_type(eeg; channel, type)
```

Change EEG channel type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`
  * `type::String`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_channel_type!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_channel_type!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_channel_type!`** &mdash; *Method*.



```julia
eeg_channel_type!(eeg; channel, new_name)
```

Change EEG channel type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`
  * `type::String`

<a id='NeuroAnalyzer.eeg_edit_electrode-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_electrode-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_electrode`** &mdash; *Method*.



```julia
eeg_edit_electrode(eeg; <keyword arguments>)
```

Edit EEG electrode.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{String, Int64}`: channel number or name
  * `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
  * `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
  * `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
  * `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
  * `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
  * `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
  * `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
  * `name::String=""`: channel name
  * `type::String=""`: channel type

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_edit_electrode!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_electrode!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_electrode!`** &mdash; *Method*.



```julia
eeg_edit_electrode!(eeg; <keyword arguments>)
```

Edit EEG electrode.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{String, Int64}`: channel number or name
  * `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
  * `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
  * `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
  * `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
  * `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
  * `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
  * `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
  * `name::String=""`: channel name
  * `type::String=""`: channel type

<a id='NeuroAnalyzer.eeg_electrode_loc-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_electrode_loc-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_electrode_loc`** &mdash; *Method*.



```julia
eeg_electrode_loc(eeg; channel, output)
```

Return locations of EEG channel electrode.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, String}`: channel number or name
  * `output::Bool=true`: print output if true

**Returns**

Named tuple containing:

  * `theta::Float64`: polar planar theta coordinate
  * `radius::Float64`: polar planar radius coordinate
  * `x::Float64`: Cartesian X spherical coordinate
  * `y::Float64`: Cartesian Y spherical coordinate
  * `z::Float64`: Cartesian Z spherical coordinate
  * `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Float64`: spherical radius, the distance from the origin to the point
  * `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

<a id='NeuroAnalyzer.loc_swapxy-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_swapxy-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_swapxy`** &mdash; *Method*.



```julia
loc_swapxy(locs; planar, spherical)
```

Swap channel locations x and y axes.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.loc_swapxy!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_swapxy!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_swapxy!`** &mdash; *Method*.



```julia
loc_swapxy!(locs; planar, spherical)
```

Swap channel locations x and y axes.

**Arguments**

  * `locs::DataFrame`
  * `planar::Bool=true`: modify planar coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.loc_sph2cart-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_sph2cart-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_sph2cart`** &mdash; *Method*.



```julia
loc_sph2cart(locs)
```

Convert spherical locations to Cartesian.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.loc_sph2cart!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_sph2cart!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_sph2cart!`** &mdash; *Method*.



```julia
loc_sph2cart!(locs)
```

Convert spherical locations to Cartesian.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.loc_cart2sph-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_cart2sph-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_cart2sph`** &mdash; *Method*.



```julia
loc_cart2sph(locs)
```

Convert Cartesian locations to spherical.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.loc_cart2sph!-Tuple{DataFrame}' href='#NeuroAnalyzer.loc_cart2sph!-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.loc_cart2sph!`** &mdash; *Method*.



```julia
loc_cart2sph!(locs)
```

Convert Cartesian locations to spherical.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.eeg_view_markers-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_view_markers-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_view_markers`** &mdash; *Method*.



```julia
eeg_view_markers(eeg)
```

Show markers.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_marker-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_marker-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_marker`** &mdash; *Method*.



```julia
eeg_delete_marker(eeg; n)
```

Delete marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `n::Int64`: marker number

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_delete_marker!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_delete_marker!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_delete_marker!`** &mdash; *Method*.



```julia
eeg_delete_marker!(eeg; n)
```

Delete marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `n::Int64`: marker number

<a id='NeuroAnalyzer.eeg_add_marker-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_marker-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_marker`** &mdash; *Method*.



```julia
eeg_add_marker(eeg; id, start, len, desc)
```

Add marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `id::String`: marker ID
  * `start::Int64`: marker time in samples
  * `len::Int64`: marker length in samples
  * `desc::String`: marker description
  * `channel::Int64`: channel number, if 0 then marker is related to all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_marker!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_marker!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_marker!`** &mdash; *Method*.



```julia
eeg_add_marker!(eeg; id, start, len, desc)
```

Add marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `id::String`: marker ID
  * `start::Int64`: marker time in samples
  * `len::Int64`: marker length in samples
  * `desc::String`: marker description
  * `channel::Int64`: channel number, if 0 then marker is related to all channels

<a id='NeuroAnalyzer.eeg_channel_idx-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_channel_idx-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_channel_idx`** &mdash; *Method*.



```julia
eeg_channel_idx(eeg; type=:eeg)
```

Return index of EEG channels of `type` type.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Vector{Symbol}=:all`: channel type

**Returns**

  * `channel_n::Int64`

<a id='NeuroAnalyzer.eeg_vch-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_vch-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_vch`** &mdash; *Method*.



```julia
eeg_vch(eeg; f)
```

Calculate virtual channel using formula `f`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the EEG object

**Returns**

  * `vc::Array{Float64, 3}`: single channel × time × epochs

<a id='NeuroAnalyzer.eeg_edit_marker-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_marker-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_marker`** &mdash; *Method*.



```julia
eeg_edit_marker(eeg; n, id, start, len, desc)
```

Edit EEG marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `n::Int64`: marker number
  * `id::String`: marker ID
  * `start::Int64`: marker time in samples
  * `len::Int64`: marker length in samples
  * `desc::String`: marker description
  * `channel::Int64`: channel number, if 0 then marker is related to all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_edit_marker!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_edit_marker!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_edit_marker!`** &mdash; *Method*.



```julia
eeg_edit_marker!(eeg; n, id, start, len, desc)
```

Edit EEG marker.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `n::Int64`: marker number
  * `id::String`: marker ID
  * `start::Int64`: marker time in samples
  * `len::Int64`: marker length in samples
  * `desc::String`: marker description
  * `channel::Int64`: channel number, if 0 then marker is related to all channels


<a id='EEG-process'></a>

<a id='EEG-process-1'></a>

## EEG process

<a id='NeuroAnalyzer.eeg_reference_ch-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_ch-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_ch`** &mdash; *Method*.



```julia
eeg_reference_ch(eeg; channel, med)
```

Reference to selected channel(s). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reference_ch!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_ch!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_ch!`** &mdash; *Method*.



```julia
eeg_reference_ch!(eeg; channel, med)
```

Reference to selected channel(s). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.eeg_reference_car-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_car-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_car`** &mdash; *Method*.



```julia
eeg_reference_car(eeg; exclude_fpo, exclude_current, med)
```

Reference to common average reference. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR calculation
  * `exclude_current::Bool=true`: exclude current channel from CAR calculation
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reference_car!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_car!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_car!`** &mdash; *Method*.



```julia
eeg_reference_car!(eeg; exclude_fpo, exclude_current, med)
```

Reference to common average reference. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `exclude_fpo::Bool=false`: exclude Fp1, Fp2, O1, O2 from CAR mean calculation
  * `exclude_current::Bool=true`: exclude current channel from CAR mean calculation
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.eeg_derivative-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_derivative-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_derivative`** &mdash; *Method*.



```julia
eeg_derivative(eeg; channel)
```

Return the derivative of EEG channel(s) with length same as the signal.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_derivative!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_derivative!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_derivative!`** &mdash; *Method*.



```julia
eeg_derivative!(eeg; channel)
```

Return the derivative of EEG channel(s) with length same as the signal.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_detrend-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_detrend-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_detrend`** &mdash; *Method*.



```julia
eeg_detrend(eeg; channel, type, offset, order, f)
```

Perform piecewise detrending of EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `type::Symbol=:linear`

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` is subtracted
      * `:loess`: fit and subtract LOESS approximation
      * `:hp`: use HP filter
  * `offset::Real=0`: constant for `:constant` detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_detrend!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_detrend!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_detrend!`** &mdash; *Method*.



```julia
eeg_detrend!(eeg; channel, type, offset, order, span)
```

Perform piecewise detrending of EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `type::Symbol=:linear`

      * `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
      * `:linear`: linear trend is subtracted from `signal`
      * `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
      * `:poly`: polynomial of `order` order is subtracted
      * `:loess`: fit and subtract LOESS approximation
      * `:hp`: use HP filter
  * `offset::Real=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

<a id='NeuroAnalyzer.eeg_taper-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_taper-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_taper`** &mdash; *Method*.



```julia
eeg_taper(eeg; channel, taper)
```

Taper EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `taper::Union{Vector{Real, Vector{ComplexF64}}`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_taper!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_taper!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_taper!`** &mdash; *Method*.



```julia
eeg_taper!(eeg; channel, taper)
```

Taper EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `taper::Union{Vector{<:Real}, Vector{ComplexF64}}`

<a id='NeuroAnalyzer.eeg_demean-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_demean-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_demean`** &mdash; *Method*.



```julia
eeg_demean(eeg; channel)
```

Remove mean value (DC offset).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_demean!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_demean!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_demean!`** &mdash; *Method*.



```julia
eeg_demean!(eeg; channel)
```

Remove mean value (DC offset).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_normalize-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_normalize-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_normalize`** &mdash; *Method*.



```julia
eeg_normalize(eeg; channel, method)
```

Normalize EEG channel(s)

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `method::Symbol`: method for normalization, see `s_normalize()` for details

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_normalize!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_normalize!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_normalize!`** &mdash; *Method*.



```julia
eeg_normalize!(eeg; channel, method)
```

Normalize EEG channel(s)

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `method::Symbol`: method for normalization, see `s_normalize()` for details

<a id='NeuroAnalyzer.eeg_add_noise-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_noise-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_noise`** &mdash; *Method*.



```julia
eeg_add_noise(eeg; channel)
```

Add random noise to EEG channel(s)

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_add_noise!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_add_noise!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_add_noise!`** &mdash; *Method*.



```julia
eeg_add_noise!(eeg; channel)
```

Add random noise to EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_filter-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_filter-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_filter`** &mdash; *Method*.



```julia
eeg_filter(eeg; <keyword arguments>)
```

Apply filtering to EEG channel(s).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order`
      * `:conv`: convolution
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
  * `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, `:twopass`), for causal filter use `:onepass`
  * `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
  * `window::Union{AbstractVector, Nothing} - window, required for FIR filter, weighting window for`:mavg`and`:mmed`

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_filter!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_filter!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_filter!`** &mdash; *Method*.



```julia
eeg_filter!(eeg; <keyword arguments>)
```

Apply filtering to EEG channel(s).

**Arguments**

  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
      * `:mavg`: moving average (with threshold)
      * `:mmed`: moving median (with threshold)
      * `:poly`: polynomial of `order`
  * `ftype::Symbol`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `order::Int64=8`: filter order, number of taps for `:remez` filter, k-value for `:mavg` and `:mmed` (window length = 2 × k + 1)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez filters
  * `dir:Symbol=:twopass`: filter direction (`:onepass`, :onepass_reverse, `:twopass`), for causal filter use `:onepass`
  * `t::Real`: threshold for `:mavg` and `:mmed` filters; threshold = threshold * std(signal) + mean(signal) for `:mavg` or threshold = threshold * std(signal) + median(signal) for `:mmed` filter
  * `window::Union{Vector{<:Real}, Nothing} - window, required for`:fir` filter
  * `preview::Bool=false`: plot filter response

<a id='NeuroAnalyzer.eeg_pca-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pca-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pca`** &mdash; *Method*.



```julia
eeg_pca(eeg; channel, n)
```

Perform principal component analysis (PCA).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `n::Int64`: number of PCs to calculate

**Returns**

Named tuple containing:

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
  * `pc_m::Vector{Float64}`: PC mean
  * `pca::MultivariateStats.PCA{Float64}`: PC model

<a id='NeuroAnalyzer.eeg_pca_reconstruct-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pca_reconstruct-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pca_reconstruct`** &mdash; *Method*.



```julia
eeg_pca_reconstruct(eeg; channel)
```

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_pca_reconstruct!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pca_reconstruct!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pca_reconstruct!`** &mdash; *Method*.



```julia
eeg_pca_reconstruct!(eeg; channel)
```

Reconstruct EEG signals using embedded PCA components (`:pc`) and model (`:pca`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_pca_reconstruct-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, PCA{Float64}}' href='#NeuroAnalyzer.eeg_pca_reconstruct-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, PCA{Float64}}'>#</a>
**`NeuroAnalyzer.eeg_pca_reconstruct`** &mdash; *Method*.



```julia
eeg_pca_reconstruct(eeg, pc, pca; channel)
```

Reconstruct EEG signals using external PCA components (`pc` and `pca`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pca::MultivariateStats.PCA{Float64}`: PC model
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_pca_reconstruct!-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, PCA{Float64}}' href='#NeuroAnalyzer.eeg_pca_reconstruct!-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, PCA{Float64}}'>#</a>
**`NeuroAnalyzer.eeg_pca_reconstruct!`** &mdash; *Method*.



```julia
eeg_pca_reconstruct!(eeg, pc, pca; channel)
```

Reconstruct EEG signals using external PCA components (`pc` and `pca`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pca::MultivariateStats.PCA{Float64}`: PC model
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_ica-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ica-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ica`** &mdash; *Method*.



```julia
eeg_ica(eeg; <keyword arguments>)
```

Perform independent component analysis (ICA).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `n::Int64`: number of ICs
  * `tol::Float64=1.0e-6`: tolerance for ICA
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

**Returns**

Named tuple containing:

  * `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
  * `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)

<a id='NeuroAnalyzer.eeg_ica_reconstruct-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ica_reconstruct-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ica_reconstruct`** &mdash; *Method*.



```julia
eeg_ica_reconstruct(eeg; channel, ic)
```

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `ic::Union{Int64, Vector{Int64}, AbstractRange}`: list of ICs to remove

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_ica_reconstruct!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ica_reconstruct!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ica_reconstruct!`** &mdash; *Method*.



```julia
eeg_ica_reconstruct!(eeg; channel, ic)
```

Reconstruct EEG signals using embedded ICA components (`:ica` and `:ica_mw`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

<a id='NeuroAnalyzer.eeg_ica_reconstruct-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroAnalyzer.eeg_ica_reconstruct-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroAnalyzer.eeg_ica_reconstruct`** &mdash; *Method*.



```julia
eeg_ica_reconstruct(eeg, ica, ica_mw; channel, ic)
```

Reconstruct EEG signals using external ICA components (`ica` and `ica_mw`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
  * `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_ica_reconstruct!-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, Array{Float64, 3}}' href='#NeuroAnalyzer.eeg_ica_reconstruct!-Tuple{NeuroAnalyzer.EEG, Array{Float64, 3}, Array{Float64, 3}}'>#</a>
**`NeuroAnalyzer.eeg_ica_reconstruct!`** &mdash; *Method*.



```julia
eeg_ica_reconstruct!(eeg, ica, ica_mw; channel, ic)
```

Reconstruct EEG signals using external ICA components (`ica` and `ica_mw`).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `ica::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
  * `ica_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `ic::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

<a id='NeuroAnalyzer.eeg_average' href='#NeuroAnalyzer.eeg_average'>#</a>
**`NeuroAnalyzer.eeg_average`** &mdash; *Function*.



```julia
eeg_average(eeg; channel)
```

Return the average signal of EEG channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_average!' href='#NeuroAnalyzer.eeg_average!'>#</a>
**`NeuroAnalyzer.eeg_average!`** &mdash; *Function*.



```julia
eeg_average!(eeg; channel)
```

Return the average signal of EEG channels.  

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_average-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_average-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_average`** &mdash; *Method*.



```julia
eeg_average(eeg1, eeg2)
```

Return the average signal of all `eeg1` and `eeg2` channels.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_invert_polarity-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_invert_polarity-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_invert_polarity`** &mdash; *Method*.



```julia
eeg_invert_polarity(eeg; channel)
```

Invert polarity.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_invert_polarity!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_invert_polarity!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_invert_polarity!`** &mdash; *Method*.



```julia
eeg_invert_polarity!(eeg; channel)
```

Invert polarity.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: channel(s) to invert

<a id='NeuroAnalyzer.eeg_resample-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_resample-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_resample`** &mdash; *Method*.



```julia
eeg_resample(eeg; new_sr)
```

Resample (up- or down-sample).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_resample!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_resample!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_resample!`** &mdash; *Method*.



```julia
eeg_resample!(eeg; new_sr)
```

Resample (up- or down-sample).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroAnalyzer.eeg_upsample-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_upsample-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_upsample`** &mdash; *Method*.



```julia
eeg_upsample(eeg; new_sr)
```

Upsample EEG.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_upsample!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_upsample!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_upsample!`** &mdash; *Method*.



```julia
eeg_upsample!(eeg; new_sr)
```

Upsample EEG.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_downsample-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_downsample-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_downsample`** &mdash; *Method*.



```julia
eeg_downsample(eeg; new_sr)
```

Downsample EEG.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_downsample!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_downsample!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_downsample!`** &mdash; *Method*.



```julia
eeg_downsample!(eeg; new_sr)
```

Downsample EEG.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroAnalyzer.eeg_wdenoise-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.eeg_wdenoise-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.eeg_wdenoise`** &mdash; *Method*.



```julia
eeg_wdenoise(eeg; channel, wt)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_wdenoise!-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.eeg_wdenoise!-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.eeg_wdenoise!`** &mdash; *Method*.



```julia
eeg_wdenoise!(eeg; channel, wt)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

<a id='NeuroAnalyzer.eeg_reference_a-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_a-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_a`** &mdash; *Method*.



```julia
eeg_reference_a(eeg; type, med)
```

Reference to auricular (A1, A2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:l`:

      * `:l`: linked - average of A1 and A2
      * `:i`: ipsilateral - A1 for left channels, A2 for right channels
      * `:c`: contraletral - A1 for right channels, A2 for left channels
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reference_a!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_a!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_a!`** &mdash; *Method*.



```julia
eeg_reference_a!(eeg; type, med)
```

Reference to auricular (A1, A2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:l`:

      * `:l`: linked - average of A1 and A2
      * `:i`: ipsilateral - A1 for left channels, A2 for right channels
      * `:c`: contraletral - A1 for right channels, A2 for left channels
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.eeg_reference_m-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_m-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_m`** &mdash; *Method*.



```julia
eeg_reference_m(eeg; type, med)
```

Reference to mastoid (M1, M2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:l`:

      * `:l`: linked - average of M1 and M2
      * `:i`: ipsilateral - M1 for left channels, M2 for right channels
      * `:c`: contraletral - M1 for right channels, M2 for left channels
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reference_m!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_m!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_m!`** &mdash; *Method*.



```julia
eeg_reference_m!(eeg; type, med)
```

Reference to mastoid (M1, M2) channels. Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `type::Symbol=:l`:

      * `:l`: linked - average of M1 and M2
      * `:i`: ipsilateral - M1 for left channels, M2 for right channels
      * `:c`: contraletral - M1 for right channels, M2 for left channels
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.eeg_fftdenoise-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_fftdenoise-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_fftdenoise`** &mdash; *Method*.



```julia
eeg_fftdenoise(eeg; channel, pad, threshold)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `threshold::Int64=100`: PSD threshold for keeping frequency components

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_fftdenoise!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_fftdenoise!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_fftdenoise!`** &mdash; *Method*.



```julia
eeg_fftdenoise!(eeg; channel, pad, threshold)
```

Perform wavelet denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `threshold::Int64=100`: PSD threshold for keeping frequency components

<a id='NeuroAnalyzer.eeg_reference_plap-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_plap-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_plap`** &mdash; *Method*.



```julia
eeg_reference_plap(eeg; nn, weights)
```

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `nn::Int64=4`: use `nn` adjacent electrodes
  * `weights::Bool=false`: use mean of `nn` nearest channels if false; if true, mean of `nn` nearest channels is weighted by distance to the referenced channel
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_reference_plap!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_reference_plap!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_reference_plap!`** &mdash; *Method*.



```julia
eeg_reference_plap!(eeg; nn, weights)
```

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal (EEG/MEG, depending on `eeg.eeg_header[:signal_type]`) channels are processed.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `nn::Int64=4`: use `nn` adjacent electrodes
  * `weights::Bool=false`: use distance weights; use mean of nearest channels if false
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.eeg_zero-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_zero-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_zero`** &mdash; *Method*.



```julia
eeg_zero(eeg)
```

Zero EEG channels at the beginning and at the end.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_zero!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_zero!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_zero!`** &mdash; *Method*.



```julia
eeg_zero!(eeg)
```

Zero EEG channels at the beginning and at the end.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_wbp-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_wbp-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_wbp`** &mdash; *Method*.



```julia
eeg_wbp(eeg; channel, pad, frq, ncyc, demean)
```

Perform wavelet bandpass filtering.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_wbp!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_wbp!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_wbp!`** &mdash; *Method*.



```julia
eeg_wbp!(eeg; channel, pad, frq, ncyc, demean)
```

Perform wavelet bandpass filtering.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet
  * `demean::Bool=true`: demean signal prior to analysis

<a id='NeuroAnalyzer.eeg_cbp-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cbp-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cbp`** &mdash; *Method*.



```julia
eeg_cbp(eeg; channel, pad, frq, demean)
```

Perform convolution bandpass filtering.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `demean::Bool=true`: demean signal prior to analysis

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_cbp!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cbp!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cbp!`** &mdash; *Method*.



```julia
eeg_cbp!(eeg; channel, pad, frq, demean)
```

Perform convolution bandpass filtering.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Tuple{Real, Real}`: filter frequency
  * `demean::Bool=true`: demean signal prior to analysis

<a id='NeuroAnalyzer.eeg_denoise_wien-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_denoise_wien-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_denoise_wien`** &mdash; *Method*.



```julia
eeg_denoise_wien(eeg; channel)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_denoise_wien!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_denoise_wien!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_denoise_wien!`** &mdash; *Method*.



```julia
eeg_denoise_wien!(eeg; channel)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels

<a id='NeuroAnalyzer.eeg_scale-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_scale-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_scale`** &mdash; *Method*.



```julia
eeg_scale(eeg; channel, factor)
```

Multiply EEG channel(s) by `factor`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `factor::Real`: channel signal is multiplied by `factor`

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`

<a id='NeuroAnalyzer.eeg_scale!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_scale!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_scale!`** &mdash; *Method*.



```julia
eeg_scale!(eeg; channel, factor)
```

Multiply EEG channel(s) by `factor`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
  * `factor::Real`: channel signal is multiplied by `factor`


<a id='EEG-analyze'></a>

<a id='EEG-analyze-1'></a>

## EEG analyze

<a id='NeuroAnalyzer.eeg_total_power-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_total_power-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_total_power`** &mdash; *Method*.



```julia
eeg_total_power(eeg, channel, mt)
```

Calculate total power.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

  * `stp::Matrix{Float64}`: total power for each channel per epoch

<a id='NeuroAnalyzer.eeg_band_power-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_band_power-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_band_power`** &mdash; *Method*.



```julia
eeg_band_power(eeg; channel, f, mt)
```

Calculate absolute band power between two frequencies

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

  * `sbp::Matrix{Float64}`: band power for each channel per epoch

<a id='NeuroAnalyzer.eeg_cov-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cov-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cov`** &mdash; *Method*.



```julia
eeg_cov(eeg; channel, norm)
```

Calculate covariance matrix.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `norm::Bool=true`: normalize matrix

**Returns**

  * `cov_mat::Array{Float64, 3}`: covariance matrix for each epoch

<a id='NeuroAnalyzer.eeg_cor-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cor-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cor`** &mdash; *Method*.



```julia
eeg_cor(eeg; channel, norm)
```

Calculate correlation coefficients.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `norm::Bool=true`: normalize matrix

**Returns**

  * `cov_mat::Array{Float64, 3}`: correlation matrix for each epoch

<a id='NeuroAnalyzer.eeg_xcov-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_xcov-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_xcov`** &mdash; *Method*.



```julia
eeg_xcov(eeg; channel, lag, demean, norm)
```

Calculate cross-covariance.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize cross-covariance

**Returns**

Named tuple containing:

  * `xcov::Matrix{Float64}`: ch1-ch1, ch1-ch2, ch1-ch3, etc.
  * `lags::Vector{Float64}`: lags in ms

<a id='NeuroAnalyzer.eeg_xcov-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_xcov-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_xcov`** &mdash; *Method*.



```julia
eeg_xcov(eeg1, eeg2; channel1, channel2, epoch1, epoch2, lag, demean, norm)
```

Calculate cross-covariance between two EEG signals.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
  * `lag::Int64=1`: lags range is `-lag:lag`
  * `demean::Bool=false`: demean signal prior to analysis
  * `norm::Bool=false`: normalize cross-covariance

**Returns**

Named tuple containing:

  * `xcov::Array{Float64, 3}`
  * `lags::Vector{Float64}`: lags in ms

<a id='NeuroAnalyzer.eeg_psd-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_psd-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_psd`** &mdash; *Method*.



```julia
eeg_psd(eeg; channel, norm, mt)
```

Calculate power spectrum density.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `psd_pow::Array{Float64, 3}`:powers
  * `psd_frq::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.eeg_stationarity-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_stationarity-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_stationarity`** &mdash; *Method*.



```julia
eeg_stationarity(eeg; channel, window, method)
```

Calculate stationarity.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `window::Int64=10`: time window in samples
  * `method::Symbol=:euclid`: stationarity method:

      * `:mean`
      * `:var`: variance stationarity
      * `:euclid`
      * `:hilbert`
      * `:adf`

**Returns**

  * `stationarity::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_mi-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_mi-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_mi`** &mdash; *Method*.



```julia
eeg_mi(eeg; channel)
```

Calculate mutual information between EEG channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_mi-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_mi-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_mi`** &mdash; *Method*.



```julia
eeg_mi(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate mutual information between two EEG channels.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

  * `mi::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_entropy-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_entropy-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_entropy`** &mdash; *Method*.



```julia
eeg_entropy(eeg; channel)
```

Calculate entropy.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

Named tuple containing:

  * `ent::Array{Float64, 2}`
  * `s_ent::Array{Float64, 2}`: Shanon entropy
  * `le_ent::Array{Float64, 2}`: log energy entropy

<a id='NeuroAnalyzer.eeg_negentropy-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_negentropy-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_negentropy`** &mdash; *Method*.



```julia
eeg_negentropy(eeg; channel)
```

Calculate negentropy.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `ne::Matrix{Float64}`

<a id='NeuroAnalyzer.eeg_band-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_band-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_band`** &mdash; *Method*.



```julia
eeg_band(eeg, band)
```

Return frequency limits for a `band`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `band::Symbol`: band range name: `:list`, `:total`, `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower`, `:gamma_higher`.

**Returns**

  * `band_frequency::Tuple{Real, Real}`

<a id='NeuroAnalyzer.eeg_band-Tuple{Int64}' href='#NeuroAnalyzer.eeg_band-Tuple{Int64}'>#</a>
**`NeuroAnalyzer.eeg_band`** &mdash; *Method*.



```julia
eeg_band(fs, band)
```

Return frequency limits for a `band`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `band::Symbol`: band range name: `:list`, `:total`, `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower`, `:gamma_higher`.

**Returns**

  * `band_frequency::Tuple{Real, Real}`

<a id='NeuroAnalyzer.eeg_tcoherence-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tcoherence-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tcoherence`** &mdash; *Method*.



```julia
eeg_tcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate coherence (mean over time) and MSC (magnitude-squared coherence).

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
  * `pad::Int64=0`: number of zeros to add signal for FFT

**Returns**

Named tuple containing:

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `ic::Array{Float64, 3}`: imaginary part of coherence

<a id='NeuroAnalyzer.eeg_freqs-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_freqs-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_freqs`** &mdash; *Method*.



```julia
eeg_freqs(eeg)
```

Return vector of frequencies and Nyquist frequency.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

Named tuple containing:

  * `hz::Vector{Float64}`
  * `nyquist::Float64`

<a id='NeuroAnalyzer.eeg_difference-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_difference-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_difference`** &mdash; *Method*.



```julia
eeg_difference(eeg; channel, n, method)
```

Calculate mean difference and its 95% CI between EEG channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:absdiff`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `signals_statistic::Matrix{Float64}`
  * `signals_statistic_single::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroAnalyzer.eeg_pick-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pick-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pick`** &mdash; *Method*.



```julia
eeg_picks(eeg; pick)
```

Return set of channel indices corresponding with `pick` of electrodes

**Arguments**

  * `pick::Vector{Symbol}`: pick of electrodes, e.g. `:list`, `:central` (or `:c`), `:left` (or `:l`), `:right` (or `:r`), `:frontal` (or `:f`), `:temporal` (or `:t`), `:parietal` (or `:p`), `:occipital` (or `:o`); picks may be combined, e.g. `[:left, :frontal]`

**Returns**

  * `channels::Vector{Int64}`: channel numbers

<a id='NeuroAnalyzer.eeg_epochs_stats-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_epochs_stats-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_epochs_stats`** &mdash; *Method*.



```julia
eeg_epochs_stats(eeg)
```

Calculate epochs statistics.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

Named tuple containing:

  * `e_mean::Vector(Float64)`: mean
  * `e_median::Vector(Float64)`: median
  * `e_std::Vector(Float64)`: standard deviation
  * `e_var::Vector(Float64)`: variance
  * `e_kurt::Vector(Float64)`: kurtosis
  * `e_skew::Vector(Float64)`: skewness
  * `e_mean_diff::Vector(Float64)`: mean diff value
  * `e_median_diff::Vector(Float64)`: median diff value
  * `e_max_dif::Vector(Float64)`: max difference
  * `e_dev_mean::Vector(Float64)`: deviation from channel mean

<a id='NeuroAnalyzer.eeg_spectrogram-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:CWT' href='#NeuroAnalyzer.eeg_spectrogram-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:CWT'>#</a>
**`NeuroAnalyzer.eeg_spectrogram`** &mdash; *Method*.



```julia
eeg_spectrogram(eeg; channel, norm, mt, st, demean)
```

Calculate spectrogram.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `method::Symbol=:standard`: method of calculating spectrogram:

      * `:standard`: standard
      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `pad::Int64=0`: number of zeros to add
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limits
  * `frq_n::Int64=0`: number of frequencies, default is length(frq*lim[1]:frq*lim[2])
  * `norm::Bool=true`: normalize powers to dB
  * `demean::Bool=true`: demean signal prior to analysis
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq*n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq*n) for frq === :lin
  * `wt<:CWT=wavelet(Morlet(π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

Named tuple containing:

  * `s_pow::Array{Float64, 3}`
  * `s_frq::Vector{Float64}`
  * `s_t::Vector{Float64}`

<a id='NeuroAnalyzer.eeg_spectrum-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_spectrum-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_spectrum`** &mdash; *Method*.



```julia
eeg_spectrum(eeg; channel, pad, h)
```

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `h::Bool=false`: use Hilbert transform for calculations instead of FFT
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
  * `amp::Array{Float64, 3}`: amplitudes
  * `pow::Array{Float64, 3}`: powers
  * `pha::Array{Float64, 3}: phase angles

<a id='NeuroAnalyzer.eeg_s2t-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_s2t-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_s2t`** &mdash; *Method*.



```julia
eeg_s2t(eeg; t)
```

Convert time in samples to seconds.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `t::Int64`: time in samples

**Returns**

  * `t_s::Float64`: time in seconds

<a id='NeuroAnalyzer.eeg_t2s-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_t2s-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_t2s`** &mdash; *Method*.



```julia
eeg_t2s(eeg; t)
```

Convert time in seconds to samples.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `t::Real`: time in seconds

**Returns**

  * `t_s::Int64`: time in samples

<a id='NeuroAnalyzer.eeg_channels_stats-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_channels_stats-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_channels_stats`** &mdash; *Method*.



```julia
eeg_channels_stats(eeg)
```

Calculate channels statistics per epoch.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`

**Returns**

Named tuple containing:

  * `c_mean::Matrix(Float64)`: mean
  * `c_median::Matrix(Float64)`: median
  * `c_std::Matrix(Float64)`: standard deviation
  * `c_var::Matrix(Float64)`: variance
  * `c_kurt::Matrix(Float64)`: kurtosis
  * `c_skew::Matrix(Float64)`: skewness
  * `c_mean_diff::Matrix(Float64)`: mean diff value
  * `c_median_diff::Matrix(Float64)`: median diff value
  * `c_max_dif::Matrix(Float64)`: max difference
  * `c_dev_mean::Matrix(Float64)`: deviation from channel mean

<a id='NeuroAnalyzer.eeg_snr-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_snr-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_snr`** &mdash; *Method*.



```julia
eeg_snr(eeg; channel)
```

Calculate SNR of `eeg` channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `snr::Matrix(Float64)`: SNR for each channel per epoch

<a id='NeuroAnalyzer.eeg_standardize-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_standardize-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_standardize`** &mdash; *Method*.



```julia
eeg_standardize(eeg; channel)
```

Standardize channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `eeg_new::NeuroAnalyzer.EEG`: standardized EEG
  * `scaler::Matrix{Float64}`: standardizing matrix

<a id='NeuroAnalyzer.eeg_standardize!-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_standardize!-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_standardize!`** &mdash; *Method*.



```julia
eeg_standardize!(eeg; channel)
```

Standardize channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `scaler::Matrix{Float64}`: standardizing matrix

<a id='NeuroAnalyzer.eeg_fconv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_fconv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_fconv`** &mdash; *Method*.



```julia
eeg_fconv(eeg; channel, kernel, norm)
```

Perform convolution in the frequency domain.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
  * `norm::Bool=false`: normalize kernel
  * `pad::Int64=0`: number of zeros to add signal for FFT

**Returns**

  * `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal

<a id='NeuroAnalyzer.eeg_tconv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tconv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tconv`** &mdash; *Method*.



```julia
eeg_tconv(eeg; channel, kernel)
```

Perform convolution in the time domain.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

**Returns**

  * `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal

<a id='NeuroAnalyzer.eeg_dft-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_dft-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_dft`** &mdash; *Method*.



```julia
eeg_dft(eeg; channel)
```

Return FFT and DFT sample frequencies for a DFT.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `pad::Int64=0`: number of zeros to add signal for FFT

**Returns**

Named tuple containing:

  * `sfft::Array{ComplexF64, 3}`: FFT
  * `sf::Vector{Float64}`: sample frequencies

<a id='NeuroAnalyzer.eeg_msci95-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_msci95-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_msci95`** &mdash; *Method*.



```julia
eeg_msci95(eeg; channel, n, method)
```

Calculate mean, standard deviation and 95% confidence interval for EEG channels.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

**Returns**

Named tuple containing:

  * `s_m::Matrix{Float64}`: mean
  * `s_s::Matrix{Float64}`: standard deviation
  * `s_u::Matrix{Float64}`: upper 95% CI
  * `s_l::Matrix{Float64}`: lower 95% CI

<a id='NeuroAnalyzer.eeg_mean-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_mean-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_mean`** &mdash; *Method*.



```julia
eeg_mean(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculates mean, standard deviation and 95% confidence interval for two EEG channels.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2:NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

Named tuple containing:

  * `s_m::Matrix{Float64}`: mean by epochs
  * `s_s::Matrix{Float64}`: standard deviation by epochs
  * `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
  * `s_l::Matrix{Float64}`: lower 95% CI bound by epochs

<a id='NeuroAnalyzer.eeg_difference-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_difference-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_difference`** &mdash; *Method*.



```julia
eeg_difference(eeg1, eeg2; channel1, channel2, epoch1, epoch2, n, method)
```

Calculates mean difference and 95% confidence interval for two EEG channels.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2:NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:absdiff, :diff2int]`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `s_stat::Matrix{Float64}`
  * `s_stat_single::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroAnalyzer.eeg_acov-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_acov-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_acov`** &mdash; *Method*.



eeg_acov(eeg; channel, lag=1, demean=false, norm=false)

Calculate autocovariance.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `lag::Int64`: lags range is `-lag:lag`
  * `demean::Bool`: demean eeg prior to analysis
  * `norm::Bool`: normalize autocovariance

**Returns**

Named tuple containing:

  * `acov::Matrix{Float64}`
  * `lags::Vector{Float64}`: lags in ms

<a id='NeuroAnalyzer.eeg_tenv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tenv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tenv`** &mdash; *Method*.



```julia
eeg_tenv(eeg; channel, d)
```

Calculate temporal envelope.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env::Array{Float64, 3}`: temporal envelope
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_tenv_mean-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tenv_mean-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tenv_mean`** &mdash; *Method*.



```julia
eeg_tenv_mean(eeg; channel, dims, d)
```

Calculate temporal envelope: mean and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_tenv_median-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tenv_median-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tenv_median`** &mdash; *Method*.



```julia
eeg_tenv_median(eeg; channel, dims, d)
```

Calculate temporal envelope of `eeg`: median and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_penv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_penv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_penv`** &mdash; *Method*.



```julia
eeg_penv(eeg; channel, d)
```

Calculate power (in dB) envelope.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env::Array{Float64, 3}`: power spectrum envelope
  * `p_env_frq::Vector{Float64}`: frequencies for each envelope

<a id='NeuroAnalyzer.eeg_penv_mean-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_penv_mean-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_penv_mean`** &mdash; *Method*.



```julia
eeg_penv_mean(eeg; channel, dims, d)
```

Calculate power (in dB) envelope: mean and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroAnalyzer.eeg_penv_median-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_penv_median-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_penv_median`** &mdash; *Method*.



```julia
eeg_penv_median(eeg; channel, dims, d)
```

Calculate power (in dB) envelope: median and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroAnalyzer.eeg_senv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_senv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_senv`** &mdash; *Method*.



```julia
eeg_senv(eeg; channel, d, mt, t)
```

Calculate spectral envelope.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

**Returns**

Named tuple containing:

  * `s_env::Array{Float64, 3}`: spectral envelope
  * `s_env_t::Vector{Float64}`: spectrogram time

<a id='NeuroAnalyzer.eeg_senv_mean-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_senv_mean-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_senv_mean`** &mdash; *Method*.



```julia
eeg_senv_mean(eeg; channel, dims, d, mt, t)
```

Calculate spectral envelope: mean and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: mean
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroAnalyzer.eeg_senv_median-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_senv_median-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_senv_median`** &mdash; *Method*.



```julia
eeg_senv_median(eeg; channel, dims, d, mt)
```

Calculate spectral envelope: median and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over chan (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `mt::Bool=false`: if true use multi-tapered spectrogram
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: median
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroAnalyzer.eeg_ispc-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ispc-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ispc`** &mdash; *Method*.



```julia
eeg_ispc(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate ISPC (Inter-Site-Phase Clustering).

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

Named tuple containing:

  * `ispc::Array{Float64, 2}`: ISPC value
  * `ispc_angle::Array{Float64, 2}`: ISPC angle
  * `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
  * `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
  * `s1_phase::Array{Float64, 3}`: signal 1 phase
  * `s2_phase::Array{Float64, 3}`: signal 2 phase

<a id='NeuroAnalyzer.eeg_itpc-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_itpc-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_itpc`** &mdash; *Method*.



```julia
eeg_itpc(eeg; channel, t, w)
```

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs/trials.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc::Vector{Float64}`: ITPC or wITPC value
  * `itpcz::Vector{Float64}`: Rayleigh's ITPC Z value
  * `itpc_angle::Vector{Float64}`: ITPC angle
  * `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)

<a id='NeuroAnalyzer.eeg_pli-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pli-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pli`** &mdash; *Method*.



```julia
eeg_pli(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate PLI (Phase Lag Index) between `eeg1` and `eeg2`.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

Named tuple containing:

  * `pli::Array{Float64, 2}`: PLI value
  * `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
  * `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
  * `s1_phase::Array{Float64, 3}`: signal 1 phase
  * `s2_phase::Array{Float64, 3}`: signal 2 phase

<a id='NeuroAnalyzer.eeg_pli-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_pli-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_pli`** &mdash; *Method*.



```julia
eeg_pli(eeg; channel)
```

Calculate PLIs (Phase Lag Index).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `pli_m::Array{Float64, 3}`: PLI value matrices over epochs

<a id='NeuroAnalyzer.eeg_ispc-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ispc-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ispc`** &mdash; *Method*.



```julia
eeg_ispc(eeg; channel)
```

Calculate ISPCs (Inter-Site-Phase Clustering).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `ispc_m::Array{Float64, 3}`: ISPC value matrices over epochs

<a id='NeuroAnalyzer.eeg_aec-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_aec-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_aec`** &mdash; *Method*.



```julia
eeg_aec(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate amplitude envelope correlation.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Int64`
  * `channel2::Int64`
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

Named tuple containing:

  * `aec::Vector{Float64}`: power correlation value
  * `aec_p::Vector{Float64}`: power correlation p-value

<a id='NeuroAnalyzer.eeg_ged-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ged-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ged`** &mdash; *Method*.



```julia
eeg_ged(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Perform generalized eigendecomposition..

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`: signal data to be analyzed
  * `eeg2::NeuroAnalyzer.EEG`: original signal data
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

  * `sged::Array{Float64, 3}`
  * `ress::Matrix{Float64}`
  * `ress_normalized::Matrix{Float64}`

<a id='NeuroAnalyzer.eeg_frqinst-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_frqinst-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_frqinst`** &mdash; *Method*.



```julia
eeg_frqinst(eeg; channel)
```

Calculate instantaneous frequency of `eeg`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `frqinst::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_itpc_s-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_itpc_s-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_itpc_s`** &mdash; *Method*.



```julia
eeg_itpc_s(eeg; <keyword arguments>)
```

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Int64`
  * `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
  * `frq_n::Int64`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
  * `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
  * `itpc_frq::Vector{Float64}`: frequencies list

<a id='NeuroAnalyzer.eeg_tkeo-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_tkeo-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_tkeo`** &mdash; *Method*.



```julia
eeg_tkeo(eeg; channel)
```

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

  * `tkeo::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_mwpsd-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_mwpsd-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_mwpsd`** &mdash; *Method*.



```julia
eeg_mwpsd(eeg; channel, pad, norm, frq_lim, frq_n, frq, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `norm::Bool`=true: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
  * `frq_n::Int64=10`: number of frequencies
  * `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq*n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq*n) for frq === :lin

**Returns**

Named tuple containing:

  * `w_pow::Array{Float64, 4}`
  * `w_frq::Matrix{Float64}`

<a id='NeuroAnalyzer.eeg_fcoherence-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_fcoherence-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_fcoherence`** &mdash; *Method*.



```julia
eeg_fcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2, frq_lim)
```

Calculate coherence (mean over frequencies) and MSC (magnitude-squared coherence).

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
  * `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

**Returns**

Named tuple containing:

  * `c::Array{Float64, 3}`: coherence
  * `msc::Array{Float64, 3}`: MSC
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.eeg_vartest-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_vartest-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_vartest`** &mdash; *Method*.



```julia
eeg_vartest(eeg; channel)
```

Calculate variance F-test.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels

**Returns**

Named tuple containing:

  * `f::Array{Float64, 3}`
  * `p::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_vartest-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_vartest-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_vartest`** &mdash; *Method*.



```julia
eeg_vartest(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Calculate variance F-test.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

Named tuple containing:

  * `f::Array{Float64, 3}`
  * `p::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_band_mpower-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_band_mpower-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_band_mpower`** &mdash; *Method*.



```julia
eeg_band_mpower(eeg; channel, f, mt)
```

Calculate mean and maximum band power and its frequency.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `mbp::Matrix{Float64}`: mean band power [μV^2/Hz] per channel per epoch
  * `maxfrq::Matrix{Float64}`: frequency of maximum band power [Hz] per channel per epoch
  * `maxbp::Matrix{Float64}`: power at maximum band frequency [μV^2/Hz] per channel per epoch

<a id='NeuroAnalyzer.eeg_rel_psd-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_rel_psd-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_rel_psd`** &mdash; *Method*.



```julia
eeg_rel_psd(eeg; channel, norm, mt, f)
```

Calculate relative power spectrum density.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram
  * `f::Union{Tuple{Real, Real}, Nothing}=nothing`: calculate power relative to frequency range or total power

**Returns**

Named tuple containing:

  * `psd_pow::Array{Float64, 3}`:powers
  * `psd_frq::Array{Float64, 3}`: frequencies

<a id='NeuroAnalyzer.eeg_fbsplit-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_fbsplit-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_fbsplit`** &mdash; *Method*.



```julia
eeg_fbsplit(eeg; channel, order)
```

Split EEG signal into frequency bands.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `order::Int64=8`: band-pass filter order

**Returns**

Named tuple containing:

  * `band_names::Vector{Symbol}`
  * `band_frq::Vector{Tuple{Real, Real}}`
  * `signal_split::Array{Float64, 4}`

<a id='NeuroAnalyzer.eeg_chdiff-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_chdiff-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_chdiff`** &mdash; *Method*.



```julia
eeg_chdiff(eeg1, eeg2; channel1, channel2, epoch1, epoch2)
```

Subtract channels.

**Arguments**

  * `eeg1::NeuroAnalyzer.EEG`
  * `eeg2::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

**Returns**

  * `ch_diff::Matrix{Float64}`

<a id='NeuroAnalyzer.eeg_cps-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cps-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cps`** &mdash; *Method*.



```julia
eeg_cps(eeg; channel, norm)
```

Calculate cross power spectrum.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `norm::Bool=true`: normalize do dB

**Returns**

Named tuple containing:

  * `cps_pw::Array{Float64, 4}`: cross power spectrum power
  * `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
  * `cps_fq::Vector{Float64}`: cross power spectrum frequencies

<a id='NeuroAnalyzer.eeg_cps-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_cps-Tuple{NeuroAnalyzer.EEG, NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_cps`** &mdash; *Method*.



```julia
eeg_cps(eeg1, eeg2; channel1, channel2, epoch1, epoch2, norm)
```

Calculate cross power spectrum.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
  * `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
  * `norm::Bool=true`: normalize do dB

**Returns**

Named tuple containing:

  * `cps_pw::Vector{Float64}`: cross power spectrum power
  * `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
  * `cps_fq::Vector{Float64}`: cross power spectrum frequencies

<a id='NeuroAnalyzer.eeg_phdiff-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_phdiff-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_phdiff`** &mdash; *Method*.



```julia
eeg_phdiff(eeg; channel, pad, h)
```

Calculate phase difference between EEG channels and mean phase of reference `channel`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of reference channels, default is all EEG/MEG channels except the analyzed one
  * `avg::Symbol=:phase`: method of averaging: `:phase` or `:signal`; for `:signal` `channel` signals are averaged prior to phase calculation; for `:phase` phase is calculated for each reference channel separately and then averaged
  * `pad::Int64=0`: pad signals with 0s
  * `h::Bool=false`: use FFT or Hilbert transformation

**Returns**

  * `ph_diff::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_ampdiff-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_ampdiff-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_ampdiff`** &mdash; *Method*.



```julia
eeg_ampdiff(eeg; channel)
```

Calculate amplitude difference between each `eeg` channel and mean amplitude of `channel`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of reference channels, default is all EEG/MEG channels except the analyzed one

**Returns**

  * `amp_diff::Array{Float64, 3}`

<a id='NeuroAnalyzer.eeg_dwt-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet' href='#NeuroAnalyzer.eeg_dwt-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:DiscreteWavelet'>#</a>
**`NeuroAnalyzer.eeg_dwt`** &mdash; *Method*.



```julia
eeg_dwt(eeg; channel, wt, type, l)
```

Perform discrete wavelet transformation (DWT).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: Stationary Wavelet Transforms (`:sdwt`) or Autocorrelation Wavelet Transforms (`:acdwt`)
  * `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

**Returns**

  * `dwt_c::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)

<a id='NeuroAnalyzer.eeg_cwt-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:CWT' href='#NeuroAnalyzer.eeg_cwt-Union{Tuple{NeuroAnalyzer.EEG}, Tuple{T}} where T<:CWT'>#</a>
**`NeuroAnalyzer.eeg_cwt`** &mdash; *Method*.



```julia
eeg_cwt(eeg; channel, wt)
```

Perform continuous wavelet transformation (CWT).

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

  * `cwt_c::Array{Float64, 4}`: CWT coefficients (by rows)

<a id='NeuroAnalyzer.eeg_psdslope-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_psdslope-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_psdslope`** &mdash; *Method*.



```julia
eeg_psdslope(eeg; channel, f, norm, mt)
```

Calculate PSD linear fit and slope.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `f::Union{Real, Real}=(0, 0)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
  * `norm::Bool=false`: normalize do dB
  * `mt::Bool=false`: if true use multi-tapered periodogram

**Returns**

Named tuple containing:

  * `lf::Array{Float64, 3}`: linear fit for each channel and epoch
  * `psd_slope::Array{Float64, 2}`: slopes of each linear fit
  * `frq::Vector{Float64}`: range of frequencies for the linear fits

<a id='NeuroAnalyzer.eeg_henv-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_henv-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_henv`** &mdash; *Method*.



```julia
eeg_henv(eeg; channel, d)
```

Calculate Hilbert spectrum amplitude envelope of `eeg`.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_henv_mean-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_henv_mean-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_henv_mean`** &mdash; *Method*.



```julia
eeg_henv_mean(eeg; channel, dims, d)
```

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
  * `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  * `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_henv_median-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_henv_median-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_henv_median`** &mdash; *Method*.



```julia
eeg_henv_median(eeg; channel, dims, d)
```

Calculate Hilbert spectrum amplitude envelope of `eeg`: median and 95% CI.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
  * `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  * `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.eeg_apply-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_apply-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_apply`** &mdash; *Method*.



```julia
eeg_apply(eeg; channel, f)
```

Apply any function.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `f::String`: function to be applied, e.g. `f="mean(eeg, dims=3)"; EEG signal is given using variable`eeg` here.

**Returns**

  * `out::Array{Float64, 3}`


<a id='EEG-plots'></a>

<a id='EEG-plots-1'></a>

## EEG plots

<a id='NeuroAnalyzer.eeg_plot_save-Tuple{Union{Figure, Plots.Plot{Plots.GRBackend}}}' href='#NeuroAnalyzer.eeg_plot_save-Tuple{Union{Figure, Plots.Plot{Plots.GRBackend}}}'>#</a>
**`NeuroAnalyzer.eeg_plot_save`** &mdash; *Method*.



```julia
eeg_plot_save(p; file_name::String)
```

Saves plot as file (PDF/PNG/TIFF). File format is determined using `file_name` extension.

**Arguments**

  * `p::Union{Plots.Plot{Plots.GRBackend}, GLMakie.Figure}`
  * `file_name::String`

<a id='NeuroAnalyzer.plot_signal-Tuple{AbstractVector, AbstractArray}' href='#NeuroAnalyzer.plot_signal-Tuple{AbstractVector, AbstractArray}'>#</a>
**`NeuroAnalyzer.plot_signal`** &mdash; *Method*.



```julia
plot_signal(t, signal; <keyword arguments>)
```

Plot amplitude of single- or multi-channel `signal`.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `signal::Union{AbstractVector, AbstractArray}`: data to plot
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `scale::Bool=true`: draw scale
  * `units::String="μV"`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_signal-Tuple{AbstractVector, AbstractArray, Vector{Bool}}' href='#NeuroAnalyzer.plot_signal-Tuple{AbstractVector, AbstractArray, Vector{Bool}}'>#</a>
**`NeuroAnalyzer.plot_signal`** &mdash; *Method*.



```julia
plot_signal(t, signal, bad; <keyword arguments>)
```

Plot amplitude of single- or multi-channel `signal`.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `signal::Union{AbstractVector, AbstractArray}`: data to plot
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `scale::Bool=true`: draw scale
  * `units::String="μV"`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_signal_avg-Tuple{AbstractVector, AbstractArray}' href='#NeuroAnalyzer.plot_signal_avg-Tuple{AbstractVector, AbstractArray}'>#</a>
**`NeuroAnalyzer.plot_signal_avg`** &mdash; *Method*.



```julia
plot_signal_avg(t, signal; <keyword arguments>)
```

Plot amplitude mean and ±95% CI of averaged `signal` channels.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `signal::AbstractArray`: data to plot
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `scale::Bool=true`: draw scale
  * `units::String="μV"`: units of the scale
  * `norm::Bool=false`: normalize to -1 .. +1
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_signal_butterfly-Tuple{AbstractVector, AbstractArray}' href='#NeuroAnalyzer.plot_signal_butterfly-Tuple{AbstractVector, AbstractArray}'>#</a>
**`NeuroAnalyzer.plot_signal_butterfly`** &mdash; *Method*.



```julia
plot_signal_butterfly(t, signal; <keyword arguments>)
```

Butterfly plot of `signal` channels.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `signal::AbstractArray`: data to plot
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `scale::Bool=true`: draw scale
  * `units::String="μV"`: units of the scale
  * `mono::Bool=false`: use color or grey palette
  * `norm::Bool=false`: normalize to -1 .. +1
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot`** &mdash; *Method*.



```julia
eeg_plot(eeg; <keyword arguments>)
```

Plot signal.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: channel(s) to plot, default is all channels
  * `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is no label
  * `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or grey palette
  * `emarkers::Bool`: draw epoch markers if available
  * `markers::Bool`: draw markers if available
  * `scale::Bool=true`: draw scale
  * `units::String="μV"`: units of the scale
  * `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
  * `norm::Bool=false`: normalize signal for butterfly and averaged plots
  * `bad::Union{Bool, Matrix{Bool}}=false`: list of bad channels; if not empty - plot bad channels using this list
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}' href='#NeuroAnalyzer.eeg_plot-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}'>#</a>
**`NeuroAnalyzer.eeg_plot`** &mdash; *Method*.



```julia
eeg_plot(eeg, c; <keyword arguments>)
```

Plot embedded or external component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `epoch::Union{Int64, AbstractRange}=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
  * `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is no label
  * `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or grey palette
  * `emarkers::Bool`: draw epoch markers if available
  * `markers::Bool`: draw markers if available
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
  * `norm::Bool=false`: normalize signal for butterfly and averaged plots
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd-Tuple{Vector{Float64}, Vector{Float64}}' href='#NeuroAnalyzer.plot_psd-Tuple{Vector{Float64}, Vector{Float64}}'>#</a>
**`NeuroAnalyzer.plot_psd`** &mdash; *Method*.



```julia
plot_psd(s_frq, s_pow; <keyword arguments>)
```

Plot PSD (power spectrum density).

**Arguments**

  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Vector{Float64}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for the Y-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_avg-Tuple{Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_psd_avg-Tuple{Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_psd_avg`** &mdash; *Method*.



```julia
plot_psd_avg(s_frq, s_pow; <keyword arguments>)
```

Plot PSD mean and ±95% CI of averaged channels.

**Arguments**

  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 3}`: powers
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_butterfly-Tuple{Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_psd_butterfly-Tuple{Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_psd_butterfly`** &mdash; *Method*.



```julia
plot_psd_butterfly(s_frq, s_pow; <keyword arguments>)
```

Butterfly PSD plot.

**Arguments**

  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 3}`: powers
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_3d-Tuple{Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_psd_3d-Tuple{Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_psd_3d`** &mdash; *Method*.



```julia
plot_psd_w3d(s_frq, s_pow; <keyword arguments>)
```

Plot 3-d waterfall PSD plot.

**Arguments**

  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 3}`: powers
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_topo-Tuple{DataFrame, Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_psd_topo-Tuple{DataFrame, Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_psd_topo`** &mdash; *Method*.



```julia
plot_psd_topo(s_frq, s_pow; <keyword arguments>)
```

Plot topographical map `eeg` PSD. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*theta, loc*radius, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 3}`: powers
  * `Union{Vector{Int64}, AbstractRange}`: which channels to plot
  * `labels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `fig::GLMakie.Figure`

<a id='NeuroAnalyzer.eeg_plot_psd-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_psd-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_psd`** &mdash; *Method*.



```julia
eeg_plot_psd(eeg; <keyword arguments>)
```

Plot power spectrum density.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:welch`: method of calculating PSD: Welch's periodogram, (`:welch`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
  * `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to EEG band: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [μV^2/Hz]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [μV^2/Hz]
  * `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or grey palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend} | GLMakie.Figure`

<a id='NeuroAnalyzer.eeg_plot_psd-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}' href='#NeuroAnalyzer.eeg_plot_psd-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}'>#</a>
**`NeuroAnalyzer.eeg_plot_psd`** &mdash; *Method*.



```julia
eeg_plot_psd(eeg::NeuroAnalyzer.EEG, c::Union{Symbol, AbstractArray}; epoch::Int64, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, norm::Bool=true, method::Symbol=:welch, frq_lim::Tuple{Real, Real}=(0, 0), ncyc::Union{Int64, Tuple{Int64, Int64}}=6, ref::Symbol=:abs, ax::Symbol=:linlin, xlabel::String="default", ylabel::String="default", zlabel::String="default", title::String="default", mono::Bool=false, type::Symbol=:normal, kwargs...)
```

Plot power spectrum density of embedded or external component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:welch`: method of calculating PSD: Welch's periodogram, (`:welch`), multi-tapered periodogram (`:mt`), Morlet wavelet convolution (`:mw`)
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: x-axis limit
  * `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to EEG band: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`), linear-log10 (`:linlog`), log10-log10 (:loglog)
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [μV^2/Hz]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [μV^2/Hz]
  * `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or grey palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_spectrogram-Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_spectrogram-Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_spectrogram`** &mdash; *Method*.



```julia
plot_spectrogram(s_t, s_frq, s_pow; <keyword arguments>)
```

Plot single-channel spectrogram.

**Arguments**

  * `s_t::Vector{Float64}`: time
  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 2}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_spectrogram-Tuple{Vector{String}, Vector{Float64}, Matrix{Float64}}' href='#NeuroAnalyzer.plot_spectrogram-Tuple{Vector{String}, Vector{Float64}, Matrix{Float64}}'>#</a>
**`NeuroAnalyzer.plot_spectrogram`** &mdash; *Method*.



```julia
plot_spectrogram(s_ch, s_frq, s_pow; <keyword arguments>)
```

Plot multiple-channel spectrogram.

**Arguments**

  * `s_ch::Vector{String}`: channel labels
  * `s_frq::Vector{Float64}`: frequencies
  * `s_pow::Array{Float64, 2}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_spectrogram-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_spectrogram-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_spectrogram(eeg; <keyword arguments>)
```

Plots spectrogram.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `epoch::Int64`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:standard`: method of calculating spectrogram:

      * `:standard`: standard
      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
  * `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]

[channel: 1, epoch: 1, time window: 0 ms:10 s]

  * `mono::Bool=false`: use color or grey palette
  * `markers::Bool`: draw markers if available
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_spectrogram-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}' href='#NeuroAnalyzer.eeg_plot_spectrogram-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}'>#</a>
**`NeuroAnalyzer.eeg_plot_spectrogram`** &mdash; *Method*.



```julia
eeg_plot_spectrogram(eeg, c; <keyword arguments>)
```

Plots spectrogram of embedded or external component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `epoch::Int64`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:standard`: method of calculating spectrogram:

      * `:standard`: standard
      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
  * `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]

[component: 1, epoch: 1, time window: 0 ms:10 s]

  * `mono::Bool=false`: use color or grey palette
  * `markers::Bool`: draw markers if available
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_electrodes-Tuple{DataFrame}' href='#NeuroAnalyzer.plot_electrodes-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.plot_electrodes`** &mdash; *Method*.



```julia
plot_electrodes(locs; <keyword arguments>)
```

Preview of electrode locations. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*theta, loc*radius, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
  * `labels::Bool=true`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette
  * `head_details::Bool=true`: draw nose and ears
  * `plot_size::Int64=400`: plot dimensions in pixels (size × size)

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_electrodes3d-Tuple{DataFrame}' href='#NeuroAnalyzer.plot_electrodes3d-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.plot_electrodes3d`** &mdash; *Method*.



```julia
plot_electrodes3d(locs; <keyword arguments>)
```

3D interactive preview of electrode locations. It uses spherical :loc*x, :loc*y and :loc_z locations.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*theta, loc*radius, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
  * `labels::Bool=true`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette

`plot_size::Int64=800`: plot dimensions in pixels (plot*size×plot*size)

**Returns**

  * `fig::GLMakie.Figure`

<a id='NeuroAnalyzer.eeg_plot_electrodes-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_electrodes-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_electrodes`** &mdash; *Method*.



```julia
eeg_plot_electrodes(eeg; <keyword arguments>)
```

Preview of electrode locations.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: which channel should be highlighted
  * `labels::Bool=true`: plot electrode labels
  * `head::Bool`=true: plot head
  * `head_labels::Bool=false`: plot head labels
  * `plot_size::Int64=400`: plot dimensions in pixels (plot*size×plot*size)
  * `head_details::Bool=true`: draw nose and ears
  * `mono::Bool=false`: use color or grey palette
  * `threed::Bool=false`: 3-dimensional plot
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_matrix-Tuple{Matrix{<:Real}}' href='#NeuroAnalyzer.plot_matrix-Tuple{Matrix{<:Real}}'>#</a>
**`NeuroAnalyzer.plot_matrix`** &mdash; *Method*.



```julia
plot_matrix(m; <keyword arguments>)
```

Plot matrix.

**Arguments**

  * `m::Array{<:Real, 2}`
  * `xlabels::Vector{String}`
  * `ylabels::Vector{String}`
  * `xlabel::String=""`
  * `ylabel::String=""`
  * `title::String=""`
  * `cb_title::String=""`: color bar title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_covmatrix-Tuple{AbstractVector, AbstractVector}' href='#NeuroAnalyzer.plot_covmatrix-Tuple{AbstractVector, AbstractVector}'>#</a>
**`NeuroAnalyzer.plot_covmatrix`** &mdash; *Method*.



```julia
plot_covmatrix(m, lags; <keyword arguments>)
```

Plot cross/auto-covariance matrix.

**Arguments**

  * `m::Abstractvector`: covariance matrix
  * `lags::AbstractVector`: covariance lags, lags will be displayed in ms
  * `xlabel::String="lag"`
  * `ylabel::String=""`
  * `title::String=""`
  * `cb_title::String=""`: color bar title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_histogram-Tuple{AbstractVector}' href='#NeuroAnalyzer.plot_histogram-Tuple{AbstractVector}'>#</a>
**`NeuroAnalyzer.plot_histogram`** &mdash; *Method*.



```julia
plot_histogram(signal; <keyword arguments>)
```

Plot histogram of `signal`.

**Arguments**

  * `signal::AbstractVector`
  * `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
  * `label::String=""`: channel label
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_filter_response-Tuple{}' href='#NeuroAnalyzer.plot_filter_response-Tuple{}'>#</a>
**`NeuroAnalyzer.plot_filter_response`** &mdash; *Method*.



```julia
plot_filter_response(<keyword arguments>)
```

Plot filter response.

**Arguments**

  * `fs::Int64`: sampling rate
  * `fprototype::Symbol`: filter class: `:fir`, `:butterworth`, `:chebyshev1`, `:chebyshev2`, `:elliptic`
  * `ftype::Symbol`: filter type: low-pass (`:lp`), high-pass (`:hp`), band-pass (`:bp`), band-stop (`:bs`)
  * `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
  * `order::Int64`: filter order
  * `rp::Real`: dB ripple in the passband
  * `rs::Real`: dB attenuation in the stopband
  * `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez filters
  * `window::window::Union{Vector{Float64}, Nothing}`: window, required for FIR filter
  * `mono::Bool=false`: use color or grey palette
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_weights-Tuple{DataFrame}' href='#NeuroAnalyzer.plot_weights-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.plot_weights`** &mdash; *Method*.



```julia
plot_weights(locs; <keyword arguments>)
```

Plot weights at electrode positions. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*theta, loc*radius, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
  * `weights::Vector{<:Real}=[]`: weights vector
  * `labels::Bool=true`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette
  * `head_details::Bool=true`: draw nose and ears
  * `plot_size::Int64=400`: plot dimensions in pixels (size × size)

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_weights-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_weights-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_weights`** &mdash; *Method*.



```julia
eeg_plot_weights(eeg; <keyword arguments>)
```

Plot weights at electrode positions. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `weights::Matrix{<:Real}`: matrix of weights
  * `labels::Bool=false`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette
  * `head_details::Bool=true`: draw nose and ears
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_connections-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_connections-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_connections`** &mdash; *Method*.



```julia
eeg_plot_connections(eeg; <keyword arguments>)
```

Plot weights at electrode positions. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `connections::Matrix{<:Real}`: matrix of connections weights
  * `threshold::Real`: plot all connection above threshold
  * `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
  * `weights::Bool=true`: weight line widths and alpha based on connection value
  * `labels::Bool=false`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette
  * `head_details::Bool=true`: draw nose and ears
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_connections-Tuple{DataFrame}' href='#NeuroAnalyzer.plot_connections-Tuple{DataFrame}'>#</a>
**`NeuroAnalyzer.plot_connections`** &mdash; *Method*.



```julia
eeg_plot_connections(eeg; <keyword arguments>)
```

Plot connections between channels. It uses polar :loc*radius and :loc*theta locations, which are translated into Cartesian x and y positions.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`
  * `channel::Union{Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `connections::Matrix{<:Real}`: matrix of connections weights
  * `threshold::Real`: plot all connection above threshold
  * `threshold_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)
  * `weights::Bool=true`: weight line widths and alpha based on connection value
  * `labels::Bool=false`: plot electrode labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or grey palette
  * `head_details::Bool=true`: draw nose and ears
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_stats-Tuple{NeuroAnalyzer.EEG, Union{Symbol, Dict, VecOrMat{<:Real}}}' href='#NeuroAnalyzer.eeg_plot_stats-Tuple{NeuroAnalyzer.EEG, Union{Symbol, Dict, VecOrMat{<:Real}}}'>#</a>
**`NeuroAnalyzer.eeg_plot_stats`** &mdash; *Method*.



```julia
eeg_plot(eeg, c; <keyword arguments>)
```

Plot channel/epoch data.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `c::Union{Vector{<:Real}, Matrix{<:Real}, Symbol, Dict}`: component to plot
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: epoch to display
  * `epoch::Union{Int64, Vector{Int64}, AbstractRange}=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: channel to display
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `plot_by::Symbol`: c values refer to: :labels, :channels or :epochs
  * `type::Symbol`: plot type:

      * `:hist`: histogram
      * `:kd`: kernel density
      * `:bar`: bar plot
      * `:heat`: heatmap; rows of `c` must be time points, columns must be either channels or epochs
      * `:box`: box plot; `c` must contain ≥ 2 values per group
      * `:violin`: violin plot; `c` must contain ≥ 2 values per group
      * `:paired`: paired
      * `:polar`: polar
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

**Notes**

Labeled matrix is a dictionary of labels and vectors associated with these labels. This way plotting of vectors non-equal length is possible. Labeled matrix is created using `_labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{<:Real}})` and converted back to keys and values using `_dict2labeled_matrix(d::Dict)`.

For `:polar plot` if `c` is a vector, than it contains phases in radians. If `c` is a two column matrix, than first column contains phases in radians and second column contains lengths.

<a id='NeuroAnalyzer.plot_topo-Tuple{Vector{<:Real}}' href='#NeuroAnalyzer.plot_topo-Tuple{Vector{<:Real}}'>#</a>
**`NeuroAnalyzer.plot_topo`** &mdash; *Method*.



```julia
plot_topo(c; <keyword arguments>)
```

Plot topographical view.

**Arguments**

  * `signal::Vector{<:Real}`: values to plot (one value per channel)
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
  * `locs::DataFrame`: columns: channel, labels, loc*theta, loc*radius, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or grey palette
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `head_labels::Bool=false`: plot head labels
  * `head_details::Bool=true`: draw nose and ears
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_topo-Tuple{NeuroAnalyzer.EEG}' href='#NeuroAnalyzer.eeg_plot_topo-Tuple{NeuroAnalyzer.EEG}'>#</a>
**`NeuroAnalyzer.eeg_plot_topo`** &mdash; *Method*.



```julia
eeg_plot_topo(eeg; <keyword arguments>)
```

Topographical plot.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `epoch::Union{Int64, AbstractRange}=0`: epoch to display
  * `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
  * `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
  * `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or grey palette
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method: `:mean`, `:median`
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `head_labels::Bool=false`: plot head labels
  * `head_details::Bool=true`: draw nose and ears
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_topo-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}' href='#NeuroAnalyzer.eeg_plot_topo-Tuple{NeuroAnalyzer.EEG, Union{Symbol, AbstractArray}}'>#</a>
**`NeuroAnalyzer.eeg_plot_topo`** &mdash; *Method*.



```julia
eeg_plot_topo(eeg; <keyword arguments>)
```

Topographical plot of embedded or external component.

**Arguments**

  * `eeg::NeuroAnalyzer.EEG`: EEG object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `epoch::Union{Int64, AbstractRange}=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
  * `segment::Tuple{Int64, Int64}=(1, 10*eeg_sr(eeg))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
  * `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or grey palette
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method: `:mean`, `:median`
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `nmethod::Symbol=:minmax`: method for normalization, see `s_normalize()`
  * `plot_size::Int64=800`: plot dimensions in pixels (size × size)
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `head_labels::Bool=false`: plot head labels
  * `head_details::Bool=true`: draw nose and ears
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}' href='#NeuroAnalyzer.eeg_plot_compose-Tuple{Vector{Plots.Plot{Plots.GRBackend}}}'>#</a>
**`NeuroAnalyzer.eeg_plot_compose`** &mdash; *Method*.



```julia
eeg_plot_compose(p; <keyword arguments>)
```

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:

  * `(2, 2)`: 2 × 2 plots, regular layout
  * `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

**Arguments**

  * `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
  * `layout::Union(Matrix{Any}, Tuple{Int64, Int64}}`: layout
  * `mono::Bool=false`: use color or grey palette
  * `title::String=""`: plot title
  * `kwargs`: optional arguments for `p` vector plots

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_empty-Tuple{}' href='#NeuroAnalyzer.plot_empty-Tuple{}'>#</a>
**`NeuroAnalyzer.plot_empty`** &mdash; *Method*.



```julia
plot_empty()
```

Return an empty plot, useful for filling matrices of plots. 

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


<a id='Neurostimulation'></a>

<a id='Neurostimulation-1'></a>

## Neurostimulation

<a id='NeuroAnalyzer.tes_dose-Tuple{Real, Real, Int64}' href='#NeuroAnalyzer.tes_dose-Tuple{Real, Real, Int64}'>#</a>
**`NeuroAnalyzer.tes_dose`** &mdash; *Method*.



```julia
tes_dose(current, pad_area, duration)
```

Converts `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

**Arguments**

  * `current::Real`: stimulation current [mA]
  * `pad_area::Real`: electrode pad area [cm^2]
  * `duration::Int64`: stimulation duration [s]

**Returns**

  * `charge::Float64`: charge [C]
  * `current_density::Float64`: current density [A/m^2]
  * `charge_density::Float64`: delibvered charge density [kC/m^2]

**Source**

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.


<a id='NIRS'></a>

<a id='NIRS-1'></a>

## NIRS


<a id='MRI'></a>

<a id='MRI-1'></a>

## MRI


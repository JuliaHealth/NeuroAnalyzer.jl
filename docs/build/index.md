


![](assets/neuroanalyzer.png)


<a id='NeuroAnalyzer.jl-documentation'></a>

<a id='NeuroAnalyzer.jl-documentation-1'></a>

# NeuroAnalyzer.jl documentation


This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).


<a id='NeuroAnalyzer'></a>

<a id='NeuroAnalyzer-1'></a>

## NeuroAnalyzer

<a id='NeuroAnalyzer.na_info' href='#NeuroAnalyzer.na_info'>#</a>
**`NeuroAnalyzer.na_info`** &mdash; *Function*.



```julia
na_info()
```

Show NeuroAnalyzer and imported packages versions.

<a id='NeuroAnalyzer.na_plugins_install' href='#NeuroAnalyzer.na_plugins_install'>#</a>
**`NeuroAnalyzer.na_plugins_install`** &mdash; *Function*.



```julia
na_plugins_install(plugin)
```

Install NeuroAnalyzer plugin from remote Git repository.

**Arguments**

  * `plugin::String`: plugin Git repository URL

<a id='NeuroAnalyzer.na_plugins_list' href='#NeuroAnalyzer.na_plugins_list'>#</a>
**`NeuroAnalyzer.na_plugins_list`** &mdash; *Function*.



```julia
na_plugins_list()
```

List NeuroAnalyzer plugins.

<a id='NeuroAnalyzer.na_plugins_reload' href='#NeuroAnalyzer.na_plugins_reload'>#</a>
**`NeuroAnalyzer.na_plugins_reload`** &mdash; *Function*.



```julia
na_plugins_reload()
```

Reload NeuroAnalyzer plugins.

<a id='NeuroAnalyzer.na_plugins_remove' href='#NeuroAnalyzer.na_plugins_remove'>#</a>
**`NeuroAnalyzer.na_plugins_remove`** &mdash; *Function*.



```julia
na_plugins_remove(plugin)
```

Remove NeuroAnalyzer plugin.

**Arguments**

  * `plugin::String`: plugin name

<a id='NeuroAnalyzer.na_plugins_update' href='#NeuroAnalyzer.na_plugins_update'>#</a>
**`NeuroAnalyzer.na_plugins_update`** &mdash; *Function*.



```julia
na_plugins_update(plugin)
```

Install NeuroAnalyzer plugin.

**Arguments**

  * `plugin::String`: plugin to update; if empty, update all

<a id='NeuroAnalyzer.na_set_prefs' href='#NeuroAnalyzer.na_set_prefs'>#</a>
**`NeuroAnalyzer.na_set_prefs`** &mdash; *Function*.



```julia
na_set_prefs(use_cuda, progress_bar, verbose)
```

Save NeuroAnalyzer preferences.

**Arguments**

  * `use_cuda::Bool`
  * `progress_bar::Bool`
  * `verbose::Bool`

<a id='NeuroAnalyzer.na_set_progress_bar' href='#NeuroAnalyzer.na_set_progress_bar'>#</a>
**`NeuroAnalyzer.na_set_progress_bar`** &mdash; *Function*.



```julia
na_set_progress_bar(value)
```

Change `progress_bar` preference.

**Arguments**

  * `value::Bool`: value

<a id='NeuroAnalyzer.na_set_use_cuda' href='#NeuroAnalyzer.na_set_use_cuda'>#</a>
**`NeuroAnalyzer.na_set_use_cuda`** &mdash; *Function*.



```julia
na_set_use_cuda(value)
```

Change `use_cuda` preference.

**Arguments**

  * `value::Bool`: value

<a id='NeuroAnalyzer.na_set_verbose' href='#NeuroAnalyzer.na_set_verbose'>#</a>
**`NeuroAnalyzer.na_set_verbose`** &mdash; *Function*.



```julia
na_set_verbose(value)
```

Change `verbose` preference.

**Arguments**

  * `value::Bool`: value

<a id='NeuroAnalyzer.na_version' href='#NeuroAnalyzer.na_version'>#</a>
**`NeuroAnalyzer.na_version`** &mdash; *Function*.



```julia
na_version()
```

Convert NeuroAnalyzer version to string.

**Returns**

  * `VER::String`


<a id='Utils'></a>

<a id='Utils-1'></a>

## Utils

<a id='NeuroAnalyzer.add_component' href='#NeuroAnalyzer.add_component'>#</a>
**`NeuroAnalyzer.add_component`** &mdash; *Function*.



```julia
add_component(obj; c, v)
```

Add component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name
  * `v::Any`: component value

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_component!' href='#NeuroAnalyzer.add_component!'>#</a>
**`NeuroAnalyzer.add_component!`** &mdash; *Function*.



```julia
add_component!(obj; c, v)
```

Add component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name
  * `v::Any`: component value

<a id='NeuroAnalyzer.add_note' href='#NeuroAnalyzer.add_note'>#</a>
**`NeuroAnalyzer.add_note`** &mdash; *Function*.



```julia
add_note(obj; note)
```

Add recording note to the object header.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `note::String`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_note!' href='#NeuroAnalyzer.add_note!'>#</a>
**`NeuroAnalyzer.add_note!`** &mdash; *Function*.



```julia
add_note!(obj; note)
```

Add recording note to the object header.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `note::String`

<a id='NeuroAnalyzer.apply' href='#NeuroAnalyzer.apply'>#</a>
**`NeuroAnalyzer.apply`** &mdash; *Function*.



```julia
apply(obj; ch, f)
```

Apply custom function.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `f::String`: function to be applied, e.g. `f="mean(obj, dims=3)"`; OBJ signal is given using variable `obj` here.

**Returns**

  * `out::Array{Float64, 3}`

<a id='NeuroAnalyzer.areduce' href='#NeuroAnalyzer.areduce'>#</a>
**`NeuroAnalyzer.areduce`** &mdash; *Function*.



```julia
areduce(a, f; n)
```

Reduce an array at indices of a vector being multiplications of a constant. Useful e.g. for simplifying values across frequencies, when the number of frequencies (and thus values) is high.

**Arguments**

  * `a::AbstractArray`: e.g. signal data
  * `f::AbstractVector`: e.g. frequencies
  * `n::Float64=0.5`: reduce at multiplications of this value

**Returns**

  * `a_new::Array{eltype(a), ndims(a)}`
  * `f_new::Vector{eltype(f)}`

<a id='NeuroAnalyzer.arr2mat' href='#NeuroAnalyzer.arr2mat'>#</a>
**`NeuroAnalyzer.arr2mat`** &mdash; *Function*.



```julia
arr2mat(x)
```

Reshape array into matrix.

**Arguments**

  * `x::AbstractArray`

**Returns**

  * `m::Matrix{eltype(x)}`

<a id='NeuroAnalyzer.band_frq' href='#NeuroAnalyzer.band_frq'>#</a>
**`NeuroAnalyzer.band_frq`** &mdash; *Function*.



```julia
band_frq(obj, band)
```

Return band frequency limits.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `band::Symbol`: band range name:

      * `:list`
      * `:total`
      * `:delta`: 0.1 - 4.0 Hz
      * `:theta`: 4.0 - 8.0 Hz
      * `:alpha`: 8.0 - 13.0 Hz
      * `:alpha_lower`: 8.0 - 10.5 Hz
      * `:alpha_higher`: 10.5 - 13.0 Hz
      * `:beta`: 14.0 - 30.0 Hz
      * `:beta_lower`: 14.0 - 25.0 Hz
      * `:beta_higher`: 25.0 - 30.0 Hz
      * `:gamma`: 30.0 - 150.0 Hz
      * `:gamma_1`: 30.0 - 40.0 Hz
      * `:gamma_2`: 40.0 - 50.0 Hz
      * `:gamma_lower`: 30.0 - 80.0 Hz
      * `:gamma_higher`: 80.0 - 150.0 Hz

**Returns**

  * `band_frequency::Tuple{Real, Real}`


```
band_frq(fs, band)
```

Return band frequency limits.

**Arguments**

  * `fs::Int64`: sampling rate
  * `band::Symbol`: band range name:

      * `:list`
      * `:total`
      * `:delta`: 0.1 - 4.0 Hz
      * `:theta`: 4.0 - 8.0 Hz
      * `:alpha`: 8.0 - 13.0 Hz
      * `:alpha_lower`: 8.0 - 10.5 Hz
      * `:alpha_higher`: 10.5 - 13.0 Hz
      * `:beta`: 14.0 - 30.0 Hz
      * `:beta_lower`: 14.0 - 25.0 Hz
      * `:beta_higher`: 25.0 - 30.0 Hz
      * `:gamma`: 30.0 - 150.0 Hz
      * `:gamma_1`: 30.0 - 40.0 Hz
      * `:gamma_2`: 40.0 - 50.0 Hz
      * `:gamma_lower`: 30.0 - 80.0 Hz
      * `:gamma_higher`: 80.0 - 150.0 Hz

**Returns**

  * `band_frq::Tuple{Real, Real}`

<a id='Base.size' href='#Base.size'>#</a>
**`Base.size`** &mdash; *Function*.



```julia
size(A::AbstractArray, [dim])
```

Return a tuple containing the dimensions of `A`. Optionally you can specify a dimension to just get the length of that dimension.

Note that `size` may not be defined for arrays with non-standard indices, in which case [`axes`](@ref) may be useful. See the manual chapter on [arrays with custom indices](@ref man-custom-indices).

See also: [`length`](@ref), [`ndims`](@ref), [`eachindex`](@ref), [`sizeof`](@ref).

**Examples**

```julia-repl
julia> A = fill(1, (2,3,4));

julia> size(A)
(2, 3, 4)

julia> size(A, 2)
3
```


<a target='_blank' href='https://github.com/JuliaLang/julia/blob/bd47eca2c8aacd145b6c5c02e47e2b9ec27ab456/base/abstractarray.jl#L20-L41' class='documenter-source'>source</a><br>


```
size(cb::CircularBuffer)
```

Return a tuple with the size of the buffer.


```
size(df::AbstractDataFrame[, dim])
```

Return a tuple containing the number of rows and columns of `df`. Optionally a dimension `dim` can be specified, where `1` corresponds to rows and `2` corresponds to columns.

See also: [`nrow`](@ref), [`ncol`](@ref)

**Examples**

```julia-repl
julia> df = DataFrame(a=1:3, b='a':'c');

julia> size(df)
(3, 2)

julia> size(df, 1)
3
```


```
size(dfr::DataFrameRow[, dim])
```

Return a 1-tuple containing the number of elements of `dfr`. If an optional dimension `dim` is specified, it must be `1`, and the number of elements is returned directly as a number.

See also: [`length`](@ref)

**Examples**

```julia-repl
julia> dfr = DataFrame(a=1:3, b='a':'c')[1, :]
DataFrameRow
 Row │ a      b
     │ Int64  Char
─────┼─────────────
   1 │     1  a

julia> size(dfr)
(2,)

julia> size(dfr, 1)
2
```


```
size(p::Plan, [dim])
```

Return the size of the input of a plan `p`, optionally at a specified dimenion `dim`.


```
size(::AbstractPolynomial, [i])
```

Returns the size of the polynomials coefficients, along axis `i` if provided.


```
size(model::AbstractDimensionalityReduction, d::Int)
```

Returns the dimension of the input data if `d == 1`, the dimension of the output data if `d == 2`, otherwise throws error.


```
size(f)
```

Dimensions of the coefficient matrix of the whitening transform `f`.


```
size(M)
```

Returns a tuple with the dimensions of input (the dimension of the observation space) and output (the dimension of the principal subspace).


```
size(M::PPCA)
```

Returns a tuple with values of the input dimension $d$, *i.e* the dimension of the observation space, and the output dimension $p$, *i.e* the dimension of the principal subspace.


```
size(M::KernelPCA)
```

Returns a tuple with the input dimension $d$, *i.e* the dimension of the observation space, and the output dimension $p$, *i.e* the dimension of the principal subspace.


```
size(M:CCA)
```

Return a tuple with the dimension of `X`, `Y`, and the output dimension.


```
size(M::MDS)
```

Returns tuple where the first value is the MDS model `M` input dimension, *i.e* the dimension of the observation space, and the second value is the output dimension, *i.e* the dimension of the embedding.


```
size(M::MetricMDS)
```

Returns tuple where the first value is the MDS model `M` input dimension, *i.e* the dimension of the observation space, and the second value is the output dimension, *i.e* the dimension of the embedding.


```
size(M::MulticlassLDA)
```

Get the input (*i.e* the dimension of the observation space) and output (*i.e* the dimension of the transformed features) dimensions of the model `M`.


```
size(M)
```

Get the input (*i.e* the dimension of the observation space) and output (*i.e* the dimension of the subspace projection) dimensions of the model `M`.


```
size(M::ICA)
```

Returns a tuple with the input dimension, *i.e* the number of observed mixtures, and the output dimension, *i.e* the number of independent components.


```
size(M::PPCA)
```

Returns a tuple with values of the input dimension $d$, *i.e* the dimension of the observation space, and the output dimension $p$, *i.e* the dimension of the principal subspace.


```
size(s::Sampleable)
```

The size (i.e. shape) of each sample. Always returns `()` when `s` is univariate, and `(length(s),)` when `s` is multivariate.


```
size(d::MultivariateDistribution)
```

Return the sample size of distribution `d`, *i.e* `(length(d),)`.


```
size(d::MatrixDistribution)
```

Return the size of each sample from distribution `d`.


```
size(g, i)
```

Return the number of vertices in `g` if `i`=1 or `i`=2, or `1` otherwise.

**Examples**

```julia-repl
julia> using Graphs

julia> g = cycle_graph(4);

julia> size(g, 1)
4

julia> size(g, 2)
4

julia> size(g, 3)
1
```


```
size(m::MixedModel)
```

Returns the size of a mixed model as a tuple of length four: the number of observations, the number of (non-singular) fixed-effects parameters, the number of conditional modes (random effects), the number of grouping variables


```
size(obj)
```

Return size of the object data.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `size::Tuple{Int64, Int64, Int64}`


```
size(obj, n)
```

Return size of the object data.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64`

**Returns**

  * `size::Int64`

<a id='NeuroAnalyzer.cextrema' href='#NeuroAnalyzer.cextrema'>#</a>
**`NeuroAnalyzer.cextrema`** &mdash; *Function*.



```julia
cextrema(x)
```

Return extreme values of the complex vector.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

Tuple containing:

  * `cmax::ComplexF64`
  * `cmin::ComplexF64`

<a id='NeuroAnalyzer.channel_cluster' href='#NeuroAnalyzer.channel_cluster'>#</a>
**`NeuroAnalyzer.channel_cluster`** &mdash; *Function*.



```julia
channels_cluster(obj, cluster)
```

Return channels belonging to a cluster of channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `cluster::Symbol`: available clusters are:

      * `:f1`: left frontal (F1, F3, F5, F7, AF3, AF7)
      * `:f2`: right frontal (F2, F4, F6, F8, AF4, AF8)
      * `:t1`: left temporal (C3, C5, T7, FC3, FC5, FT7)
      * `:t2`: right temporal (C4, C6, T8, FC4, FC6, FT8)
      * `:c1`: anterior central (Cz, C1, C2, FC1, FC2, FCz)
      * `:c2`: posterior central (Pz, P1, P2, CP1, CP2, CPz)
      * `:p1`: left parietal (P3, P5, P7, CP3, CP5, TP7)
      * `:p2`: right parietal (P4, P6, P8, CP4, CP6, TP8)
      * `:o`: occipital (Oz, O1, O2, POz, PO3, PO4)

**Returns**

  * `ch::Vector{Int64}`: list of channel numbers belonging to a given cluster of channels

<a id='NeuroAnalyzer.channel_pick' href='#NeuroAnalyzer.channel_pick'>#</a>
**`NeuroAnalyzer.channel_pick`** &mdash; *Function*.



```julia
channel_pick(obj; p)
```

Return set of channel indices corresponding to a set of electrodes ("pick", e.g. left or frontal electrodes).

**Arguments**

  * `p::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`

      * `:list`
      * `:central` (or `:c`)
      * `:left` (or `:l`)
      * `:right` (or `:r`)
      * `:frontal` (or `:f`)
      * `:temporal` (or `:t`)
      * `:parietal` (or `:p`)
      * `:occipital` (or `:o`)

**Returns**

  * `channels::Vector{Int64}`: channel numbers

<a id='NeuroAnalyzer.chtypes' href='#NeuroAnalyzer.chtypes'>#</a>
**`NeuroAnalyzer.chtypes`** &mdash; *Function*.



```julia
chtypes(obj)
```

Return channel chtypes.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `chtypes::Vector{String}`

<a id='NeuroAnalyzer.cmax' href='#NeuroAnalyzer.cmax'>#</a>
**`NeuroAnalyzer.cmax`** &mdash; *Function*.



```julia
cmax(x)
```

Return maximum value of the complex vector.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmax::ComplexF64`

<a id='NeuroAnalyzer.cmin' href='#NeuroAnalyzer.cmin'>#</a>
**`NeuroAnalyzer.cmin`** &mdash; *Function*.



```julia
cmin(x)
```

Return minimum value of the complex vector.

**Arguments**

  * `x::Vector{ComplexF64}`

**Returns**

  * `cmin::ComplexF64`

<a id='NeuroAnalyzer.component_type' href='#NeuroAnalyzer.component_type'>#</a>
**`NeuroAnalyzer.component_type`** &mdash; *Function*.



```julia
component_type(obj, c)
```

Return component data type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name

**Return**

  * `c_type::DataType`

<a id='NeuroAnalyzer.cums' href='#NeuroAnalyzer.cums'>#</a>
**`NeuroAnalyzer.cums`** &mdash; *Function*.



```julia
cums(signal)
```

Calculate cumulative sum of a 3-dimensional array.

**Arguments**

  * `signal::Array{<:Real, 3}`

**Returns**

  * `signal_cs::Array{Float64, 3}`

<a id='NeuroAnalyzer.datatype' href='#NeuroAnalyzer.datatype'>#</a>
**`NeuroAnalyzer.datatype`** &mdash; *Function*.



```julia
datatype(obj)
```

Return data type of the object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `data_type::String`

<a id='NeuroAnalyzer.delete_component' href='#NeuroAnalyzer.delete_component'>#</a>
**`NeuroAnalyzer.delete_component`** &mdash; *Function*.



```julia
delete_component(obj; c)
```

Delete component. 

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_component!' href='#NeuroAnalyzer.delete_component!'>#</a>
**`NeuroAnalyzer.delete_component!`** &mdash; *Function*.



```julia
delete_component!(obj; c)
```

Delete component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name

<a id='NeuroAnalyzer.delete_note' href='#NeuroAnalyzer.delete_note'>#</a>
**`NeuroAnalyzer.delete_note`** &mdash; *Function*.



```julia
delete_note(obj)
```

Delete recording note from the object header.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_note!' href='#NeuroAnalyzer.delete_note!'>#</a>
**`NeuroAnalyzer.delete_note!`** &mdash; *Function*.



```julia
delete_note!(obj)
```

Delete recording note from the object header.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delmean' href='#NeuroAnalyzer.delmean'>#</a>
**`NeuroAnalyzer.delmean`** &mdash; *Function*.



```julia
delmean(s)
```

Demean signal.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::AbstractArray`

<a id='NeuroAnalyzer.describe' href='#NeuroAnalyzer.describe'>#</a>
**`NeuroAnalyzer.describe`** &mdash; *Function*.



```julia
describe(obj)
```

Return basic descriptive statistics of the object data.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.detector_labels' href='#NeuroAnalyzer.detector_labels'>#</a>
**`NeuroAnalyzer.detector_labels`** &mdash; *Function*.



```julia
detector_labels(obj)
```

Return NIRS detector labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroAnalyzer.epoch_len' href='#NeuroAnalyzer.epoch_len'>#</a>
**`NeuroAnalyzer.epoch_len`** &mdash; *Function*.



```julia
epoch_len(obj)
```

Return epoch length.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `ep_len::Int64`


```
epoch_len(study)
```

Return length of epochs per NeuroAnalyzer NEURO object in the study.

**Arguments**

  * `study::NeuroAnalyzer.STUDY`

**Returns**

  * `len::Int64`

<a id='NeuroAnalyzer.extract_component' href='#NeuroAnalyzer.extract_component'>#</a>
**`NeuroAnalyzer.extract_component`** &mdash; *Function*.



```julia
extract_component(obj, c)
```

Extract component values.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Symbol`: component name

**Returns**

  * `c::Any`

<a id='NeuroAnalyzer.f2t' href='#NeuroAnalyzer.f2t'>#</a>
**`NeuroAnalyzer.f2t`** &mdash; *Function*.



```julia
f2t(f)
```

Convert frequency in Hz to cycle length in ms.

**Arguments**

  * `f::Real`: frequency in Hz

**Returns**

  * `f::Float64`: cycle length in ms

<a id='NeuroAnalyzer.fft0' href='#NeuroAnalyzer.fft0'>#</a>
**`NeuroAnalyzer.fft0`** &mdash; *Function*.



```julia
fft0(x, n)
```

Perform zeros-padded FFT.

**Arguments**

  * `x::AbstractVector`
  * `n::Int64`: number of zeros to add

**Returns**

  * `fft0::Vector{ComplexF64}`

<a id='NeuroAnalyzer.fft2' href='#NeuroAnalyzer.fft2'>#</a>
**`NeuroAnalyzer.fft2`** &mdash; *Function*.



```julia
fft2(x)
```

Perform zeros-padded FFT, so the length of padded vector is a power of 2.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `fft2::Vector{ComplexF64}`

<a id='NeuroAnalyzer.fft_transform' href='#NeuroAnalyzer.fft_transform'>#</a>
**`NeuroAnalyzer.fft_transform`** &mdash; *Function*.



```julia
fft_transform(x; fs, wlen, woverlap, w, demean, pad, mode)
```

Perform FFT transformation.

**Arguments**

  * `x::AbstractVector`
  * `fs::Int64`: sampling rate
  * `wlen::Int64=fs`: window length
  * `woverlap::Int64=round(Int64, wlen * 0.97)`:
  * `w::Bool=false`: if true, apply Hanning window per segment
  * `demean::Bool=false`: if true, demean each segment
  * `nfft::Int64=0`: length of input vector to the FFT; if nfft > n_samples, then the input signal will be zero-padded until it is of length nfft
  * `mode::Symbol=:r`:

      * `:r`: use one-sided FFT (rfft)
      * `:f`: use two-sided FFT (fft)

**Returns**

  * `mf::Vector{ComplexF64}`: Fourier coefficients
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.findpeaks' href='#NeuroAnalyzer.findpeaks'>#</a>
**`NeuroAnalyzer.findpeaks`** &mdash; *Function*.



```julia
findpeaks(signal; d)
```

Find peaks.

**Arguments**

  * `signal::AbstractVector`
  * `d::Int64=32`: distance between peeks in samples

**Returns**

  * `p_idx::Vector{Int64}`

<a id='NeuroAnalyzer.f_nearest' href='#NeuroAnalyzer.f_nearest'>#</a>
**`NeuroAnalyzer.f_nearest`** &mdash; *Function*.



```julia
f_nearest(m, pos)
```

Find nearest position tuple in a matrix of positions.

**Arguments**

  * `m::Matrix{Tuple{Float64, Float64}}`: matrix of positions
  * `p::Tuple{Float64, Float64}`: position tuple

**Returns**

  * `pos::Tuple{Int64, Int64}`: row and column in m

<a id='NeuroAnalyzer.freqs' href='#NeuroAnalyzer.freqs'>#</a>
**`NeuroAnalyzer.freqs`** &mdash; *Function*.



```julia
freqs(t)
```

Return vector of frequencies and Nyquist frequency for time vector.

**Arguments**

  * `t::AbstractVector, AbstractRange}`: time vector

**Returns**

  * `hz::Vector{Float64}`
  * `nf::Float64`


```
freqs(s, fs)
```

Return vector of frequencies and Nyquist frequency for signal.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`

**Returns**

  * `hz::Vector{Float64`: signal vector
  * `nf::Float64`


```
freqs(obj)
```

Return vector of frequencies and Nyquist frequency.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

Named tuple containing:

  * `hz::Vector{Float64}`
  * `nf::Float64`

<a id='NeuroAnalyzer.generate_csine' href='#NeuroAnalyzer.generate_csine'>#</a>
**`NeuroAnalyzer.generate_csine`** &mdash; *Function*.



```julia
generate_csine(f, t, a)
```

Generates complex sine wave.

**Arguments**

  * `f::Real`: frequency [Hz]
  * `t::Union{AbstractVector, AbstractRange}`: time vector
  * `a::Real`: amplitude

**Returns**

  * cs::Vector{ComplexF64}`

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

  * `g::Vector{Float64}`

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

  * `mw::Vector{ComplexF64}`

**Source**

Cohen MX. A better way to define and describe Morlet wavelets for time-frequency analysis. NeuroImage. 2019 Oct;199:81–6. 

<a id='NeuroAnalyzer.generate_noise' href='#NeuroAnalyzer.generate_noise'>#</a>
**`NeuroAnalyzer.generate_noise`** &mdash; *Function*.



```julia
generate_noise(n, amp; type)
```

Generate noise.

**Arguments**

  * `n::Int64`: length (in samples)
  * `amp::Real=1.0`: amplitude, signal amplitude will be in `[-amp, +amp]`
  * `type::Symbol=:whiten`: noise type:

      * `:whiten`: normal distributed
      * `:whiteu`: uniformly distributed
      * `:pink`

**Returns**

  * `s::Float64`

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

  * `s::Vector{Float64}`

<a id='NeuroAnalyzer.generate_sine' href='#NeuroAnalyzer.generate_sine'>#</a>
**`NeuroAnalyzer.generate_sine`** &mdash; *Function*.



```julia
generate_sine(f, t, a, p)
```

Generates sine wave.

**Arguments**

  * `f::Real`: frequency [Hz]
  * `t::Union{AbstractVector, AbstractRange}`: time vector
  * `a::Real`: amplitude
  * `p::Real`: initial phase

**Returns**

  * s::Vector{Float64}`

<a id='NeuroAnalyzer.generate_square' href='#NeuroAnalyzer.generate_square'>#</a>
**`NeuroAnalyzer.generate_square`** &mdash; *Function*.



```julia
generate_square(t, a, p, w, offset)
```

Generates square wave.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: time vector
  * `a::Real`: amplitude
  * `p::Real`: duty cycle
  * `w::Real`: width
  * `offset::Real`: amplitude offset

**Returns**

  * `s::Vector{Float64}`

<a id='NeuroAnalyzer.generate_triangle' href='#NeuroAnalyzer.generate_triangle'>#</a>
**`NeuroAnalyzer.generate_triangle`** &mdash; *Function*.



```julia
generate_triangle(t, a)
```

Generates triangle wave.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: time vector
  * `a::Real`: amplitude

**Returns**

  * `s::Vector{Float64}`

<a id='NeuroAnalyzer.generate_window' href='#NeuroAnalyzer.generate_window'>#</a>
**`NeuroAnalyzer.generate_window`** &mdash; *Function*.



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
  * `even::Bool=false`: if true, make the window of even length (increase length by 1 for odd value of `n`)

**Returns**

  * `w::Vector{Float64}`:: generated window

<a id='NeuroAnalyzer.gradient' href='#NeuroAnalyzer.gradient'>#</a>
**`NeuroAnalyzer.gradient`** &mdash; *Function*.



gradient(x; rev)

Calculate gradient of a 1-dimensional scalar field.

**Arguments**

  * `x::AbstractVector`
  * `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

**Returns**

Named tuple containing:

  * `g::Vector{Vector{Float64}}`: vector field of gradients
  * `g_len::Vector{Float64}`: scalar field of gradient lengths


gradient(x; rev)

Calculate gradient of a 2-dimensional scalar field.

**Arguments**

  * `x::AbstractMatrix`
  * `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

**Returns**

Named tuple containing:

  * `g::Matrix{Vector{Float64}}`: vector field of gradients
  * `g_len::Matrix{Float64}`: scalar field of gradient lengths


gradient(x; rev)

Calculate gradient of a ≥3-dimensional scalar field.

**Arguments**

  * `x::AbstractArray`
  * `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

**Returns**

Named tuple containing:

  * `g::Array{Vector{Float64}, 3}`: vector field of gradients
  * `g_len::Array{Float64, 3}`: scalar field of gradient lengths

<a id='NeuroAnalyzer.history' href='#NeuroAnalyzer.history'>#</a>
**`NeuroAnalyzer.history`** &mdash; *Function*.



```julia
history(obj)
```

Show processing history.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `history::Vector{String}`

<a id='NeuroAnalyzer.hz2rads' href='#NeuroAnalyzer.hz2rads'>#</a>
**`NeuroAnalyzer.hz2rads`** &mdash; *Function*.



```julia
hz2rads(f)
```

Convert frequency in Hz to rad/s.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroAnalyzer.ifft0' href='#NeuroAnalyzer.ifft0'>#</a>
**`NeuroAnalyzer.ifft0`** &mdash; *Function*.



```julia
ifft0(x, n)
```

Perform IFFT of zero-padded vector.

**Arguments**

  * `x::AbstractVector`
  * `n::Int64`: number of zeros added to `x`

**Returns**

  * `ifft0::Vector{ComplexF64}`: reconstructed signal trimmed to original length

<a id='NeuroAnalyzer.info' href='#NeuroAnalyzer.info'>#</a>
**`NeuroAnalyzer.info`** &mdash; *Function*.



```julia
info(obj)
```

Show info.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.l1' href='#NeuroAnalyzer.l1'>#</a>
**`NeuroAnalyzer.l1`** &mdash; *Function*.



```julia
l1(a1, a2)
```

Compare two arrays (e.g. two spectrograms), using L1 (Manhattan) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l1::Float64`

<a id='NeuroAnalyzer.l2' href='#NeuroAnalyzer.l2'>#</a>
**`NeuroAnalyzer.l2`** &mdash; *Function*.



```julia
l2(a1, a2)
```

Compare two arrays (e.g. two spectrograms), using L2 (Euclidean) distance.

**Arguments**

  * `a1::AbstractArray`: first array
  * `a2::AbstractArray`: second array

**Returns**

  * `l2::Float64`

<a id='NeuroAnalyzer.labels' href='#NeuroAnalyzer.labels'>#</a>
**`NeuroAnalyzer.labels`** &mdash; *Function*.



```julia
labels(obj)
```

Return channel labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroAnalyzer.linspace' href='#NeuroAnalyzer.linspace'>#</a>
**`NeuroAnalyzer.linspace`** &mdash; *Function*.



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

<a id='NeuroAnalyzer.list_component' href='#NeuroAnalyzer.list_component'>#</a>
**`NeuroAnalyzer.list_component`** &mdash; *Function*.



```julia
list_component(obj)
```

List component names.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `components::Vector{Symbol}`

<a id='NeuroAnalyzer.logspace' href='#NeuroAnalyzer.logspace'>#</a>
**`NeuroAnalyzer.logspace`** &mdash; *Function*.



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

<a id='NeuroAnalyzer.make_table' href='#NeuroAnalyzer.make_table'>#</a>
**`NeuroAnalyzer.make_table`** &mdash; *Function*.



```julia
make_table(; header, data)
```

Display data as a table.

**Arguments**

  * `header::Matrix{String}`: table header, e.g. `header = ["Group" "A" "B"]`
  * `data::Matrix{Any}`: table data, e.g. `data = ["var1" 1.0 2.0; "var2" 3.0 4.0]`

<a id='NeuroAnalyzer.maxat' href='#NeuroAnalyzer.maxat'>#</a>
**`NeuroAnalyzer.maxat`** &mdash; *Function*.



```julia
maxat(x, y)
```

Find maximum value of one vector and return value at its index from another vector.

**Argument**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `value::eltype(y)`
  * `idx::Int64`

<a id='NeuroAnalyzer.minat' href='#NeuroAnalyzer.minat'>#</a>
**`NeuroAnalyzer.minat`** &mdash; *Function*.



```julia
minat(x, y)
```

Find minimum value of one vector and return value at its index from another vector.

**Argument**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `value::eltype(y)`
  * `idx::Int64`

<a id='NeuroAnalyzer.m_norm' href='#NeuroAnalyzer.m_norm'>#</a>
**`NeuroAnalyzer.m_norm`** &mdash; *Function*.



```julia
m_norm(m)
```

Normalize matrix.

**Arguments**

  * `m::AbstractArray`

**Returns**

  * `m_norm::AbstractArray`

<a id='NeuroAnalyzer.m_pad0' href='#NeuroAnalyzer.m_pad0'>#</a>
**`NeuroAnalyzer.m_pad0`** &mdash; *Function*.



```julia
m_pad0(m)
```

Pad matrix with zeros to make it square.

**Arguments**

  * `m::Matrix{<:Number}`

**Returns**

  * `m::Matrix{<:Number}`

<a id='NeuroAnalyzer.m_sort' href='#NeuroAnalyzer.m_sort'>#</a>
**`NeuroAnalyzer.m_sort`** &mdash; *Function*.



```julia
m_sort(m, m_idx; rev, dims)
```

Sorts matrix using sorting index.

**Arguments**

  * `m::Matrix`
  * `m_idx::Vector{Int64}`: sorting index
  * `rev::Bool=false`: reverse sort
  * `dims::Int64=1`: sort by columns (`dims=1`) or by rows (`dims=2`)

**Returns**

  * `m_sorted::Matrix`

<a id='NeuroAnalyzer.m_sortperm' href='#NeuroAnalyzer.m_sortperm'>#</a>
**`NeuroAnalyzer.m_sortperm`** &mdash; *Function*.



```julia
m_sortperm(m; rev, dims)
```

Generates matrix sorting index.

**Arguments**

  * `m::AbstractMatrix`
  * `rev::Bool`: reverse sort
  * `dims::Int64=1`: sort by columns (`dims=1`) or by rows (`dims=2`)

**Returns**

  * `idx::Matrix{Int64}`

<a id='NeuroAnalyzer.nchannels' href='#NeuroAnalyzer.nchannels'>#</a>
**`NeuroAnalyzer.nchannels`** &mdash; *Function*.



```julia
nchannels(obj; type)
```

Return number of channels of `type`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::String="all"`: channel type (stored in the global channel_types constant variable)

**Returns**

  * `ch_n::Int64`


```
nchannels(study)
```

Return number of channels per NeuroAnalyzer NEURO object in the study.

**Arguments**

  * `study::NeuroAnalyzer.STUDY`

**Returns**

  * `n::Int64`

<a id='NeuroAnalyzer.nepochs' href='#NeuroAnalyzer.nepochs'>#</a>
**`NeuroAnalyzer.nepochs`** &mdash; *Function*.



```julia
nepochs(obj)
```

Return number of epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `ep_n::Int64`


```
nepochs(study)
```

Return number of epochs per NeuroAnalyzer NEURO object in the study.

**Arguments**

  * `study::NeuroAnalyzer.STUDY`

**Returns**

  * `n::Int64`

<a id='NeuroAnalyzer.nextpow2' href='#NeuroAnalyzer.nextpow2'>#</a>
**`NeuroAnalyzer.nextpow2`** &mdash; *Function*.



```julia
nextpow2(x)
```

Return the next power of 2 for a given number.

**Argument**

  * `x::Int64`

**Returns**

  * `nextpow2::Int64`

<a id='NeuroAnalyzer.optode_labels' href='#NeuroAnalyzer.optode_labels'>#</a>
**`NeuroAnalyzer.optode_labels`** &mdash; *Function*.



```julia
optode_labels(obj)
```

Return optode labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroAnalyzer.pad0' href='#NeuroAnalyzer.pad0'>#</a>
**`NeuroAnalyzer.pad0`** &mdash; *Function*.



```julia
pad0(x, n)
```

Pad row(s) with zeros. Works with 1-, 2- and 3-dimensional arrays.

**Arguments**

  * `x::Union{AbstractVector, AbstractArray}`
  * `n::Int64`: padding length (number of zeros to add)

**Returns**

  * `pad0::Union{AbstractVector, AbstractArray`

<a id='NeuroAnalyzer.pad2' href='#NeuroAnalyzer.pad2'>#</a>
**`NeuroAnalyzer.pad2`** &mdash; *Function*.



```julia
pad2(x)
```

Pad row(s) with zeros to the nearest power of 2 length. Works with 1-, 2- and 3-dimensional arrays.

**Arguments**

  * `x::Union{AbstractVector, AbstractArray}`

**Returns**

  * `pad2::Union{AbstractVector, AbstractArray`

<a id='NeuroAnalyzer.padm' href='#NeuroAnalyzer.padm'>#</a>
**`NeuroAnalyzer.padm`** &mdash; *Function*.



```julia
padm(x, n)
```

Pad row(s) with mean value(s). Works with 1-, 2- and 3-dimensional arrays.

**Arguments**

  * `x::Union{AbstractVector, AbstractArray}`
  * `n::Int64`: padding length (number of values to add)
  * `mode::Symbol=:row`: how the mean is calculated:

      * `:all`: mean of all rows
      * `:row`: separate mean per each row

**Returns**

  * `padm::Union{AbstractVector, AbstractArray`

<a id='NeuroAnalyzer.paired_labels' href='#NeuroAnalyzer.paired_labels'>#</a>
**`NeuroAnalyzer.paired_labels`** &mdash; *Function*.



```julia
paired_labels(l; unq)
```

Return paired labels.

**Arguments**

  * `l::Vector{String}`
  * `unq::Bool=true`: if true, do not add pairs of the same labels, e.g. "label1-label1"

**Returns**

  * `l_paired::Vector{String}`: paired labels


```
paired_labels(l1, l2)
```

Return paired labels.

**Arguments**

  * `l1::Vector{String}`
  * `l2::Vector{String}`

**Returns**

  * `l_paired::Vector{String}`: paired labels

<a id='NeuroAnalyzer.perm_cmp' href='#NeuroAnalyzer.perm_cmp'>#</a>
**`NeuroAnalyzer.perm_cmp`** &mdash; *Function*.



```julia
perm_cmp(a1, a2; p, perm_n)
```

Compare two 3-dimensional arrays (e.g. two spectrograms), using permutation based statistic.

**Arguments**

  * `a1::Array{<:Real, 3}`: first array
  * `a2::Array{<:Real, 3}`: second array
  * `p::Float64=0.05`: p-value
  * `perm_n::Int64=1000`: number of permutations

**Returns**

Named tuple containing:

  * `zmap::Array{Float64, 3}`: array of Z-values
  * `bm::Array{Float64, 3}`: binarized mask of statistically significant positions

<a id='NeuroAnalyzer.phases' href='#NeuroAnalyzer.phases'>#</a>
**`NeuroAnalyzer.phases`** &mdash; *Function*.



```julia
phases(s)
```

Calculate phases.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `phases::Vector{Float64}`

<a id='NeuroAnalyzer.rads2hz' href='#NeuroAnalyzer.rads2hz'>#</a>
**`NeuroAnalyzer.rads2hz`** &mdash; *Function*.



```julia
rads2hz(f)
```

Convert frequency in rad/s to Hz.

**Arguments**

  * `f::Real`

**Returns**

  * `f_rads::Float64`

<a id='NeuroAnalyzer.rename_component' href='#NeuroAnalyzer.rename_component'>#</a>
**`NeuroAnalyzer.rename_component`** &mdash; *Function*.



```julia
rename_component(obj, c_old, c_new)
```

Rename component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

**Return**

  * `obj_new::NEURO`

<a id='NeuroAnalyzer.rename_component!' href='#NeuroAnalyzer.rename_component!'>#</a>
**`NeuroAnalyzer.rename_component!`** &mdash; *Function*.



```julia
rename_component!(obj, c_old, c_new)
```

Rename component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c_old::Symbol`: old component name
  * `c_new::Symbol`: new component name

<a id='NeuroAnalyzer.reset_components' href='#NeuroAnalyzer.reset_components'>#</a>
**`NeuroAnalyzer.reset_components`** &mdash; *Function*.



```julia
reset_components(obj)
```

Remove all components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reset_components!' href='#NeuroAnalyzer.reset_components!'>#</a>
**`NeuroAnalyzer.reset_components!`** &mdash; *Function*.



```julia
reset_components!(obj)
```

Remove all components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.rfft0' href='#NeuroAnalyzer.rfft0'>#</a>
**`NeuroAnalyzer.rfft0`** &mdash; *Function*.



```julia
rfft0(x, n)
```

Perform zeros-padded single-sided FFT.

**Arguments**

  * `x::AbstractVector`
  * `n::Int64`: number of zeros to add

**Returns**

  * `rfft0::Vector{ComplexF64}`

<a id='NeuroAnalyzer.rfft2' href='#NeuroAnalyzer.rfft2'>#</a>
**`NeuroAnalyzer.rfft2`** &mdash; *Function*.



```julia
rfft2(x)
```

Perform zeros-padded single-sided FFT, so the length of padded vector is a power of 2.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `rfft2::Vector{ComplexF64}`

<a id='NeuroAnalyzer.s2t' href='#NeuroAnalyzer.s2t'>#</a>
**`NeuroAnalyzer.s2t`** &mdash; *Function*.



```julia
s2t(s, fs)
```

Convert sample number to time in seconds.

**Arguments**

  * `t::Int64`: sample number
  * `fs::Int64`: sampling rate

**Returns**

  * `s2t::Float64`: time in s


```
s2t(obj; s)
```

Convert time in samples to seconds.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `s::Int64`: time in samples

**Returns**

  * `s2t::Float64`: time in seconds

<a id='NeuroAnalyzer.signal_len' href='#NeuroAnalyzer.signal_len'>#</a>
**`NeuroAnalyzer.signal_len`** &mdash; *Function*.



```julia
signal_len(obj)
```

Return signal length.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `s_len::Int64`

<a id='NeuroAnalyzer.source_labels' href='#NeuroAnalyzer.source_labels'>#</a>
**`NeuroAnalyzer.source_labels`** &mdash; *Function*.



```julia
source_labels(obj)
```

Return NIRS source labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `labels::Vector{String}`

<a id='NeuroAnalyzer.sr' href='#NeuroAnalyzer.sr'>#</a>
**`NeuroAnalyzer.sr`** &mdash; *Function*.



```julia
sr(obj)
```

Return sampling rate.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `sr::Int64`


```
sr(study)
```

Return sampling rate of NeuroAnalyzer NEURO objects in the study.

**Arguments**

  * `study::NeuroAnalyzer.STUDY`

**Returns**

  * `sr::Int64`

<a id='NeuroAnalyzer.t2f' href='#NeuroAnalyzer.t2f'>#</a>
**`NeuroAnalyzer.t2f`** &mdash; *Function*.



```julia
t2f(t)
```

Convert cycle length in ms to frequency.

**Arguments**

  * `t::Real`: cycle length in ms

**Returns**

  * `f::Float64`: frequency in Hz

<a id='NeuroAnalyzer.t2s' href='#NeuroAnalyzer.t2s'>#</a>
**`NeuroAnalyzer.t2s`** &mdash; *Function*.



```julia
t2s(t, fs)
```

Convert time in seconds to sample number.

**Arguments**

  * `t::T`: time in s
  * `fs::Int64`: sampling rate

**Returns**

  * `t2s::Int64`: sample number


```
t2s(obj; t)
```

Convert time in seconds to sample number.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `t::T`: time in seconds

**Returns**

  * `t2s::Int64`: time in samples

<a id='NeuroAnalyzer.tavg' href='#NeuroAnalyzer.tavg'>#</a>
**`NeuroAnalyzer.tavg`** &mdash; *Function*.



```julia
tavg(s)
```

Average signal across trials.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::AbstractArray`

<a id='NeuroAnalyzer.to_df' href='#NeuroAnalyzer.to_df'>#</a>
**`NeuroAnalyzer.to_df`** &mdash; *Function*.



```julia
to_df(obj)
```

Export object data as DataFrame.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `df::DataFrame`: DataFrame containing time points and channels

<a id='NeuroAnalyzer.tuple_order' href='#NeuroAnalyzer.tuple_order'>#</a>
**`NeuroAnalyzer.tuple_order`** &mdash; *Function*.



```julia
tuple_order(t, rev)
```

Order tuple elements in ascending or descending (`rev=true`) order.

**Arguments**

  * `t::Tuple{Real, Real}`
  * `rev::Bool=false`

**Returns**

  * `t::Tuple{Real, Real}`

<a id='NeuroAnalyzer.vec2mat' href='#NeuroAnalyzer.vec2mat'>#</a>
**`NeuroAnalyzer.vec2mat`** &mdash; *Function*.



```julia
vec2mat(x; wlen, woverlap)
```

Reshape vector into matrix using fixed segment length and overlapping.

**Arguments**

  * `x::AbstractVector`
  * `wlen::Int64`: window length (in samples)
  * `woverlap::Int64`: overlap with the previous window (in samples)

**Returns**

  * `m::Matrix{eltype(x)}`

<a id='NeuroAnalyzer.view_header' href='#NeuroAnalyzer.view_header'>#</a>
**`NeuroAnalyzer.view_header`** &mdash; *Function*.



```julia
header(obj)
```

Show keys and values of the object header.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.view_note' href='#NeuroAnalyzer.view_note'>#</a>
**`NeuroAnalyzer.view_note`** &mdash; *Function*.



```julia
view_note(obj)
```

Return the object recording note.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.vreduce' href='#NeuroAnalyzer.vreduce'>#</a>
**`NeuroAnalyzer.vreduce`** &mdash; *Function*.



```julia
vreduce(x, f; n)
```

Reduce two vectors at indices of the second vector being multiplications of a constant. Useful e.g. for simplifying values across frequencies, when the number of frequencies (and thus values) is high.

**Arguments**

  * `x::AbstractVector`: e.g. signal data
  * `f::AbstractVector`: e.g. frequencies
  * `n::Float64=0.5`: reduce at multiplications of this value

**Returns**

  * `x_new::Vector{eltype(x)}`
  * `f_new::Vector{eltype(f)}`

<a id='NeuroAnalyzer.vsearch' href='#NeuroAnalyzer.vsearch'>#</a>
**`NeuroAnalyzer.vsearch`** &mdash; *Function*.



```julia
vsearch(y, x; acc)
```

Return the positions of the value in the vector.

**Arguments**

  * `y::T`: value of interest
  * `x::AbstractVector`: vector to search within
  * `acc::Bool=false`: if true, return the difference between `y` and `x[idx]`

**Returns**

  * `idx::Int64`
  * `d::Real`: the difference between `y` and `x[idx]`


```
vsearch(y, x; acc)
```

Return the positions of the value in the vector.

**Arguments**

  * `y::AbstractVector`: vector of interest
  * `x::AbstractVector`: vector to search within
  * `acc::Bool=false`: if true, return the difference between `y` and `x[idx:idx + length(y)]`

**Returns**

  * `idx::Int64`
  * `d::Real`: the difference between `y` and `x[idx:idx + length(y)]`

<a id='NeuroAnalyzer.vsplit' href='#NeuroAnalyzer.vsplit'>#</a>
**`NeuroAnalyzer.vsplit`** &mdash; *Function*.



```julia
vsplit(x, n)
```

Splits vector into pieces.

**Argument**

  * `x::AbstractVector`
  * `n::Int64`: length of one piece

**Returns**

  * `x::Vector{AbstractVector}`


<a id='IO'></a>

<a id='IO-1'></a>

## IO

<a id='NeuroAnalyzer.export_csv' href='#NeuroAnalyzer.export_csv'>#</a>
**`NeuroAnalyzer.export_csv`** &mdash; *Function*.



```julia
export_csv(obj; file_name, header, components, markers, overwrite)
```

Export `NeuroAnalyzer.NEURO` object to CSV.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`
  * `header::Bool=false`: export header
  * `epoch_time::Bool=false`: export epoch time points
  * `components::Bool=false`: export components
  * `markers::Bool=false`: export event markers
  * `locs::Bool=false`: export channel locations
  * `history::Bool=false`: export history
  * `overwrite::Bool=false`

<a id='NeuroAnalyzer.export_locs' href='#NeuroAnalyzer.export_locs'>#</a>
**`NeuroAnalyzer.export_locs`** &mdash; *Function*.



```julia
export_locs(obj; file_name, overwrite)
```

Export channel locations data, format is based on `file_name` extension (.ced, .locs or .tsv)

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`
  * `overwrite::Bool=false`


```
export_locs(locs; file_name, overwrite)
```

Export channel locations, format is based on `file_name` extension (.ced, .locs, .tsv)

**Arguments**

  * `locs::DataFrame`
  * `file_name::String`
  * `overwrite::Bool=false`

**Returns**

  * `success::Bool`

<a id='NeuroAnalyzer.export_markers' href='#NeuroAnalyzer.export_markers'>#</a>
**`NeuroAnalyzer.export_markers`** &mdash; *Function*.



```julia
export_markers(obj; file_name, overwrite)
```

Export `NeuroAnalyzer.NEURO` object markers to CSV.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`
  * `overwrite::Bool=false`

<a id='NeuroAnalyzer.import_alice4' href='#NeuroAnalyzer.import_alice4'>#</a>
**`NeuroAnalyzer.import_alice4`** &mdash; *Function*.



```julia
import_alice4(file_name; detect_type)
```

Load EDF exported from Alice 4 Polysomnography System and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

  * EDF files exported from Alice 4 have incorrect value of `data_records` (-1) and multiple sampling rate; channels are upsampled to the highest rate.

<a id='NeuroAnalyzer.import_bdf' href='#NeuroAnalyzer.import_bdf'>#</a>
**`NeuroAnalyzer.import_bdf`** &mdash; *Function*.



```julia
import_bdf(file_name; default_type)
```

Load BDF/BDF+ file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on channel label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

  * sampling_rate = n.samples ÷ data.record.duration
  * gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
  * value = (value - digital minimum ) × gain + physical minimum

**Source**

https://www.biosemi.com/faq/file_format.htm

<a id='NeuroAnalyzer.import_bv' href='#NeuroAnalyzer.import_bv'>#</a>
**`NeuroAnalyzer.import_bv`** &mdash; *Function*.



```julia
import_bv(file_name; detect_type)
```

Load BrainVision BVCDF file and return `NeuroAnalyzer.NEURO` object. At least two files are required: .vhdr (header) and .eeg (signal data). If available, markers are loaded from .vmrk file.

**Arguments**

  * `file_name::String`: name of the file to load, should point to .vhdr file.
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_cnt' href='#NeuroAnalyzer.import_cnt'>#</a>
**`NeuroAnalyzer.import_cnt`** &mdash; *Function*.



```julia
import_cnt(file_name; data_format, detect_type)
```

Load Neuroscan continuous signal file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label
  * `data_format::Symbol=:i32`: Neuroscan stores data either in 16-bit (`:i16`) or 32-bit (`:i32`) representation, but does not say so in the file format

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

Based on loadcnt.m by Sean Fitzgibbon and Arnaud Delorme (https://cnl.salk.edu/~arno/cntload/index.html)

<a id='NeuroAnalyzer.import_csv' href='#NeuroAnalyzer.import_csv'>#</a>
**`NeuroAnalyzer.import_csv`** &mdash; *Function*.



```julia
import_csv(file_name; detect_type)
```

Load CSV file (e.g. exported from EEGLAB) and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

CSV first row or first column must contain channel names. Shape of data array will be detected automatically. Sampling rate will be detected. If file is gzip-ed, it will be uncompressed automatically while reading.

<a id='NeuroAnalyzer.import_dat' href='#NeuroAnalyzer.import_dat'>#</a>
**`NeuroAnalyzer.import_dat`** &mdash; *Function*.



```julia
import_dat(file_name)
```

Load Neuroscan DAT file.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `dat::DataFrame`

<a id='NeuroAnalyzer.import_digitrack' href='#NeuroAnalyzer.import_digitrack'>#</a>
**`NeuroAnalyzer.import_digitrack`** &mdash; *Function*.



```julia
import_digitrack(file_name; detect_type)
```

Load Digitrack ASCII file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_duomag' href='#NeuroAnalyzer.import_duomag'>#</a>
**`NeuroAnalyzer.import_duomag`** &mdash; *Function*.



```julia
import_duomag(file_name)
```

Load DuoMAG TMS MEP recording file (.ascii or .m) and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_edf' href='#NeuroAnalyzer.import_edf'>#</a>
**`NeuroAnalyzer.import_edf`** &mdash; *Function*.



```julia
import_edf(file_name; detect_type)
```

Load EDF/EDF+ file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

  * sampling_rate = n.samples ÷ data.record.duration
  * gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
  * value = (value - digital minimum ) × gain + physical minimum

**Source**

1. Kemp B, Varri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992; 82(5): 391–3
2. Kemp B, Olivan J. European data format ‘plus’(EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003; 114: 1755–61
3. https://www.edfplus.info/specs/

<a id='NeuroAnalyzer.import_edf_annotations' href='#NeuroAnalyzer.import_edf_annotations'>#</a>
**`NeuroAnalyzer.import_edf_annotations`** &mdash; *Function*.



```julia
import_edf_annotations(file_name; detect_type)
```

Load annotations from EDF+ file and return `markers` DataFrame. This function is used for EDF+ files containing annotations only.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `markers::DataFrame`

<a id='NeuroAnalyzer.import_fiff' href='#NeuroAnalyzer.import_fiff'>#</a>
**`NeuroAnalyzer.import_fiff`** &mdash; *Function*.



```julia
import_fiff(file_name; detect_type)
```

Load FIFF (Functional Image File Format) file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

Elekta Neuromag: Functional Image File Format Description. FIFF version 1.3. March 2011

<a id='NeuroAnalyzer.import_gdf' href='#NeuroAnalyzer.import_gdf'>#</a>
**`NeuroAnalyzer.import_gdf`** &mdash; *Function*.



```julia
import_gdf(file_name; detect_type)
```

Load GDF file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Notes**

  * sampling_rate = n.samples ÷ data.record.duration
  * gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
  * value = (value - digital minimum ) × gain + physical minimum

**Source**

1. Schlögl A, Filz O, Ramoser H, Pfurtscheller G. GDF - A General Dataformat for Biosignals Version 1.25. 1998
2. Schlögl, A. GDF - A General Dataformat for Biosignals Version 2.12. 2009
3. Schlögl, A. GDF - A General Dataformat for Biosignals Version 2.51. 2013

<a id='NeuroAnalyzer.import_locs' href='#NeuroAnalyzer.import_locs'>#</a>
**`NeuroAnalyzer.import_locs`** &mdash; *Function*.



```julia
import_locs(file_name)
```

Load channel locations. Supported formats:

  * CED
  * ELC
  * LOCS
  * TSV
  * SFP
  * CSD
  * GEO
  * MAT
  * TXT
  * DAT
  * ASC

This is a meta-function that triggers appropriate `import_locs_*()` function. File format is detected based on file extension.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_asc' href='#NeuroAnalyzer.import_locs_asc'>#</a>
**`NeuroAnalyzer.import_locs_asc`** &mdash; *Function*.



```julia
import_locs_asc(file_name)
```

Load channel locations from ASC file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_ced' href='#NeuroAnalyzer.import_locs_ced'>#</a>
**`NeuroAnalyzer.import_locs_ced`** &mdash; *Function*.



```julia
import_locs_ced(file_name)
```

Load channel locations from CED file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_csd' href='#NeuroAnalyzer.import_locs_csd'>#</a>
**`NeuroAnalyzer.import_locs_csd`** &mdash; *Function*.



```julia
import_locs_csd(file_name)
```

Load channel locations from CSD file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_dat' href='#NeuroAnalyzer.import_locs_dat'>#</a>
**`NeuroAnalyzer.import_locs_dat`** &mdash; *Function*.



```julia
import_locs_dat(file_name)
```

Load channel locations from DAT file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_elc' href='#NeuroAnalyzer.import_locs_elc'>#</a>
**`NeuroAnalyzer.import_locs_elc`** &mdash; *Function*.



```julia
import_locs_elc(file_name)
```

Load channel locations from ELC file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_geo' href='#NeuroAnalyzer.import_locs_geo'>#</a>
**`NeuroAnalyzer.import_locs_geo`** &mdash; *Function*.



```julia
import_locs_geo(file_name)
```

Load channel locations from GEO file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_locs' href='#NeuroAnalyzer.import_locs_locs'>#</a>
**`NeuroAnalyzer.import_locs_locs`** &mdash; *Function*.



```julia
import_locs_locs(file_name)
```

Load channel locations from LOCS file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_mat' href='#NeuroAnalyzer.import_locs_mat'>#</a>
**`NeuroAnalyzer.import_locs_mat`** &mdash; *Function*.



```julia
import_locs_mat(file_name)
```

Load channel locations from MAT file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_sfp' href='#NeuroAnalyzer.import_locs_sfp'>#</a>
**`NeuroAnalyzer.import_locs_sfp`** &mdash; *Function*.



```julia
import_locs_sfp(file_name)
```

Load channel locations from SFP file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_tsv' href='#NeuroAnalyzer.import_locs_tsv'>#</a>
**`NeuroAnalyzer.import_locs_tsv`** &mdash; *Function*.



```julia
import_locs_tsv(file_name)
```

Load channel locations from TSV file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_locs_txt' href='#NeuroAnalyzer.import_locs_txt'>#</a>
**`NeuroAnalyzer.import_locs_txt`** &mdash; *Function*.



```julia
import_locs_txt(file_name)
```

Load channel locations from TXT file.

**Arguments**

  * `file_name::String`

**Returns**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.import_montage' href='#NeuroAnalyzer.import_montage'>#</a>
**`NeuroAnalyzer.import_montage`** &mdash; *Function*.



```julia
import_montage(file_name)
```

Load montage from a text file. Example montage files are located in the `montages/` folder. The structure of the file is:

  * first line: name of the montage, e.g. `longitudinal-BIP`
  * next lines: channel pairs or individual channels, e.g. `Fz-Cz` or `Fp1`

Each channel/channel pair must be in a separate line

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

Named tuple containing:

  * `ref_list::Vector{String}`: list of channel pairs
  * `ref_name::Vector{String}`: name of the montage

<a id='NeuroAnalyzer.import_ncs' href='#NeuroAnalyzer.import_ncs'>#</a>
**`NeuroAnalyzer.import_ncs`** &mdash; *Function*.



```julia
import_ncs(file_name)
```

Load Neuralinx Continuously Sampled Channels (CSC) and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_nirs' href='#NeuroAnalyzer.import_nirs'>#</a>
**`NeuroAnalyzer.import_nirs`** &mdash; *Function*.



```julia
import_nirs(file_name)
```

Load NIRS file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

https://github.com/BUNPC/Homer3/wiki/HOMER3-file-formats

<a id='NeuroAnalyzer.import_nirx' href='#NeuroAnalyzer.import_nirx'>#</a>
**`NeuroAnalyzer.import_nirx`** &mdash; *Function*.



```julia
import_nirx(file_name)
```

Load NIRX file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load, should point to .hdr file.

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

https://nirx.net/file-formats

<a id='NeuroAnalyzer.import_npy' href='#NeuroAnalyzer.import_npy'>#</a>
**`NeuroAnalyzer.import_npy`** &mdash; *Function*.



```julia
import_npy(file_name; sampling_rate)
```

Load NPY file (exported from MNE) and return `NeuroAnalyzer.NEURO` object. Data type and channel types are set as  is EEG.

**Arguments**

  * `file_name::String`: name of the file to load
  * `sampling_rate::Int64`: NPY file contains only signal data, therefore its sampling rate must be provided upon importing

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_nwb' href='#NeuroAnalyzer.import_nwb'>#</a>
**`NeuroAnalyzer.import_nwb`** &mdash; *Function*.



```julia
import_nwb(file_name; detect_type)
```

Load EEG data from Neurodata Without Borders (NWB) file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

1. https://www.biorxiv.org/content/10.1101/523035v1

<a id='NeuroAnalyzer.import_recording' href='#NeuroAnalyzer.import_recording'>#</a>
**`NeuroAnalyzer.import_recording`** &mdash; *Function*.



```julia
import_recording(file_name; detect_type)
```

Load recording file and return `NeuroAnalyzer.NEURO` object. Supported formats:

  * EDF/EDF+
  * BDF/BDF+
  * GDF
  * BrainVision
  * CSV
  * SET (EEGLAB dataset)
  * NWB (Neurodata Without Borders)
  * FIFF
  * SNIRF
  * NIRS

This is a meta-function that triggers appropriate `import_*()` function. File format is detected based on file extension (.edf|.bdf|.gdf|.vhdr|.ahdr|.csv|.csv.gz|.set|.nwb|.fif|.fiff|.snirf|.nirs). Signal data type (e.g. EEG or MEG is auto-detected) and stored in the `obj.header.recording[:data_type]` field.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on channel label
  * `n::Int64=0`: subject number to extract in case of multi-subject file

**Returns**

  * `::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.import_set' href='#NeuroAnalyzer.import_set'>#</a>
**`NeuroAnalyzer.import_set`** &mdash; *Function*.



```julia
import_set(file_name; detect_type)
```

Load SET file (exported from EEGLAB) and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `detect_type::Bool=true`: detect channel type based on its label

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

1. https://eeglab.org/tutorials/ConceptsGuide/Data_Structures.html

<a id='NeuroAnalyzer.import_snirf' href='#NeuroAnalyzer.import_snirf'>#</a>
**`NeuroAnalyzer.import_snirf`** &mdash; *Function*.



```julia
import_snirf(file_name; n)
```

Load Shared Near Infrared Spectroscopy Format (SNIRF) file and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load
  * `n::Int64=0`: subject number to extract in case of multi-subject file

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

https://github.com/fNIRS/snirf/blob/v1.1/snirf_specification.md

<a id='NeuroAnalyzer.import_thymatron' href='#NeuroAnalyzer.import_thymatron'>#</a>
**`NeuroAnalyzer.import_thymatron`** &mdash; *Function*.



```julia
import_thymatron(file_name)
```

Import scanned images of EEG printed by Thymatron ECT equipment and return `NeuroAnalyzer.NEURO` object.

Usually there are three images: baseline, during seizure activity and after seizure activity. If more then one image is provided, images are added as consecutive channels in the same object.

Image properties:

  * 100 DPI
  * 490 x 100 px
  * height 2.5 cm
  * width 12.5 cm

**Arguments**

  * `file_name::String`: name of the file to load
  * `dpi::Int64=100`: DPI of the scanned images

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

**Source**

Wysokiński A. EEG_ADC: Digitizer and Analyzer of Electroconvulsive Therapy Paper Electroencephalogram Recordings. JECT 2022; 4: 255-256

<a id='NeuroAnalyzer.import_xdf' href='#NeuroAnalyzer.import_xdf'>#</a>
**`NeuroAnalyzer.import_xdf`** &mdash; *Function*.



```julia
import_xdf(file_name)
```

Load Extensible Data Format (XDF) and return `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.load' href='#NeuroAnalyzer.load'>#</a>
**`NeuroAnalyzer.load`** &mdash; *Function*.



```julia
load(file_name)
```

Load `NeuroAnalyzer.NEURO` object from `file_name` file (HDF5-based).

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.load_locs' href='#NeuroAnalyzer.load_locs'>#</a>
**`NeuroAnalyzer.load_locs`** &mdash; *Function*.



```julia
load_locs(obj; file_name)
```

Load channel locations from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame. 

Accepted formats:

  * CED
  * LOCS
  * ELC
  * TSV
  * SFP
  * CSD
  * GEO
  * MAT
  * TXT
  * DAT
  * ASC

Channel locations:

  * `loc_theta`: polar angle
  * `loc_radius`: polar radius
  * `loc_x`: spherical Cartesian x
  * `loc_y`: spherical Cartesian y
  * `loc_z`: spherical Cartesian z
  * `loc_radius_sph`: spherical radius
  * `loc_theta_sph`: spherical horizontal angle
  * `loc_phi_sph`: spherical azimuth angle

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.load_locs!' href='#NeuroAnalyzer.load_locs!'>#</a>
**`NeuroAnalyzer.load_locs!`** &mdash; *Function*.



```julia
load_locs!(obj; file_name, normalize)
```

Load channel locations from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame. 

Accepted formats:

  * CED
  * LOCS
  * ELC
  * TSV
  * SFP
  * CSD
  * GEO
  * MAT
  * TXT
  * DAT
  * ASC

Channel locations:

  * `loc_theta`: polar angle
  * `loc_radius`: polar radius
  * `loc_x`: spherical Cartesian x
  * `loc_y`: spherical Cartesian y
  * `loc_z`: spherical Cartesian z
  * `loc_radius_sph`: spherical radius
  * `loc_theta_sph`: spherical horizontal angle
  * `loc_phi_sph`: spherical azimuth angle

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`

<a id='NeuroAnalyzer.load_study' href='#NeuroAnalyzer.load_study'>#</a>
**`NeuroAnalyzer.load_study`** &mdash; *Function*.



```julia
load_study(file_name)
```

Load `NeuroAnalyzer.STUDY` object from `file_name` file (HDF5-based).

**Arguments**

  * `file_name::String`: name of the file to load

**Returns**

  * `obj::NeuroAnalyzer.STUDY`

<a id='NeuroAnalyzer.save' href='#NeuroAnalyzer.save'>#</a>
**`NeuroAnalyzer.save`** &mdash; *Function*.



```julia
save(obj; file_name, overwrite)
```

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5-based).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `file_name::String`: name of the file to save to
  * `overwrite::Bool=false`

<a id='NeuroAnalyzer.save_study' href='#NeuroAnalyzer.save_study'>#</a>
**`NeuroAnalyzer.save_study`** &mdash; *Function*.



```julia
save_study(obj; file_name, overwrite)
```

Save `NeuroAnalyzer.STUDY` object to `file_name` file (HDF5-based).

**Arguments**

  * `obj::NeuroAnalyzer.STUDY`
  * `file_name::String`: name of the file to save to
  * `overwrite::Bool=false`


<a id='Edit'></a>

<a id='Edit-1'></a>

## Edit

<a id='NeuroAnalyzer.add_channel' href='#NeuroAnalyzer.add_channel'>#</a>
**`NeuroAnalyzer.add_channel`** &mdash; *Function*.



```julia
add_channel(obj; data, label, type)
```

Add channel(s) data to empty `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `data::Array{<:Number, 3}`: channel(s) data
  * `label::Union{String, Vector{String}}=string.(_c(size(data, 1)))`: channel(s) label(s)
  * `type::Union{String, Vector{String}}`: channel(s) type(s)

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_channel!' href='#NeuroAnalyzer.add_channel!'>#</a>
**`NeuroAnalyzer.add_channel!`** &mdash; *Function*.



```julia
add_channel!(obj; data, label, type)
```

Add channel(s) data to empty `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `data::Array{<:Number, 3}`: channel(s) data
  * `label::Union{String, Vector{String}}=string.(_c(size(data, 1)))`: channel(s) label(s)
  * `type::Union{String, Vector{String}}`: channel(s) type(s)

<a id='NeuroAnalyzer.add_labels' href='#NeuroAnalyzer.add_labels'>#</a>
**`NeuroAnalyzer.add_labels`** &mdash; *Function*.



```julia
add_labels(obj; labels)
```

Add channel labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `clabels::Vector{String}`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_labels!' href='#NeuroAnalyzer.add_labels!'>#</a>
**`NeuroAnalyzer.add_labels!`** &mdash; *Function*.



```julia
add_labels!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})
```

Add OBJ channel labels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `clabels::Vector{String}`

<a id='NeuroAnalyzer.add_marker' href='#NeuroAnalyzer.add_marker'>#</a>
**`NeuroAnalyzer.add_marker`** &mdash; *Function*.



```julia
add_marker(obj; id, start, len, desc, ch)
```

Add marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `id::String`: marker ID
  * `start::Real`: marker time in seconds
  * `len::Real=1.0`: marker length in seconds
  * `desc::String`: marker description
  * `ch::Int64=0`: channel number, if 0 then marker is related to all channels

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_marker!' href='#NeuroAnalyzer.add_marker!'>#</a>
**`NeuroAnalyzer.add_marker!`** &mdash; *Function*.



```julia
add_marker!(obj; id, start, len, desc, ch)
```

Add marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `id::String`: marker ID
  * `start::Real`: marker time in seconds
  * `len::Real=1.0`: marker length in seconds
  * `desc::String`: marker description
  * `ch::Int64=0`: channel number, if 0 then marker is related to all channels

<a id='NeuroAnalyzer.channel2marker' href='#NeuroAnalyzer.channel2marker'>#</a>
**`NeuroAnalyzer.channel2marker`** &mdash; *Function*.



```julia
channel2marker(obj; ch, v, id, desc)
```

Convert event channel to markers.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`: event channel number
  * `v::Real=1.0`: event channel value interpreted as an event
  * `id::String`: prefix for marker ID; default is based on event channel name (e.g. "stim1_")
  * `desc::String=""`: marker description; default is based on event channel name (e.g. "stim1")

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.channel2marker!' href='#NeuroAnalyzer.channel2marker!'>#</a>
**`NeuroAnalyzer.channel2marker!`** &mdash; *Function*.



```julia
channel2marker!(obj; ch, v, id, desc)
```

Convert event channel to markers.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`: event channel number
  * `v::Real=1.0`: event channel value interpreted as an event
  * `id::String`: prefix for marker ID; default is "mrk_"
  * `desc::String=""`: prefix for marker description; default is based on event channel name (e.g. "stim1_")

<a id='NeuroAnalyzer.chop' href='#NeuroAnalyzer.chop'>#</a>
**`NeuroAnalyzer.chop`** &mdash; *Function*.



```julia
chop(obj; n)
```

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64=sr(obj)`: number of samples to remove, default is 1 second

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.chop!' href='#NeuroAnalyzer.chop!'>#</a>
**`NeuroAnalyzer.chop!`** &mdash; *Function*.



```julia
chop!(obj; v)
```

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64=sr(obj)`: number of samples to remove, default is 1 second

<a id='NeuroAnalyzer.create_data' href='#NeuroAnalyzer.create_data'>#</a>
**`NeuroAnalyzer.create_data`** &mdash; *Function*.



```julia
create_data(obj; data, fs)
```

Create data, channel labels, types and units and time points for `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `data::Array{Float64, 3}`
  * `fs::Int64`
  * `type::String`: channel types of imported data channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.create_data!' href='#NeuroAnalyzer.create_data!'>#</a>
**`NeuroAnalyzer.create_data!`** &mdash; *Function*.



```julia
create_data!(obj; data, fs)
```

Create data, channel labels, types and units and time points for `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `data::Array{Float64, 3}`
  * `fs::Int64`
  * `type::String`: channel types of imported data channels

<a id='NeuroAnalyzer.create_object' href='#NeuroAnalyzer.create_object'>#</a>
**`NeuroAnalyzer.create_object`** &mdash; *Function*.



```julia
create_object(; data_type)
```

Create an empty `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `data_type::String`: data type of the new object ("eeg", "meg", "nirs", "ecog")

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.create_time' href='#NeuroAnalyzer.create_time'>#</a>
**`NeuroAnalyzer.create_time`** &mdash; *Function*.



```julia
create_time(obj; fs)
```

Create time points vector for `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `fs::Int64`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.create_time!' href='#NeuroAnalyzer.create_time!'>#</a>
**`NeuroAnalyzer.create_time!`** &mdash; *Function*.



```julia
create_time!(obj; fs)
```

Create time points vector for `NeuroAnalyzer.NEURO` object.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `fs::Int64`

<a id='NeuroAnalyzer.delete_channel' href='#NeuroAnalyzer.delete_channel'>#</a>
**`NeuroAnalyzer.delete_channel`** &mdash; *Function*.



```julia
delete_channel(obj; ch)
```

Delete channel(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to be removed
  * `del_opt::Bool=false`: for NIRS data is set as `true` if called from `remove_optode()`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_channel!' href='#NeuroAnalyzer.delete_channel!'>#</a>
**`NeuroAnalyzer.delete_channel!`** &mdash; *Function*.



```julia
delete_channel!(obj; ch)
```

Delete channel(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to be removed
  * `del_opt::Bool=false`: for NIRS data is set as `true` if called from `remove_optode()`

<a id='NeuroAnalyzer.delete_epoch' href='#NeuroAnalyzer.delete_epoch'>#</a>
**`NeuroAnalyzer.delete_epoch`** &mdash; *Function*.



```julia
delete_epoch(obj; ep)
```

Remove epoch(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to be removed

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_epoch!' href='#NeuroAnalyzer.delete_epoch!'>#</a>
**`NeuroAnalyzer.delete_epoch!`** &mdash; *Function*.



```julia
delete_epoch!(obj; ep)
```

Remove epoch(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to be removed

<a id='NeuroAnalyzer.delete_marker' href='#NeuroAnalyzer.delete_marker'>#</a>
**`NeuroAnalyzer.delete_marker`** &mdash; *Function*.



```julia
delete_marker(obj; n)
```

Delete marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64`: marker number

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_marker!' href='#NeuroAnalyzer.delete_marker!'>#</a>
**`NeuroAnalyzer.delete_marker!`** &mdash; *Function*.



```julia
delete_marker!(obj; n)
```

Delete marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64`: marker number

<a id='NeuroAnalyzer.delete_optode' href='#NeuroAnalyzer.delete_optode'>#</a>
**`NeuroAnalyzer.delete_optode`** &mdash; *Function*.



```julia
delete_optode(obj; opt, delete_channels)
```

Delete optode(s) and channels associated with removed optodes.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `opt::Union{Int64, Vector{Int64}, <:AbstractRange}`: optode number(s) to be removed

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.delete_optode!' href='#NeuroAnalyzer.delete_optode!'>#</a>
**`NeuroAnalyzer.delete_optode!`** &mdash; *Function*.



```julia
delete_optode!(obj; opt)
```

Delete optopode(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `opt::Union{Int64, Vector{Int64}, <:AbstractRange}`: optopode number(s) to be removed

<a id='NeuroAnalyzer.detect_bad' href='#NeuroAnalyzer.detect_bad'>#</a>
**`NeuroAnalyzer.detect_bad`** &mdash; *Function*.



```julia
detect_bad(obj; method, ch_t)
```

Detect bad channels and epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p, :var]`: detection method:

      * `:flat`: flat channel(s)
      * `:rmse`: RMSE vs average channel outside of 95%CI
      * `:rmsd`: RMSD
      * `:euclid`: Euclidean distance
      * `:var`: mean signal variance outside of 95%CI and variance inter-quartile outliers
      * `:p2p`: peak-to-peak amplitude; good for detecting transient artifacts
      * `:tkeo`: z-score TKEO value outside of 95%CI
  * `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
  * `ftol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
  * `fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
  * `p::Float64=0.99`: probability threshold (0.0 to 1.0) for marking channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326)
  * `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad
  * `tkeo_method::Symbol=:pow`: method of calculating TKEO, see `tkeo()` for details

**Returns**

Named tuple containing:

  * `bm::Matrix{Bool}`: matrix of bad channels × epochs
  * `be::Vector{Int64}`: list of bad epochs

<a id='NeuroAnalyzer.edit_channel' href='#NeuroAnalyzer.edit_channel'>#</a>
**`NeuroAnalyzer.edit_channel`** &mdash; *Function*.



```julia
edit_channel(obj; ch, field, value)
```

Edit channel properties (`:channel_type` or `:labels`) in `OBJ.header.recording`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`
  * `field::Symbol`
  * `value::Any`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.edit_channel!' href='#NeuroAnalyzer.edit_channel!'>#</a>
**`NeuroAnalyzer.edit_channel!`** &mdash; *Function*.



```julia
edit_channel!(obj; ch, field, value)
```

Edit channel properties (`:channel_type` or `:labels`) in `OBJ.header.recording`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`
  * `field::Symbol`
  * `value::Any`

<a id='NeuroAnalyzer.edit_marker' href='#NeuroAnalyzer.edit_marker'>#</a>
**`NeuroAnalyzer.edit_marker`** &mdash; *Function*.



```julia
edit_marker(obj; n, id, start, len, desc)
```

Edit marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64`: marker number
  * `id::String`: marker ID
  * `start::Real`: marker time in seconds
  * `len::Real=1.0`: marker length in seconds
  * `desc::String`: marker description
  * `ch::Int64`: channel number, if 0 then marker is related to all channels

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.edit_marker!' href='#NeuroAnalyzer.edit_marker!'>#</a>
**`NeuroAnalyzer.edit_marker!`** &mdash; *Function*.



```julia
edit_marker!(obj; n, id, start, len, desc)
```

Edit marker.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64`: marker number
  * `id::String`: marker ID
  * `start::Real`: marker time in seconds
  * `len::Real=1`: marker length in seconds
  * `desc::String`: marker description
  * `ch::Int64`: channel number, if 0 then marker is related to all channels

<a id='NeuroAnalyzer.epoch' href='#NeuroAnalyzer.epoch'>#</a>
**`NeuroAnalyzer.epoch`** &mdash; *Function*.



```julia
epoch(obj; marker, offset, ep_n, ep_len)
```

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `marker::String=""`: marker description to split at
  * `offset::Real=0`: time offset (in seconds) for marker-based epoching (each epoch time will start at `marker time - offset`)
  * `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
  * `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.epoch!' href='#NeuroAnalyzer.epoch!'>#</a>
**`NeuroAnalyzer.epoch!`** &mdash; *Function*.



```julia
epoch!(obj; marker, offset, ep_n, ep_len)
```

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `marker::String=""`: marker description to split at
  * `offset::Real=0`: time offset (in seconds) for marker-based epoching (each epoch time will start at `marker time - offset`)
  * `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
  * `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds

<a id='NeuroAnalyzer.epoch_ts' href='#NeuroAnalyzer.epoch_ts'>#</a>
**`NeuroAnalyzer.epoch_ts`** &mdash; *Function*.



```julia
epoch_ts(obj; ts)
```

Edit epochs time start.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ts::Real`: time start in seconds

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.epoch_ts!' href='#NeuroAnalyzer.epoch_ts!'>#</a>
**`NeuroAnalyzer.epoch_ts!`** &mdash; *Function*.



```julia
epoch_ts!(obj; ts)
```

Edit OBJ epochs time start.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ts::Real`: time start in seconds

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.extract_channel' href='#NeuroAnalyzer.extract_channel'>#</a>
**`NeuroAnalyzer.extract_channel`** &mdash; *Function*.



```julia
extract_channel(obj; ch)
```

Extract channel data.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name

**Returns**

  * `extract_channel::Vector{Float64}`

<a id='NeuroAnalyzer.extract_data' href='#NeuroAnalyzer.extract_data'>#</a>
**`NeuroAnalyzer.extract_data`** &mdash; *Function*.



```julia
extract_data(obj; ch)
```

Extract data.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nepochs(obj)`: index of epochs, default is all epochs
  * `time::Bool=false`: return time vector
  * `etime::Bool=false`: return epoch time vector

**Returns**

  * `signal::Array{Float64, 3}`
  * `time::Vector{Float64}`
  * `etime::Vector{Float64}`

<a id='NeuroAnalyzer.extract_epoch' href='#NeuroAnalyzer.extract_epoch'>#</a>
**`NeuroAnalyzer.extract_epoch`** &mdash; *Function*.



```julia
extract_epoch(obj; ep)
```

Extract epoch.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Int64`: epoch index

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.extract_epoch!' href='#NeuroAnalyzer.extract_epoch!'>#</a>
**`NeuroAnalyzer.extract_epoch!`** &mdash; *Function*.



```julia
extract_epoch!(obj; ep)
```

Extract epoch.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Int64`: epoch index

<a id='NeuroAnalyzer.extract_eptime' href='#NeuroAnalyzer.extract_eptime'>#</a>
**`NeuroAnalyzer.extract_eptime`** &mdash; *Function*.



```julia
extract_eptime(obj)
```

Extract epochs time.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `et::Array{Float64, 3}`

<a id='NeuroAnalyzer.extract_time' href='#NeuroAnalyzer.extract_time'>#</a>
**`NeuroAnalyzer.extract_time`** &mdash; *Function*.



```julia
extract_time(obj)
```

Extract time.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `tpts::Array{Float64, 3}`

<a id='NeuroAnalyzer.get_channel' href='#NeuroAnalyzer.get_channel'>#</a>
**`NeuroAnalyzer.get_channel`** &mdash; *Function*.



```julia
get_channel(obj; ch)
```

Return channel number (if provided by name) or name (if provided by number).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name

**Returns**

  * `ch_idx::Union{Int64, String}`: channel number or name

<a id='NeuroAnalyzer.get_channel_bytype' href='#NeuroAnalyzer.get_channel_bytype'>#</a>
**`NeuroAnalyzer.get_channel_bytype`** &mdash; *Function*.



```julia
get_channel_bytype(obj; type)
```

Return channel number(s) for channel of `type` type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::Union{String, Vector{String}}="all"`: channel type

**Returns**

  * `ch_idx::Vector{Int64}`


```
get_channel_bytype(ct; type)
```

Return channel number(s) for channel of `type` type.

**Arguments**

  * `ct::Vector{String}`: list of channel types
  * `type::Union{String, Vector{String}}="all"`: channel type

**Returns**

  * `ch_idx::Vector{Int64}`

<a id='NeuroAnalyzer.get_channel_bywl' href='#NeuroAnalyzer.get_channel_bywl'>#</a>
**`NeuroAnalyzer.get_channel_bywl`** &mdash; *Function*.



```julia
get_channel_bywl(obj; wl)
```

Return NIRS channel number(s) for wavelength `wl`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `wl::Real`: wavelength (in nm)

**Returns**

  * `ch_idx::Vector{Int64}`

<a id='NeuroAnalyzer.get_channel_type' href='#NeuroAnalyzer.get_channel_type'>#</a>
**`NeuroAnalyzer.get_channel_type`** &mdash; *Function*.



```julia
get_channel_type(obj; ch, type)
```

Get channel type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.join' href='#NeuroAnalyzer.join'>#</a>
**`NeuroAnalyzer.join`** &mdash; *Function*.



```julia
join(obj1, obj2)
```

Join two NeuroAnalyzer objects. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.join!' href='#NeuroAnalyzer.join!'>#</a>
**`NeuroAnalyzer.join!`** &mdash; *Function*.



```julia
join!(obj1, obj2)
```

Join two NeuroAnalyzer objects. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.keep_channel' href='#NeuroAnalyzer.keep_channel'>#</a>
**`NeuroAnalyzer.keep_channel`** &mdash; *Function*.



```julia
keep_channel(obj; ch)
```

Keep channel(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to keep

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.keep_channel!' href='#NeuroAnalyzer.keep_channel!'>#</a>
**`NeuroAnalyzer.keep_channel!`** &mdash; *Function*.



```julia
keep_channel!(obj; ch)
```

Keep channel(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to keep

<a id='NeuroAnalyzer.keep_channel_type' href='#NeuroAnalyzer.keep_channel_type'>#</a>
**`NeuroAnalyzer.keep_channel_type`** &mdash; *Function*.



```julia
keep_channel_type(obj; type)
```

Keep channel(s) of `type` type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::String="eeg"`: type of channels to keep

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.keep_channel_type!' href='#NeuroAnalyzer.keep_channel_type!'>#</a>
**`NeuroAnalyzer.keep_channel_type!`** &mdash; *Function*.



```julia
keep_channel_type!(obj; type)
```

Keep OBJ channels of `type` type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::String="eeg"`: type of channels to keep

<a id='NeuroAnalyzer.keep_epoch' href='#NeuroAnalyzer.keep_epoch'>#</a>
**`NeuroAnalyzer.keep_epoch`** &mdash; *Function*.



```julia
keep_epoch(obj; ep)
```

Keep epoch(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to keep

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.keep_epoch!' href='#NeuroAnalyzer.keep_epoch!'>#</a>
**`NeuroAnalyzer.keep_epoch!`** &mdash; *Function*.



```julia
keep_epoch!(obj; ep)
```

Keep OBJ epoch(s).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to keep

<a id='NeuroAnalyzer.reflect' href='#NeuroAnalyzer.reflect'>#</a>
**`NeuroAnalyzer.reflect`** &mdash; *Function*.



```julia
reflect(obj; n)
```

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64=sr(obj)`: number of samples to add, default is 1 second

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reflect!' href='#NeuroAnalyzer.reflect!'>#</a>
**`NeuroAnalyzer.reflect!`** &mdash; *Function*.



```julia
reflect!(obj; n)
```

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `n::Int64=sr(obj)`: number of samples to add, default is 1 second

<a id='NeuroAnalyzer.rename_channel' href='#NeuroAnalyzer.rename_channel'>#</a>
**`NeuroAnalyzer.rename_channel`** &mdash; *Function*.



```julia
rename_channel(obj; ch, name)
```

Rename channel.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name
  * `name::String`: new name

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.rename_channel!' href='#NeuroAnalyzer.rename_channel!'>#</a>
**`NeuroAnalyzer.rename_channel!`** &mdash; *Function*.



```julia
rename_channel!(obj; ch, name)
```

Rename channel.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name
  * `name::String`: new name

<a id='NeuroAnalyzer.replace_channel' href='#NeuroAnalyzer.replace_channel'>#</a>
**`NeuroAnalyzer.replace_channel`** &mdash; *Function*.



```julia
replace_channel(obj; ch, s)
```

Replace channel.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name
  * `s::Array{Float64, 3}`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.replace_channel!' href='#NeuroAnalyzer.replace_channel!'>#</a>
**`NeuroAnalyzer.replace_channel!`** &mdash; *Function*.



```julia
replace_channel!(obj; ch, s)
```

Replace channel.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name
  * `s::Array{Float64, 3}`: signal to replace with

<a id='NeuroAnalyzer.set_channel_type' href='#NeuroAnalyzer.set_channel_type'>#</a>
**`NeuroAnalyzer.set_channel_type`** &mdash; *Function*.



```julia
set_channel_type(obj; ch, type)
```

Set channel type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`
  * `type::String`

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.set_channel_type!' href='#NeuroAnalyzer.set_channel_type!'>#</a>
**`NeuroAnalyzer.set_channel_type!`** &mdash; *Function*.



```julia
set_channel_type!(obj; ch, new_name)
```

Set channel type.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`
  * `type::String`

<a id='NeuroAnalyzer.signal_channels' href='#NeuroAnalyzer.signal_channels'>#</a>
**`NeuroAnalyzer.signal_channels`** &mdash; *Function*.



```julia
signal_channels(obj)
```

Return all signal (e.g. EEG or MEG) channels; signal is determined by `:data_type` variable in `obj.header.recording`). For MEG data type, 'meg', `grad` and `mag` channels are returned.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `chs::Vector{Int64}`


```
signal_channels(dt, ct)
```

Return all signal (e.g. EEG or MEG) channels for a given data type. For MEG data type, 'meg', `grad` and `mag` channels are returned.

**Arguments**

  * `dt::String`: data type of an recording object
  * `ct::Vector{String}`: list of channel types

**Returns**

  * `chs::Vector{Int64}`

<a id='NeuroAnalyzer.trim' href='#NeuroAnalyzer.trim'>#</a>
**`NeuroAnalyzer.trim`** &mdash; *Function*.



```julia
trim(s; seg, inverse)
```

Remove segment from the signal.

**Arguments**

  * `v::AbstractVector`
  * `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
  * `inverse::Bool=false`: if true, keep the segment

**Returns**

  * `trim::Vector{Float64}`


```
trim(m; seg, inverse)
```

Remove segment from the signal.

**Arguments**

  * `m::AbstractMatrix`
  * `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
  * `inverse::Bool=false`: if true, keep the segment

**Returns**

  * `trim::Array{Float64}`


```
trim(a; seg, inverse)
```

Remove segment from the signal.

**Arguments**

  * `a::AbstractArray`
  * `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
  * `inverse::Bool=false`: if true, keep the segment

**Returns**

  * `trim::Array{Float64}`


```
trim(obj; seg, inverse, remove_epochs)
```

Trim signal by removing parts of the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
  * `inverse::Bool=false`: if true, keep the segment
  * `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.trim!' href='#NeuroAnalyzer.trim!'>#</a>
**`NeuroAnalyzer.trim!`** &mdash; *Function*.



```julia
trim!(obj; seg, inverse, remove_epochs)
```

Trim signal by removing parts of the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
  * `inverse::Bool=false`: if true, keep the segment
  * `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

<a id='NeuroAnalyzer.vch' href='#NeuroAnalyzer.vch'>#</a>
**`NeuroAnalyzer.vch`** &mdash; *Function*.



```julia
vch(obj; f)
```

Calculate virtual channel using formula `f`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the OBJ object

**Returns**

  * `vc::Array{Float64, 3}`: single channel × time × epochs

<a id='NeuroAnalyzer.view_marker' href='#NeuroAnalyzer.view_marker'>#</a>
**`NeuroAnalyzer.view_marker`** &mdash; *Function*.



```julia
view_marker(obj)
```

Show markers.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`


<a id='Process'></a>

<a id='Process-1'></a>

## Process

<a id='NeuroAnalyzer.add_signal' href='#NeuroAnalyzer.add_signal'>#</a>
**`NeuroAnalyzer.add_signal`** &mdash; *Function*.



```julia
add_signal(s1, s2)
```

Add signal.

**Arguments**

  * `s1::AbstractVector`: target signal
  * `s2::AbstractVector`: signal to be added

**Returns**

  * `s_noisy::AbstractVector`


```
add_signal(obj; ch, s)
```

Add signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `s::AbstractVector`: signal to be added to each channel

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.add_signal!' href='#NeuroAnalyzer.add_signal!'>#</a>
**`NeuroAnalyzer.add_signal!`** &mdash; *Function*.



```julia
add_signal!(obj; ch, n)
```

Add signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `s::AbstractVector`: signal to be added to each channel

<a id='NeuroAnalyzer.average' href='#NeuroAnalyzer.average'>#</a>
**`NeuroAnalyzer.average`** &mdash; *Function*.



```julia
average(s)
```

Average all channels.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `average::AbstractArray`


```
average(s1, s2)
```

Averages two signals.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`

**Returns**

  * `average::Vector{Float64}`


```
average(obj; ch)
```

Return the average signal of channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
average(obj1, obj2)
```

Return the average signal of all `obj1` and `obj2` channels.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.average!' href='#NeuroAnalyzer.average!'>#</a>
**`NeuroAnalyzer.average!`** &mdash; *Function*.



```julia
average!(obj; ch)
```

Return the average signal of channels.  

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

<a id='NeuroAnalyzer.bpsplit' href='#NeuroAnalyzer.bpsplit'>#</a>
**`NeuroAnalyzer.bpsplit`** &mdash; *Function*.



```julia
bpsplit(obj; ch, order, window)
```

Split signal into frequency bands using IIR band-pass filter.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `order::Int64=8`: number of taps for FIR band-pass filter
  * `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb

**Returns**

Named tuple containing:

  * `s::Array{Float64, 4}`: split signal
  * `bn::Vector{Symbol}`: band names
  * `bf::Vector{Tuple{Real, Real}}`: band frequencies

<a id='NeuroAnalyzer.cbp' href='#NeuroAnalyzer.cbp'>#</a>
**`NeuroAnalyzer.cbp`** &mdash; *Function*.



```julia
cbp(s; pad, frq, fs)
```

Perform convolution band-pass filtering.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64`: pad with `pad` zeros
  * `frq::Real`: filter frequency
  * `fs::Int64`: sampling rate

**Returns**

  * `cbp::Vector{Float64}`


```
cbp(obj; ch, pad, frq)
```

Perform convolution bandpass filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.cbp!' href='#NeuroAnalyzer.cbp!'>#</a>
**`NeuroAnalyzer.cbp!`** &mdash; *Function*.



```julia
cbp!(obj; ch, pad, frq)
```

Perform convolution bandpass filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Tuple{Real, Real}`: filter frequency

<a id='NeuroAnalyzer.ch_zero' href='#NeuroAnalyzer.ch_zero'>#</a>
**`NeuroAnalyzer.ch_zero`** &mdash; *Function*.



```julia
ch_zero(obj)
```

Zero channels at the beginning and at the end.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.ch_zero!' href='#NeuroAnalyzer.ch_zero!'>#</a>
**`NeuroAnalyzer.ch_zero!`** &mdash; *Function*.



```julia
ch_zero!(obj)
```

Zero channels at the beginning and at the end.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.csd' href='#NeuroAnalyzer.csd'>#</a>
**`NeuroAnalyzer.csd`** &mdash; *Function*.



```julia
csd(obj; m, n, lambda)
```

Transform data using Current Source Density (CSD) transformation based on spherical spline surface Laplacian.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
  * `n::Int64=8`: Legendre polynomial order
  * `lambda::Float64=10^-5`: smoothing factor

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: with `csd` channel types and `µV/m²` units

**Source**

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-187 Kayser J, Tenke CE. Principal components analysis of Laplacian waveforms as a generic method for identifying ERP generator patterns: I. Evaluation with auditory oddball tasks. Clin Neurophysiol 2006;117(2):348-368

<a id='NeuroAnalyzer.csd!' href='#NeuroAnalyzer.csd!'>#</a>
**`NeuroAnalyzer.csd!`** &mdash; *Function*.



```julia
csd!(obj; m, n, lambda)
```

Transform data using Current Source Density (CSD) transformation based on spherical spline surface Laplacian.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
  * `n::Int64=8`: Legendre polynomial order
  * `lambda::Float64=10^-5`: smoothing factor

**Returns**

  * `G::Matrix{Float64}`: transformation matrix (SP spline)
  * `H::Matrix{Float64}`: transformation matrix (CSD spline)

**Source**

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7

<a id='NeuroAnalyzer.cw_trans' href='#NeuroAnalyzer.cw_trans'>#</a>
**`NeuroAnalyzer.cw_trans`** &mdash; *Function*.



```julia
cw_trans(s; wt, type, l)
```

Perform continuous wavelet transformation (CWT).

**Arguments**

  * `s::AbstractVector`
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

  * `ct::Array{Float64, 2}`: CWT coefficients (by rows)


```
cw_trans(s; wt)
```

Perform continuous wavelet transformation (CWT).

**Arguments**

  * `s::AbstractArray`
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

  * `ct::Array{Float64, 4}`: CWT coefficients (by rows)


```
cw_trans(obj; ch, wt)
```

Perform continuous wavelet transformation (CWT).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

  * `ct::Array{Float64, 4}`: CWT coefficients (by rows)

<a id='NeuroAnalyzer.denoise_fft' href='#NeuroAnalyzer.denoise_fft'>#</a>
**`NeuroAnalyzer.denoise_fft`** &mdash; *Function*.



denoise_fft(s; pad, t)

Perform FFT denoising.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64=0`: number of zeros to add
  * `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

**Returns**

Named tuple containing:

  * `s::Vector{Float64}`
  * `f_idx::BitVector`: index of components zeroed


```
denoise_fft(s; pad, t)
```

Perform FFT denoising.

**Arguments**

  * `s::AbstractArray`
  * `pad::Int64=0`: number of zeros to add
  * `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
denoise_fft(obj; ch, pad, t)
```

Perform FFT denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `t::Int64=100`: PSD threshold for keeping frequency components

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.denoise_fft!' href='#NeuroAnalyzer.denoise_fft!'>#</a>
**`NeuroAnalyzer.denoise_fft!`** &mdash; *Function*.



```julia
denoise_fft!(obj; ch, pad, t)
```

Perform FFT denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `t::Int64=100`: PSD threshold for keeping frequency components

<a id='NeuroAnalyzer.denoise_wavelet' href='#NeuroAnalyzer.denoise_wavelet'>#</a>
**`NeuroAnalyzer.denoise_wavelet`** &mdash; *Function*.



```julia
denoise_wavelet(s; wt)
```

Perform wavelet denoising.

**Arguments**

  * `s::AbstractVector`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

**Returns**

  * `s_new::Vector{Float64}`


```
denoise_wavelet(s; wt)
```

Perform wavelet denoising.

**Arguments**

  * `s::AbstractArray`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
denoise_wavelet(obj; ch, wt)
```

Perform wavelet denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.denoise_wavelet!' href='#NeuroAnalyzer.denoise_wavelet!'>#</a>
**`NeuroAnalyzer.denoise_wavelet!`** &mdash; *Function*.



```julia
denoise_wavelet!(obj; ch, wt)
```

Perform wavelet denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

<a id='NeuroAnalyzer.denoise_wien' href='#NeuroAnalyzer.denoise_wien'>#</a>
**`NeuroAnalyzer.denoise_wien`** &mdash; *Function*.



```julia
denoise_wien(s)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::Vector{Float64}`


```
denoise_wien(obj; ch)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.denoise_wien!' href='#NeuroAnalyzer.denoise_wien!'>#</a>
**`NeuroAnalyzer.denoise_wien!`** &mdash; *Function*.



```julia
denoise_wien!(obj; ch)
```

Perform Wiener deconvolution denoising.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

<a id='NeuroAnalyzer.derivative' href='#NeuroAnalyzer.derivative'>#</a>
**`NeuroAnalyzer.derivative`** &mdash; *Function*.



```julia
derivative(s)
```

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `s_new::AbstractVector`


```
derivative(s)
```

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::Array{Float64, 3}`


```
derivative(obj; ch)
```

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.derivative!' href='#NeuroAnalyzer.derivative!'>#</a>
**`NeuroAnalyzer.derivative!`** &mdash; *Function*.



```julia
derivative!(obj; ch)
```

Return derivative (calculated using symmetric difference quotient) of a discrete signal of the same length.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

<a id='NeuroAnalyzer.detect_powerline' href='#NeuroAnalyzer.detect_powerline'>#</a>
**`NeuroAnalyzer.detect_powerline`** &mdash; *Function*.



```julia
detect_powerline(s; fs)
```

Detect power line noise frequency.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate

**Returns**

  * `noise_frq::Float64`: peak noise frequency in Hz


```
detect_powerline(obj)
```

Detect power line noise frequency.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `noise_frq::Array{Float64, 3}`: peak noise frequency in Hz

<a id='NeuroAnalyzer.detrend' href='#NeuroAnalyzer.detrend'>#</a>
**`NeuroAnalyzer.detrend`** &mdash; *Function*.



```julia
detrend(s; type, offset, order, span, fs)
```

Perform piecewise detrending.

**Arguments**

  * `s::AbstractVector`
  * `type::Symbol=:linear`:

      * `:loess`: fit loess approximation and subtract it from `s`
      * `:poly`: polynomial of `order` is subtracted from `s`
      * `:mean`: the mean of `s` is subtracted from `s`
      * `:constant`: `offset` is subtracted from `s`
      * `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
      * `:linear`: linear trend (1st order polynomial) is subtracted from `s`
  * `offset::Real=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

**Returns**

  * `s_new::Vector{Float64}`


```
detrend(s; type, offset, order, f)
```

Perform piecewise detrending.

**Arguments**

  * `s::AbstractArray`
  * `type::Symbol=:linear`: detrending method

      * `:loess`: fit loess approximation and subtract it from `s`
      * `:poly`: polynomial of `order` is subtracted from `s`
      * `:mean`: the mean of `s` is subtracted from `s`
      * `:constant`: `offset` is subtracted from `s`
      * `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
      * `:linear`: linear trend (1st order polynomial) is subtracted from `s`
  * `offset::Real=0`: constant for `:constant` detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

**Returns**

  * `s_new::Array`{Float64, 3}


```
detrend(obj; ch, type, offset, order, f)
```

Perform piecewise detrending.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `type::Symbol=:linear`: detrending method

      * `:loess`: fit loess approximation and subtract it from `s`
      * `:poly`: polynomial of `order` is subtracted from `s`
      * `:mean`: the mean of `s` is subtracted from `s`
      * `:constant`: `offset` is subtracted from `s`
      * `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
      * `:linear`: linear trend (1st order polynomial) is subtracted from `s`
  * `offset::Real=0`: constant for `:constant` detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.detrend!' href='#NeuroAnalyzer.detrend!'>#</a>
**`NeuroAnalyzer.detrend!`** &mdash; *Function*.



```julia
detrend!(obj; ch, type, offset, order, span)
```

Perform piecewise detrending.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `type::Symbol=:linear`: detrending method

      * `:loess`: fit loess approximation and subtract it from `s`
      * `:poly`: polynomial of `order` is subtracted from `s`
      * `:mean`: the mean of `s` is subtracted from `s`
      * `:constant`: `offset` is subtracted from `s`
      * `:ls`: the result of a linear least-squares fit to `s` is subtracted from `s`
      * `:linear`: linear trend (1st order polynomial) is subtracted from `s`
  * `offset::Real=0`: constant for :constant detrending
  * `order::Int64=1`: polynomial fitting order
  * `f::Float64=1.0`: smoothing factor for `:loess` or frequency for `:hp`

<a id='NeuroAnalyzer.downsample' href='#NeuroAnalyzer.downsample'>#</a>
**`NeuroAnalyzer.downsample`** &mdash; *Function*.



```julia
downsample(obj; new_sr)
```

Downsample.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.downsample!' href='#NeuroAnalyzer.downsample!'>#</a>
**`NeuroAnalyzer.downsample!`** &mdash; *Function*.



```julia
downsample!(obj; new_sr)
```

Downsample.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroAnalyzer.dw_trans' href='#NeuroAnalyzer.dw_trans'>#</a>
**`NeuroAnalyzer.dw_trans`** &mdash; *Function*.



```julia
dw_trans(s; wt, type, l)
```

Perform discrete wavelet transformation (DWT).

**Arguments**

  * `s::AbstractVector`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: 

      * `:sdwt`: Stationary Wavelet Transforms
      * `:acdwt`: Autocorrelation Wavelet Transforms
  * `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

**Returns**

  * `dt::Array{Float64, 2}`: DWT coefficients cAl, cD1, ..., cDl (by rows)


```
dw_trans(s; wt, type, l)
```

Perform discrete wavelet transformation (DWT).

**Arguments**

  * `s::AbstractArray`
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: 

      * `:sdwt`: Stationary Wavelet Transforms
      * `:acdwt`: Autocorrelation Wavelet Transforms
  * `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

**Returns**

  * `dt::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)


```
dw_trans(obj; ch, wt, type, l)
```

Perform discrete wavelet transformation (DWT).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: 

      * `:sdwt`: Stationary Wavelet Transforms
      * `:acdwt`: Autocorrelation Wavelet Transforms
  * `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

**Returns**

  * `dt::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)

<a id='NeuroAnalyzer.dwtsplit' href='#NeuroAnalyzer.dwtsplit'>#</a>
**`NeuroAnalyzer.dwtsplit`** &mdash; *Function*.



```julia
dwtsplit(obj; ch, wt, type, n)
```

Split into bands using discrete wavelet transformation (DWT).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64}`: channel number
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: 

      * `:sdwt`: Stationary Wavelet Transforms
      * `:acdwt`: Autocorrelation Wavelet Transforms
  * `n::Int64=0`: number of bands, default is maximum number of bands available or total transformation

**Returns**

  * `b::Array{Float64, 4}`: bands from lowest to highest frequency (by rows)

<a id='NeuroAnalyzer.edit_montage' href='#NeuroAnalyzer.edit_montage'>#</a>
**`NeuroAnalyzer.edit_montage`** &mdash; *Function*.



```julia
edit_montage(file_name; data_format, detect_type)
```

Edit montage file in the OS editor.

**Arguments**

  * `file_name::String`: name of the file to load

<a id='NeuroAnalyzer.erp' href='#NeuroAnalyzer.erp'>#</a>
**`NeuroAnalyzer.erp`** &mdash; *Function*.



```julia
erp(obj; bl)
```

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `bl::Real=0`: baseline is the first `bl` seconds; if `bl` is greater than 0, DC value is calculated as mean of the first `bl` seconds and subtracted from the signal

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.erp!' href='#NeuroAnalyzer.erp!'>#</a>
**`NeuroAnalyzer.erp!`** &mdash; *Function*.



```julia
erp!(obj; bl)
```

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `bl::Real=0`: baseline is the first `bl` seconds; if `bl` is greater than 0, DC value is calculated as mean of the first `n` samples and subtracted from the signal

<a id='NeuroAnalyzer.fconv' href='#NeuroAnalyzer.fconv'>#</a>
**`NeuroAnalyzer.fconv`** &mdash; *Function*.



```julia
fconv(s; kernel, norm)
```

Perform convolution in the frequency domain.

**Arguments**

  * `s::AbstractVector`
  * `kernel::AbstractVector`
  * `pad::Int64=0`: number of zeros to add
  * `norm::Bool=true`: normalize kernel

**Returns**

  * `s_new::Vector{Float64}`


```
fconv(s; kernel, norm)
```

Perform convolution in the frequency domain.

**Arguments**

  * `s::AbstractArray`
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
  * `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
  * `pad::Int64=0`: number of zeros to add signal for FFT

**Returns**

  * `s_new::Array{Float64, 3}`: convoluted signal


```
fconv(obj; ch, kernel, norm)
```

Perform convolution in the frequency domain.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
  * `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
  * `pad::Int64=0`: number of zeros to add signal for FFT

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: convoluted signal

<a id='NeuroAnalyzer.fconv!' href='#NeuroAnalyzer.fconv!'>#</a>
**`NeuroAnalyzer.fconv!`** &mdash; *Function*.



```julia
fconv!(obj; ch)
```

Perform convolution in the frequency domain.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
  * `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
  * `pad::Int64=0`: number of zeros to add signal for FFT

<a id='NeuroAnalyzer.filter' href='#NeuroAnalyzer.filter'>#</a>
**`NeuroAnalyzer.filter`** &mdash; *Function*.



```julia
filter(obj; <keyword arguments>)
```

Apply filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
  * `ftype::Union{Nothing, Symbol}=nothing`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez` filters
  * `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):

      * `:twopass`
      * `:onepass`
      * `:reverse`: one pass, reverse direction
  * `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter
  * `window::Union{Nothing, AbstractVector, Int64}=nothing`: window length for `:remez` and `:fir` filters
  * `preview::Bool=false`: plot filter response

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

**Notes**

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.

<a id='NeuroAnalyzer.filter!' href='#NeuroAnalyzer.filter!'>#</a>
**`NeuroAnalyzer.filter!`** &mdash; *Function*.



```julia
filter!(obj; <keyword arguments>)
```

Apply filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
  * `ftype::Union{Nothing, Symbol}=nothing`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and `:remez` filters
  * `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):

      * `:twopass`
      * `:onepass`
      * `:reverse`: one pass, reverse direction
  * `order::Int64=8`: filter order (6 dB/octave) for IIR filters, number of taps for `:remez` filter, attenuation (× 4 dB) for `:fir` filter
  * `window::Union{Nothing, AbstractVector, Int64}=nothing`: window length for `:remez` and `:fir` filters
  * `preview::Bool=false`: plot filter response

**Notes**

For `:poly` filter `order` and `window` have to be set experimentally, recommended initial values are: `order=4` and `window=32`.

<a id='NeuroAnalyzer.filter_apply' href='#NeuroAnalyzer.filter_apply'>#</a>
**`NeuroAnalyzer.filter_apply`** &mdash; *Function*.



```julia
filter_apply(s; <keyword arguments>)
```

Apply IIR or FIR filter.

**Arguments**

  * `s::AbstractVector`
  * `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`: filter
  * `dir:Symbol=:twopass`: filter direction (for causal filter use `:onepass`):

      * `:twopass`
      * `:onepass`
      * `:reverse`: one pass, reverse direction

**Returns**

  * `s_filtered::Vector{Float64}`

<a id='NeuroAnalyzer.filter_create' href='#NeuroAnalyzer.filter_create'>#</a>
**`NeuroAnalyzer.filter_create`** &mdash; *Function*.



```julia
filter_create(; <keyword arguments>)
```

Create IIR or FIR filter.

**Arguments**

  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
  * `ftype::Union{Nothing, Symbol}=nothing`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}=0`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `n::Int64`: signal length in samples
  * `fs::Int64`: sampling rate
  * `order::Int64=8`: filter order (6 dB/octave), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
  * `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb

**Returns**

  * `flt::Union{Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`

<a id='NeuroAnalyzer.filter_g' href='#NeuroAnalyzer.filter_g'>#</a>
**`NeuroAnalyzer.filter_g`** &mdash; *Function*.



```julia
filter_g(s, fs, pad, f, gw)
```

Filter using Gaussian in the frequency domain.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `pad::Int64=0`: number of zeros to add
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

**Returns**

Named tuple containing:

  * `s_filtered::Vector{Float64}`


```
filter_g(s; pad, f, gw)
```

Filter using Gaussian in the frequency domain.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `pad::Int64=0`: number of zeros to add
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

**Returns**

  * `s_filtered::NeuroAnalyzer.NEURO`


```
filter_g(obj; ch, pad, f, gw)
```

Filter using Gaussian in the frequency domain.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.filter_g!' href='#NeuroAnalyzer.filter_g!'>#</a>
**`NeuroAnalyzer.filter_g!`** &mdash; *Function*.



```julia
filter_g!(obj; ch, pad, f, gw)
```

Filter using Gaussian in the frequency domain.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64=0`: number of zeros to add
  * `f::Real`: filter frequency
  * `gw::Real=5`: Gaussian width in Hz

<a id='NeuroAnalyzer.filter_mavg' href='#NeuroAnalyzer.filter_mavg'>#</a>
**`NeuroAnalyzer.filter_mavg`** &mdash; *Function*.



```julia
filter_mavg(s; <keyword arguments>)
```

Filter using moving average filter (with threshold).

**Arguments**

  * `s::AbstractVector`
  * `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `s_filtered::Vector{Float64}`


```
filter_mavg(s; k, t, window)
```

Filter using moving average filter (with threshold).

**Arguments**

  * `s::AbstractArray`
  * `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `s_filtered::Array{Float64, 3}`: convoluted signal


```
filter_mavg(obj; ch, k, t, window)
```

Filter using moving average filter (with threshold).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: convoluted signal

**Source**

1. https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter

<a id='NeuroAnalyzer.filter_mavg!' href='#NeuroAnalyzer.filter_mavg!'>#</a>
**`NeuroAnalyzer.filter_mavg!`** &mdash; *Function*.



```julia
filter_mavg!(obj; ch, k, t, window)
```

Filter using moving average filter (with threshold).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`)
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

<a id='NeuroAnalyzer.filter_mmed' href='#NeuroAnalyzer.filter_mmed'>#</a>
**`NeuroAnalyzer.filter_mmed`** &mdash; *Function*.



```julia
filter_mmed(s; <keyword arguments>)
```

Filter using moving median filter (with threshold).

**Arguments**

  * `s::AbstractVector`
  * `k::Int64=8`: window length is `2 × k + 1`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `s_filtered::Vector{Float64}`


```
filter_mmed(s; k, t, window)
```

Filter using moving median filter (with threshold).

**Arguments**

  * `s::AbstractArray`
  * `k::Int64=8`: window length is `2 × k + 1`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `s_filtered::Array{Float64, 3}`: convoluted signal


```
filter_mmed(obj; ch, k, t, window)
```

Filter using moving median filter (with threshold).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `k::Int64=8`: window length is `2 × k + 1`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: convoluted signal

<a id='NeuroAnalyzer.filter_mmed!' href='#NeuroAnalyzer.filter_mmed!'>#</a>
**`NeuroAnalyzer.filter_mmed!`** &mdash; *Function*.



```julia
filter_mmed!(obj; ch, k, t, window)
```

Filter using moving median filter (with threshold).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `k::Int64=8`: window length is `2 × k + 1`
  * `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
  * `window::Union{Nothing, AbstractVector}=nothing`: weighting window

<a id='NeuroAnalyzer.filter_poly' href='#NeuroAnalyzer.filter_poly'>#</a>
**`NeuroAnalyzer.filter_poly`** &mdash; *Function*.



```julia
filter_poly(s; order, window)
```

Filter using polynomial filter.

**Arguments**

  * `s::AbstractVector`
  * `order::Int64=8`: polynomial order
  * `window::Int64=10`: window length

**Returns**

  * `s_filtered::Vector{Float64}`


```
filter_poly(s; order, window)
```

Filter using polynomial filter.

**Arguments**

  * `s::AbstractArray`
  * `order::Int64=8`: polynomial order
  * `window::Int64=10`: window length

**Returns**

  * `s_filtered::Array{Float64, 3}`: convoluted signal


```
filter_poly(obj; ch, order, window)
```

Filter using polynomial filter.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `order::Int64=8`: polynomial order
  * `window::Int64=10`: window length

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: convoluted signal

<a id='NeuroAnalyzer.filter_poly!' href='#NeuroAnalyzer.filter_poly!'>#</a>
**`NeuroAnalyzer.filter_poly!`** &mdash; *Function*.



```julia
filter_poly!(obj; ch, order, window)
```

Filter using polynomial filter.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `order::Int64=8`: polynomial order
  * `window::Int64=10`: window length

<a id='NeuroAnalyzer.filter_sg' href='#NeuroAnalyzer.filter_sg'>#</a>
**`NeuroAnalyzer.filter_sg`** &mdash; *Function*.



```julia
filter_sg(s; order, window)
```

Filter using Savitzky-Golay filter.

**Arguments**

  * `s::AbstractVector`
  * `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
  * `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

**Returns**

  * `s_filtered::Vector{Float64}`


```
filter_sg(s; order, window)
```

Filter using Savitzky-Golay filter.

**Arguments**

  * `s::AbstractArray`
  * `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
  * `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

**Returns**

  * `s_filtered::Array{Float64, 3}`: convoluted signal


```
filter_sg(obj; ch, order, window)
```

Filter using Savitzky-Golay filter.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
  * `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: convoluted signal

<a id='NeuroAnalyzer.filter_sg!' href='#NeuroAnalyzer.filter_sg!'>#</a>
**`NeuroAnalyzer.filter_sg!`** &mdash; *Function*.



```julia
filter_sg!(obj; ch, order, window)
```

Filter using Savitzky-Golay filter.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
  * `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

<a id='NeuroAnalyzer.gh' href='#NeuroAnalyzer.gh'>#</a>
**`NeuroAnalyzer.gh`** &mdash; *Function*.



```julia
gh(locs; m, n)
```

Generate G and H matrices.

**Arguments**

  * `locs::DataFrame`
  * `m::Int64=4`: positive integer constant that affects spherical spline flexibility, higher `m` values mean increased rigidity
  * `n::Int64=8`: Legendre polynomial order

**Returns**

Named tuple containing:

  * `G::Matrix{Float64}`: transformation matrix (SP spline)
  * `H::Matrix{Float64}`: transformation matrix (CSD spline)

**Source**

Perrin F, Pernier J, Bertrand O, Echallier JF. Spherical splines for scalp potential and current density mapping. Electroencephalography and Clinical Neurophysiology. 1989;72(2):184-7

<a id='NeuroAnalyzer.ica_decompose' href='#NeuroAnalyzer.ica_decompose'>#</a>
**`NeuroAnalyzer.ica_decompose`** &mdash; *Function*.



```julia
ica_decompose(s; <keyword arguments>)
```

Calculate `n` first Independent Components using FastICA algorithm.

**Arguments**

  * `s::AbstractMatrix`
  * `n::Int64`: number of ICs
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor:

      * `:tanh`
      * `:gaus`

**Returns**

Named tuple containing:

  * `ic::Matrix{Float64}`: components IC(1)..IC(n) (W * data), components are sorted by decreasing variance
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n) (inv(W))


```
ica_decompose(obj; <keyword arguments>)
```

Perform independent component analysis (ICA) using FastICA algorithm.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `n::Int64=length(ch)`: number of ICs, default is the number of channels
  * `iter::Int64=100`: maximum number of iterations
  * `f::Symbol=:tanh`: neg-entropy functor:

      * `:tanh`
      * `:gaus`

**Returns**

Named tuple containing:

  * `ic::Matrix{Float64}`: components IC(1)..IC(n) (W * data), components are sorted by decreasing variance
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n) (inv(W))
  * `ic_var::Vector{Float64}`: variance of components

<a id='NeuroAnalyzer.ica_reconstruct' href='#NeuroAnalyzer.ica_reconstruct'>#</a>
**`NeuroAnalyzer.ica_reconstruct`** &mdash; *Function*.



```julia
ica_reconstruct(; ic, ic_mw, ic_idx)
```

Reconstruct signal using ICA components.

**Arguments**

  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
  * `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

**Returns**

  * `s_new::Matrix{Float64}`: reconstructed signal


```
ica_reconstruct(obj; ch, ic_idx, keep)
```

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}`: list of ICs to remove or keep or keep
  * `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
ica_reconstruct(obj, ic, ic_mw; ch, ic_idx, keep)
```

Reconstruct signals using external ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
  * `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.ica_reconstruct!' href='#NeuroAnalyzer.ica_reconstruct!'>#</a>
**`NeuroAnalyzer.ica_reconstruct!`** &mdash; *Function*.



```julia
ica_reconstruct!(obj; ch, ic_idx, keep)
```

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
  * `keep::Bool=false`: if `true`, then the ICs are kept instead of removed


```
ica_reconstruct!(obj, ic, ic_mw; ch, ic_idx, keep)
```

Reconstruct signals using external ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
  * `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

<a id='NeuroAnalyzer.ica_remove' href='#NeuroAnalyzer.ica_remove'>#</a>
**`NeuroAnalyzer.ica_remove`** &mdash; *Function*.



```julia
ica_remove(obj, ic, ic_mw; ch, ic_idx)
```

Remove external ICA components from the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
ica_remove(obj, ic; ch, ic_idx)
```

Remove embedded ICA components (`:ic` and `:ic_mw`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.ica_remove!' href='#NeuroAnalyzer.ica_remove!'>#</a>
**`NeuroAnalyzer.ica_remove!`** &mdash; *Function*.



```julia
ica_remove!(obj, ic, ic_mw; ch, ic_idx)
```

Remove external ICA components from the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep


```
ica_remove!(obj; ch, ic_idx)
```

Remove embedded ICA components (`:ic` and `:ic_mw`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep

<a id='NeuroAnalyzer.icw_trans' href='#NeuroAnalyzer.icw_trans'>#</a>
**`NeuroAnalyzer.icw_trans`** &mdash; *Function*.



```julia
icw_trans(ct; wt, type)
```

Perform inverse continuous wavelet transformation (iCWT).

**Arguments**

  * `ct::AbstractArray`: CWT coefficients (by rows)
  * `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `type::Symbol=df`: inverse style type:

      * `:nd`: NaiveDelta
      * `:pd`: PenroseDelta
      * `:df`: DualFrames

**Returns**

  * `s::Vector{Float64}`: reconstructed signal

<a id='NeuroAnalyzer.idw_trans' href='#NeuroAnalyzer.idw_trans'>#</a>
**`NeuroAnalyzer.idw_trans`** &mdash; *Function*.



```julia
idw_trans(dwt_coefs; wt, type)
```

Perform inverse discrete wavelet transformation (iDWT) of the `dwt_coefs`.

**Arguments**

  * `dwt_coefs::AbstractArray`: DWT coefficients cAl, cD1, ..., cDl (by rows)
  * `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
  * `type::Symbol`: transformation type: 

      * `:sdwt`: Stationary Wavelet Transforms
      * `:acdwt`: Autocorrelation Wavelet Transforms

**Returns**

  * `s_new::Vector{Float64}`: reconstructed signal

<a id='NeuroAnalyzer.intensity2od' href='#NeuroAnalyzer.intensity2od'>#</a>
**`NeuroAnalyzer.intensity2od`** &mdash; *Function*.



```julia
intensity2od(s)
```

Convert NIRS intensity (RAW data) to optical density (OD).

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `od::AbstractArray`


```
intensity2od(obj; ch)
```

Convert NIRS intensity (RAW data) to optical density (OD).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))`: index of channels, default is NIRS intensity channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.intensity2od!' href='#NeuroAnalyzer.intensity2od!'>#</a>
**`NeuroAnalyzer.intensity2od!`** &mdash; *Function*.



```julia
intensity2od!(obj; ch)
```

Convert NIRS intensity (RAW data) to optical density (OD).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_int"))`: index of channels, default is NIRS intensity channels

<a id='NeuroAnalyzer.invert_polarity' href='#NeuroAnalyzer.invert_polarity'>#</a>
**`NeuroAnalyzer.invert_polarity`** &mdash; *Function*.



```julia
invert_polarity(obj; ch)
```

Invert polarity.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.invert_polarity!' href='#NeuroAnalyzer.invert_polarity!'>#</a>
**`NeuroAnalyzer.invert_polarity!`** &mdash; *Function*.



```julia
invert_polarity!(obj; ch)
```

Invert polarity.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to invert

<a id='NeuroAnalyzer.lrinterpolate_channel' href='#NeuroAnalyzer.lrinterpolate_channel'>#</a>
**`NeuroAnalyzer.lrinterpolate_channel`** &mdash; *Function*.



```julia
lrinterpolate_channel(obj; ch, ep)
```

Interpolate channel using linear regression.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number to interpolate
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.lrinterpolate_channel!' href='#NeuroAnalyzer.lrinterpolate_channel!'>#</a>
**`NeuroAnalyzer.lrinterpolate_channel!`** &mdash; *Function*.



```julia
lrinterpolate_channel!(obj; ch, ep)
```

Interpolate channel using linear regression.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number to interpolate
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.normalize' href='#NeuroAnalyzer.normalize'>#</a>
**`NeuroAnalyzer.normalize`** &mdash; *Function*.



```julia
normalize(s, n; method)
```

Normalize.

**Arguments**

  * `s::AbstractArray`
  * `m::Real=0.0`
  * `n::Real=1.0`
  * `method::Symbol`:

      * `:zscore`: by z-score
      * `:minmax`: in [-1, +1]
      * `:log`: using log-transformation
      * `:log10`: using log10-transformation
      * `:neglog`: using -log-transformation
      * `:neglog10`: using -log10-transformation
      * `:neg`: in [-∞, 0]
      * `:pos`: in [0, +∞]
      * `:perc`: in percentages
      * `:gauss`: to Gaussian
      * `:invroot`: in inverse root (1/sqrt(x))
      * `:softmax`: exp(x_i) / sum(exp(x))
      * `:n`: in [0, n], default is [0, 1]; to normalize to [n1, n2], use `normalize_n(s) .* (n2 - n1) .+ n1`
      * `:none`

**Returns**

  * `normalized::Vector{Float64}`


```
normalize(obj; ch, method)
```

Normalize channel(s)

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `method::Symbol`: method for normalization, see `normalize()` for details

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.normalize!' href='#NeuroAnalyzer.normalize!'>#</a>
**`NeuroAnalyzer.normalize!`** &mdash; *Function*.



```julia
normalize!(obj; ch, method)
```

Normalize channel(s)

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `method::Symbol`: method for normalization, see `normalize()` for details

<a id='NeuroAnalyzer.normalize_gauss' href='#NeuroAnalyzer.normalize_gauss'>#</a>
**`NeuroAnalyzer.normalize_gauss`** &mdash; *Function*.



```julia
normalize_gauss(s, dims)
```

Normalize to Gaussian.

**Arguments**

  * `s::AbstractArray`
  * `dims::Int64=1`: dimension for cumsum()

**Returns**

  * `normalize_gauss::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_invroot' href='#NeuroAnalyzer.normalize_invroot'>#</a>
**`NeuroAnalyzer.normalize_invroot`** &mdash; *Function*.



```julia
normalize_invroot(s)
```

Normalize in inverse root (1/sqrt(x)).

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_invroot::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_log' href='#NeuroAnalyzer.normalize_log'>#</a>
**`NeuroAnalyzer.normalize_log`** &mdash; *Function*.



```julia
normalize_log(s)
```

Normalize using log-transformation.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_log::AbstractArray`

<a id='NeuroAnalyzer.normalize_log10' href='#NeuroAnalyzer.normalize_log10'>#</a>
**`NeuroAnalyzer.normalize_log10`** &mdash; *Function*.



```julia
normalize_log10(s)
```

Normalize using log10-transformation.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_log10::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_minmax' href='#NeuroAnalyzer.normalize_minmax'>#</a>
**`NeuroAnalyzer.normalize_minmax`** &mdash; *Function*.



```julia
normalize_minmax(s)
```

Normalize in [-1, +1].

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_minmax::AbstractArray`

<a id='NeuroAnalyzer.normalize_n' href='#NeuroAnalyzer.normalize_n'>#</a>
**`NeuroAnalyzer.normalize_n`** &mdash; *Function*.



```julia
normalize_n(s, n)
```

Normalize in [0, n], default is [0, +1].

**Arguments**

  * `s::AbstractArray`
  * `n::Real=1.0`

**Returns**

  * `normalize_n::AbstractArray`

<a id='NeuroAnalyzer.normalize_neg' href='#NeuroAnalyzer.normalize_neg'>#</a>
**`NeuroAnalyzer.normalize_neg`** &mdash; *Function*.



```julia
normalize_neg(s)
```

Normalize in [-∞, 0].

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_neg::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_neglog' href='#NeuroAnalyzer.normalize_neglog'>#</a>
**`NeuroAnalyzer.normalize_neglog`** &mdash; *Function*.



```julia
normalize_neglog(s)
```

Normalize to using -log-transformation.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_neglog::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_neglog10' href='#NeuroAnalyzer.normalize_neglog10'>#</a>
**`NeuroAnalyzer.normalize_neglog10`** &mdash; *Function*.



```julia
normalize_neglog10(s)
```

Normalize using -log10-transformation.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_neglog::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_perc' href='#NeuroAnalyzer.normalize_perc'>#</a>
**`NeuroAnalyzer.normalize_perc`** &mdash; *Function*.



```julia
normalize_perc(s)
```

Normalize in percentages.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_perc::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_pos' href='#NeuroAnalyzer.normalize_pos'>#</a>
**`NeuroAnalyzer.normalize_pos`** &mdash; *Function*.



```julia
normalize_pos(s)
```

Normalize in [0, +∞].

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_pos::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_softmax' href='#NeuroAnalyzer.normalize_softmax'>#</a>
**`NeuroAnalyzer.normalize_softmax`** &mdash; *Function*.



```julia
normalize_softmax(s)
```

Softmax normalize: `exp(x_i) / sum(exp(x))`

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_softmax::Vector{Float64}`

<a id='NeuroAnalyzer.normalize_zscore' href='#NeuroAnalyzer.normalize_zscore'>#</a>
**`NeuroAnalyzer.normalize_zscore`** &mdash; *Function*.



```julia
normalize_zscore(s)
```

Normalize (by z-score).

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `normalize_zscore::Vector{Float64}`

<a id='NeuroAnalyzer.normpower' href='#NeuroAnalyzer.normpower'>#</a>
**`NeuroAnalyzer.normpower`** &mdash; *Function*.



```julia
normpower(s)
```

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `s_new::Vector{Float64}`


```
normpower(s)
```

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::Array{Float64, 3}`


```
normpower(obj; channel, t)
```

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.normpower!' href='#NeuroAnalyzer.normpower!'>#</a>
**`NeuroAnalyzer.normpower!`** &mdash; *Function*.



```julia
normpower!(obj; ch)
```

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

<a id='NeuroAnalyzer.npl' href='#NeuroAnalyzer.npl'>#</a>
**`NeuroAnalyzer.npl`** &mdash; *Function*.



```julia
npl(obj)
```

Calculate non-phase-locked signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: must be ERP object

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.npl!' href='#NeuroAnalyzer.npl!'>#</a>
**`NeuroAnalyzer.npl!`** &mdash; *Function*.



```julia
npl!(obj)
```

Calculate non-phase-locked signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: must be ERP object

<a id='NeuroAnalyzer.od2conc' href='#NeuroAnalyzer.od2conc'>#</a>
**`NeuroAnalyzer.od2conc`** &mdash; *Function*.



```julia
od2conc(obj; ch, ppf)
```

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"))`: index of channels, default is NIRS intensity channels
  * `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: Partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.od2conc!' href='#NeuroAnalyzer.od2conc!'>#</a>
**`NeuroAnalyzer.od2conc!`** &mdash; *Function*.



```julia
od2conc(obj; ch, ppf)
```

Convert NIRS optical density (OD) to concentration (HbO, HbR, HbT).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="nirs_od"))`: index of channels, default is NIRS intensity channels
  * `ppf::Vector{Real}=ones(length(obj.header.recording[:wavelengths]))`: Partial path length factors for each wavelength. This is a vector of factors per wavelength. Typical value is ~6 for each wavelength if the absorption change is uniform over the volume of tissue measured. To approximate the partial volume effect of a small localized absorption change within an adult human head, this value could be as small as 0.1. Convention is becoming to set `ppf=1` and to not divide by the source-detector separation such that the resultant "concentration" is in units of Molar mm (or Molar cm if those are the spatial units). This is becoming wide spread in the literature but there is no fixed citation. Use a value of 1 to choose this option.

<a id='NeuroAnalyzer.pca_decompose' href='#NeuroAnalyzer.pca_decompose'>#</a>
**`NeuroAnalyzer.pca_decompose`** &mdash; *Function*.



```julia
pca_decompose(s, n)
```

Calculate `n` first Primary Components (PCs).

**Arguments**

  * `s::AbstractArray`
  * `n::Int64`: number of PCs

**Returns**

Named tuple containing:

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pcv::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch (% of total variance)
  * `pcm::PCA{Float64}`: PC mean
  * `pc_model::MultivariateStats.PCA{Float64}`: PC model


```
pca_decompose(obj; ch, n)
```

Calculate `n` first Primary Components (PCs).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `n::Int64`: number of PCs to calculate

**Returns**

Named tuple containing:

  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pcv::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch (% of total variance)
  * `pcm::Vector{Float64}`: PC mean
  * `pc_model::MultivariateStats.PCA{Float64}`: PC model

<a id='NeuroAnalyzer.pca_reconstruct' href='#NeuroAnalyzer.pca_reconstruct'>#</a>
**`NeuroAnalyzer.pca_reconstruct`** &mdash; *Function*.



```julia
pca_reconstruct(s, pc, pca)
```

Reconstructs signal using PCA components.

**Arguments**

  * `s::AbstractArray`
  * `pc::AbstractArray:`: IC(1)..IC(n) × epoch
  * `pc_model::MultivariateStats.PCA{Float64}:`: PC model

**Returns**

  * `s_new::Array{Float64, 3}`


```
pca_reconstruct(obj; ch)
```

Reconstruct signal using embedded PCA components (`:pc`) and model (`:pc_model`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
pca_reconstruct(obj, pc, pc_model; ch)
```

Reconstruct signal using external PCA components (`pc` and `pca`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_model::MultivariateStats.PCA{Float64}`: PC model
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.pca_reconstruct!' href='#NeuroAnalyzer.pca_reconstruct!'>#</a>
**`NeuroAnalyzer.pca_reconstruct!`** &mdash; *Function*.



```julia
pca_reconstruct!(obj; ch)
```

Reconstruct signal using embedded PCA components (`:pc`) and model (`:pc_model`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels


```
pca_reconstruct!(obj, pc, pc_model; ch)
```

Reconstruct signals using external PCA components (`pc` and `pc_model`).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
  * `pc_model::MultivariateStats.PCA{Float64}`: PC model
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

<a id='NeuroAnalyzer.plinterpolate' href='#NeuroAnalyzer.plinterpolate'>#</a>
**`NeuroAnalyzer.plinterpolate`** &mdash; *Function*.



```julia
plinterpolate(s; locs, ch, imethod, nmethod, cart)
```

Interpolate channel(s) using planar interpolation.

**Arguments**

  * `s::Vector{<:Real}`: values to plot (one value per channel)
  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `interpolation_factor::Int64=100`: interpolation quality

**Returns**

  * `int_s::Matrix{Float64}`: interpolated signal
  * `int_x::Vector{Float64}`: X-axis coordinates
  * `int_y::Vector{Float64}`: Y-axis coordinates

<a id='NeuroAnalyzer.plinterpolate_channel' href='#NeuroAnalyzer.plinterpolate_channel'>#</a>
**`NeuroAnalyzer.plinterpolate_channel`** &mdash; *Function*.



```julia
plinterpolate_channel(obj; ch, ep, m, q)
```

Interpolate channel(s) using planar interpolation.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to interpolate
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `interpolation_factor::Int64=100`: interpolation quality

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.plinterpolate_channel!' href='#NeuroAnalyzer.plinterpolate_channel!'>#</a>
**`NeuroAnalyzer.plinterpolate_channel!`** &mdash; *Function*.



```julia
plinterpolate_channel!(obj; ch, ep, imethod, interpolation_factor)
```

Interpolate channel(s) using planar interpolation.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to interpolate
  * `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
  * `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
  * `interpolation_factor::Int64=100`: interpolation quality

<a id='NeuroAnalyzer.reference_a' href='#NeuroAnalyzer.reference_a'>#</a>
**`NeuroAnalyzer.reference_a`** &mdash; *Function*.



```julia
reference_a(obj; type, med)
```

Reference to auricular (A1, A2) channels. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::Symbol=:l`:

      * `:l`: linked - average of A1 and A2
      * `:i`: ipsilateral - A1 for left channels, A2 for right channels
      * `:c`: contraletral - A1 for right channels, A2 for left channels
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reference_a!' href='#NeuroAnalyzer.reference_a!'>#</a>
**`NeuroAnalyzer.reference_a!`** &mdash; *Function*.



```julia
reference_a!(obj; type, med)
```

Reference to auricular (A1, A2) channels. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::Symbol=:l`:

      * `:l`: linked - average of A1 and A2
      * `:i`: ipsilateral - A1 for left channels, A2 for right channels
      * `:c`: contraletral - A1 for right channels, A2 for left channels
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.reference_avg' href='#NeuroAnalyzer.reference_avg'>#</a>
**`NeuroAnalyzer.reference_avg`** &mdash; *Function*.



```julia
reference_avg(obj; exclude_fpo, exclude_current, average, med, weighted)
```

Reference to averaged reference. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from CAR calculation
  * `exclude_current::Bool=false`: exclude current channel from CAR calculation
  * `average::Bool=true`: average reference channels prior to subtracting, otherwise add all reference channels
  * `med::Bool=false`: use median instead of mean
  * `weighted::Bool=false`: use weighted reference channels (weights depend on the distance from the current electrode)

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reference_avg!' href='#NeuroAnalyzer.reference_avg!'>#</a>
**`NeuroAnalyzer.reference_avg!`** &mdash; *Function*.



```julia
reference_avg!(obj; exclude_fpo, exclude_current, average, med, weighted)
```

Reference to averaged reference. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `exclude_fpo::Bool=false`: exclude Fp1, Fp2 (due to eye blinks), O1, O2 (due to head movements) from CAR calculation
  * `exclude_current::Bool=false`: exclude current channel from CAR mean calculation
  * `average::Bool=true`: average reference channels prior to subtracting, otherwise add all reference channels
  * `med::Bool=false`: use median instead of mean
  * `weighted::Bool=false`: use weighted reference channels (weights depend on the distance from the current electrode)

<a id='NeuroAnalyzer.reference_ce' href='#NeuroAnalyzer.reference_ce'>#</a>
**`NeuroAnalyzer.reference_ce`** &mdash; *Function*.



```julia
reference_ce(obj; ch, med)
```

Reference to common electrode(s). Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reference_ce!' href='#NeuroAnalyzer.reference_ce!'>#</a>
**`NeuroAnalyzer.reference_ce!`** &mdash; *Function*.



```julia
reference_ce!(obj; ch, med)
```

Reference to common electrode(s). Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.reference_custom' href='#NeuroAnalyzer.reference_custom'>#</a>
**`NeuroAnalyzer.reference_custom`** &mdash; *Function*.



```julia
reference_custom(obj; ref_list, ref_name)
```

Reference using custom montage. Only signal channels are processed. Custom montage may be imported using `import_montage()`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
  * `ref_name::String="longitudinal-BIP"`: name of the montage

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

**Notes**

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data. For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 2 ("Cz") - amplitude of channel 1 ("Fz").

Examples of montages:

  * bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "transverse-BIP"
  * bipolar longitudinal: ["Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
  * bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"
  * bipolar longitudinal: ["FPz-Fz", "Fz-Cz", "Cz-Pz", "Pz-Oz", "Fp1-F7", "Fp1-F3", "F7-T7", "F7-T3", "T7-P7", "T7-T5", "T3-P7", "T3-T5", "P7-O1", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "F8-T8", "F8-T4", "T8-P8", "T8-T6", "T4-P8", "T4-T6", "T6-O2", "P8-O2", "Fp2-F4", "F4-C4", "C4-P4", "P4-O2"], "longitudinal-BIP"

<a id='NeuroAnalyzer.reference_custom!' href='#NeuroAnalyzer.reference_custom!'>#</a>
**`NeuroAnalyzer.reference_custom!`** &mdash; *Function*.



```julia
reference_custom!(obj; ref_list, ref_name)
```

Reference using custom montage. Only signal channels are processed. Custom montage may be imported using `import_montage()`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ref_list::Vector{String}=["Fz-Cz", "Cz-Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"]`: list of channel pairs
  * `ref_name::String="BIP ||"`: name of the montage

**Notes**

If the reference contains a single channel (e.g. "Fz"), than the channel is copied to the referenced data. For each reference pair (e.g. "Fz-Cz"), the referenced channel is equal to the amplitude of channel 2 ("Cz") - amplitude of channel 1 ("Fz").

Examples of montages:

  * bipolar transverse: ["Fp2-Fp1", "F8-Fp2", "F8-F4", "F4-Fz", "Fz-F3", "F3-F7", "Fp1-F7", "T4-C4", "C4-Cz", "Cz-C3", "C3-T3", "T6-P4", "P4-Pz", "Pz-P3", "P3-T5", "O2-O1"], "BIP ="
  * bipolar longitudinal: ["Fz", "Cz", "Pz", "Fp1-F7", "Fp1-F3", "F7-T3", "T3-T5", "T5-O1", "F3-C3", "C3-P3", "P3-O1", "Fp2-F8", "Fp2-F4", "F8-T4", "T4-T6", "T6-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"
  * bipolar longitudinal: ["Fp-Fz", "Fz-Cz", "Cz-Pz", "Pz-O", "Fp1-F7", "Fp1-F3", "F7-T7", "T7-P7", "P7-O1", "F3-C3", "C3-P3", "P3-O1", "Fp1-F7", "Fp2-F4", "F8-T8", "T8-P8", "P8-O2", "F4-C4", "C4-P4", "P4-O2"], "BIP ||"

<a id='NeuroAnalyzer.reference_m' href='#NeuroAnalyzer.reference_m'>#</a>
**`NeuroAnalyzer.reference_m`** &mdash; *Function*.



```julia
reference_m(obj; type, med)
```

Reference to mastoid (M1, M2) channels. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::Symbol=:l`:

      * `:l`: linked - average of M1 and M2
      * `:i`: ipsilateral - M1 for left channels, M2 for right channels
      * `:c`: contraletral - M1 for right channels, M2 for left channels
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reference_m!' href='#NeuroAnalyzer.reference_m!'>#</a>
**`NeuroAnalyzer.reference_m!`** &mdash; *Function*.



```julia
reference_m!(obj; type, med)
```

Reference to mastoid (M1, M2) channels. Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `type::Symbol=:l`:

      * `:l`: linked - average of M1 and M2
      * `:i`: ipsilateral - M1 for left channels, M2 for right channels
      * `:c`: contraletral - M1 for right channels, M2 for left channels
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.reference_plap' href='#NeuroAnalyzer.reference_plap'>#</a>
**`NeuroAnalyzer.reference_plap`** &mdash; *Function*.



```julia
reference_plap(obj; nn, weights)
```

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `nn::Int64=4`: use `nn` adjacent electrodes
  * `weighted::Bool=false`: use mean of `nn` nearest channels if false; if true, mean of `nn` nearest channels is weighted by distance to the referenced channel
  * `med::Bool=false`: use median instead of mean

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.reference_plap!' href='#NeuroAnalyzer.reference_plap!'>#</a>
**`NeuroAnalyzer.reference_plap!`** &mdash; *Function*.



```julia
reference_plap!(obj; nn, weights)
```

Reference using planar Laplacian (using `nn` adjacent electrodes). Only signal channels are processed.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `nn::Int64=4`: use `nn` adjacent electrodes
  * `weighted::Bool=false`: use distance weights; use mean of nearest channels if false
  * `med::Bool=false`: use median instead of mean

<a id='NeuroAnalyzer.remove_dc' href='#NeuroAnalyzer.remove_dc'>#</a>
**`NeuroAnalyzer.remove_dc`** &mdash; *Function*.



```julia
remove_dc(s, n)
```

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

**Arguments**

  * `s::AbstractVector`
  * `n::Int64=0`: baseline is the first `n` samples

**Returns**

  * `s_new::Vector{Float64}`


```
remove_dc(s; n)
```

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

**Arguments**

  * `s::AbstractArray`
  * `n::Int64=0`: baseline is the first `n` samples

**Returns**

  * `s::Array{Float64, 3}`


```
remove_dc(obj; ch, n)
```

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `n::Int64=0`: baseline is the first `n` samples

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.remove_dc!' href='#NeuroAnalyzer.remove_dc!'>#</a>
**`NeuroAnalyzer.remove_dc!`** &mdash; *Function*.



```julia
remove_dc!(obj; ch, n)
```

Remove mean value (DC offset). If `n` is greater than 0, mean value is calculated for the first `n` samples.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `n::Int64=0`: baseline is the first `n` samples

<a id='NeuroAnalyzer.remove_pops' href='#NeuroAnalyzer.remove_pops'>#</a>
**`NeuroAnalyzer.remove_pops`** &mdash; *Function*.



```julia
remove_pops(s; r, repair)
```

Detect and repair electrode pops (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≥2 seconds.

**Arguments**

  * `s::AbstractVector`
  * `r::Int64=20`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples
  * `repair::Bool=true`: recover the segment if `true`

**Returns**

Named tuple containing:

  * `s::Vector{Float64}`
  * `pop_location::Int64`: sample number in the signal
  * `left_segment::Int64`: length of segment before the pop that starts when signal crosses 0
  * `right_segment::Int64`: length of segment after the pop that ends when signal crosses 0


```
remove_pops(obj; <keyword arguments>)
```

Detect and repair electrode pops (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected per segment, signal length should be ≈2 seconds.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `repair::Bool=true`: recover the segment if `true`
  * `window::Real=10.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
  * `r::Int64=sr(obj)÷2`: detection segment length; pops are checked within `(pop_location - r):(pop_location + r)` samples

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`: returned if `repair=true`
  * `pop_loc::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
  * `l_seg::Vector{Int64}`: length of segment before the pop that starts when signal crosses 0
  * `r_seg::Vector{Int64}`: length of segment after the pop that ends when signal crosses 0

<a id='NeuroAnalyzer.remove_pops!' href='#NeuroAnalyzer.remove_pops!'>#</a>
**`NeuroAnalyzer.remove_pops!`** &mdash; *Function*.



```julia
remove_pops!(obj; <keyword arguments>)
```

Detect and repair electrode pops (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≈2 seconds.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `repair::Bool=true`: recover the segment if `true`
  * `window::Real=20.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
  * `r::Int64=sr(obj)÷2`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples

**Returns**

  * `pop_loc::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
  * `l_seg::Vector{Int64}`: length of segment before the pop that starts when signal crosses 0
  * `r_seg::Vector{Int64}`: length of segment after the pop that ends when signal crosses 0

<a id='NeuroAnalyzer.remove_powerline' href='#NeuroAnalyzer.remove_powerline'>#</a>
**`NeuroAnalyzer.remove_powerline`** &mdash; *Function*.



```julia
remove_powerline(obj; <keyword arguments>)
```

Remove power line noise and its peaks above power line frequency.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pl_frq::Real=50`: power line frequency
  * `method::Symbol=:iir`: use IIR filter
  * `pr::Real=2.0`: prominence of noise peaks in dB
  * `d::Real=5.0`: minimum distance between peaks in Hz
  * `q::Real=0.1`: optimization step size

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`
  * `df::DataFrame`: list of peaks detected

<a id='NeuroAnalyzer.remove_powerline!' href='#NeuroAnalyzer.remove_powerline!'>#</a>
**`NeuroAnalyzer.remove_powerline!`** &mdash; *Function*.



```julia
remove_powerline!(obj; <keyword arguments>)
```

Remove power line noise and harmonics.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pl_frq::Real=50`: power line frequency
  * `method::Symbol=:iir`: use IIR filter
  * `pr::Real=2.0`: prominence of noise peaks in dB
  * `d::Real=5.0`: minimum distance between peaks in Hz
  * `q::Real=0.1`: optimization step size

**Returns**

  * `df::DataFrame`: list of peaks detected

<a id='NeuroAnalyzer.resample' href='#NeuroAnalyzer.resample'>#</a>
**`NeuroAnalyzer.resample`** &mdash; *Function*.



```julia
resample(s; t, new_sr)
```

Resample to `new_sr` sampling frequency.

**Arguments**

  * `s::AbstractVector`
  * `old_sr::Int64`: old sampling rate
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_new::Vector{Float64}`


```
resample(s; old_sr::Int64, new_sr::Int64)
```

Resamples all channels and time vector `t` to `new_sr` sampling frequency.

**Arguments**

  * `s::AbstractArray`
  * `old_sr::Int64`: old sampling rate
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `s_new::Array{Float64, 3}`


```
resample(obj; new_sr)
```

Resample (up- or down-sample).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `old_sr::Int64`: old sampling rate - `new_sr::Int64`: new sampling rate

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.resample!' href='#NeuroAnalyzer.resample!'>#</a>
**`NeuroAnalyzer.resample!`** &mdash; *Function*.



```julia
resample!(obj; new_sr)
```

Resample (up- or down-sample).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroAnalyzer.scale' href='#NeuroAnalyzer.scale'>#</a>
**`NeuroAnalyzer.scale`** &mdash; *Function*.



```julia
scale(obj; ch, factor)
```

Multiply channel(s) by `factor`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `factor::Real`: signal is multiplied by `factor`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.scale!' href='#NeuroAnalyzer.scale!'>#</a>
**`NeuroAnalyzer.scale!`** &mdash; *Function*.



```julia
scale!(obj; ch, factor)
```

Multiply channel(s) by `factor`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `factor::Real`: signal is multiplied by `factor`

<a id='NeuroAnalyzer.sort_epochs' href='#NeuroAnalyzer.sort_epochs'>#</a>
**`NeuroAnalyzer.sort_epochs`** &mdash; *Function*.



```julia
sort_epochs(obj; s)
```

Sort epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `s::Vector{Int64}`: vector of sorting indices

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.sort_epochs!' href='#NeuroAnalyzer.sort_epochs!'>#</a>
**`NeuroAnalyzer.sort_epochs!`** &mdash; *Function*.



```julia
sort_epochs(obj; s)
```

Sort epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `s::Vector{Int64}`: vector of sorting indices

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.standardize' href='#NeuroAnalyzer.standardize'>#</a>
**`NeuroAnalyzer.standardize`** &mdash; *Function*.



```julia
standardize(s)
```

Standardize channels.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `s_new::NeuroAnalyzer.NEURO`:
  * `scaler::Vector{Any}`: standardizing matrix


```
standardize(obj; ch)
```

Standardize channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`
  * `scaler::Vector{Any}`: standardizing matrix

<a id='NeuroAnalyzer.standardize!' href='#NeuroAnalyzer.standardize!'>#</a>
**`NeuroAnalyzer.standardize!`** &mdash; *Function*.



```julia
standardize!(obj; ch)
```

Standardize channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `scaler::Matrix{Float64}`: standardizing matrix

<a id='NeuroAnalyzer.taper' href='#NeuroAnalyzer.taper'>#</a>
**`NeuroAnalyzer.taper`** &mdash; *Function*.



```julia
taper(s; taper)
```

Taper the signal.

**Arguments**

  * `s::AbstractVector`
  * `t::Union{AbstractVector, Vector{ComplexF64}}`

**Returns**

  * `s_new::Vector{Union{Float64, ComplexF64}}`


```
taper(s; t)
```

Taper the signal.

**Arguments**

  * `s::AbstractArray`
  * `t::Union{Vector{Real, Vector{ComplexF64}}`

**Returns**

  * `s_new::Array{Float64, 3`


```
taper(obj; channel, t)
```

Taper the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `t::Union{Vector{Real, Vector{ComplexF64}}`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.taper!' href='#NeuroAnalyzer.taper!'>#</a>
**`NeuroAnalyzer.taper!`** &mdash; *Function*.



```julia
taper!(obj; ch, t)
```

Taper the signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `t::Union{Vector{<:Real}, Vector{ComplexF64}}`

<a id='NeuroAnalyzer.tconv' href='#NeuroAnalyzer.tconv'>#</a>
**`NeuroAnalyzer.tconv`** &mdash; *Function*.



```julia
tconv(signal; kernel)
```

Performs convolution in the time domain.

**Arguments**

  * `s::AbstractVector`
  * `kernel::AbstractVector`

**Returns**

  * `s_new::Vector{Float64}`


```
tconv(s; kernel)
```

Perform convolution in the time domain.

**Arguments**

  * `s::AbstractArray`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

**Returns**

  * `s_new::Array{Float64, 3}`: convoluted signal


```
tconv(obj; ch, kernel)
```

Perform convolution in the time domain.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

**Returns**

  * `s_new::Array{Float64, 3}`: convoluted signal

<a id='NeuroAnalyzer.upsample' href='#NeuroAnalyzer.upsample'>#</a>
**`NeuroAnalyzer.upsample`** &mdash; *Function*.



```julia
upsample(obj; new_sr)
```

Upsample.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `new_sr::Int64`: new sampling rate

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.upsample!' href='#NeuroAnalyzer.upsample!'>#</a>
**`NeuroAnalyzer.upsample!`** &mdash; *Function*.



```julia
upsample!(obj; new_sr)
```

Upsample.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `new_sr::Int64`: new sampling rate

<a id='NeuroAnalyzer.wbp' href='#NeuroAnalyzer.wbp'>#</a>
**`NeuroAnalyzer.wbp`** &mdash; *Function*.



```julia
wbp(s; pad, frq, fs, ncyc)
```

Perform wavelet band-pass filtering.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64`: pad with `pad` zeros
  * `frq::Real`: filter frequency
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet

**Returns**

  * `s_new::Vector{Float64}`


```
wbp(s; ch, pad, frq, ncyc)
```

Perform wavelet band-pass filtering.

**Arguments**

  * `s::AbstractArray`
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `fs::Int64`: sampling rate
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet

**Returns**

  * `s_new::Array{Float64, 3`


```
wbp(obj; ch, pad, frq, ncyc)
```

Perform wavelet band-pass filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.wbp!' href='#NeuroAnalyzer.wbp!'>#</a>
**`NeuroAnalyzer.wbp!`** &mdash; *Function*.



```julia
wbp!(obj; ch, pad, frq, ncyc)
```

Perform wavelet band-pass filtering.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
  * `pad::Int64`: pad the `signal` with `pad` zeros
  * `frq::Real`: filter frequency
  * `ncyc::Int64=6`: number of cycles for Morlet wavelet


<a id='Locs'></a>

<a id='Locs-1'></a>

## Locs

<a id='NeuroAnalyzer.add_locs' href='#NeuroAnalyzer.add_locs'>#</a>
**`NeuroAnalyzer.add_locs`** &mdash; *Function*.



```julia
add_locs(obj; locs)
```

Add electrode positions from `locs`. 

Electrode locations:

  * `labels`          channel label
  * `loc_theta`       polar angle
  * `loc_radius`      polar radius
  * `loc_x`           Cartesian x
  * `loc_y`           Cartesian y
  * `loc_z`           Cartesian z
  * `loc_radius_sph`  spherical radius
  * `loc_theta_sph`   spherical horizontal angle
  * `loc_phi_sph`     spherical azimuth angle

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `locs::DataFrame`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`


```
add_locs(p1, p2; view, file_name)
```

Add locations to a plot. Locations are placed in the top right corner. If `file_name` is provided, the plot is saved as PNG file.

**Arguments**

  * `p1::Plots.Plot{Plots.GRBackend}`: signal plot
  * `p2::Plots.Plot{Plots.GRBackend}`: locations plot
  * `view::Bool=true`: view the output image
  * `file_name::String=""`: output image filename

**Returns**

  * `c::Cairo.CairoSurfaceBase{UInt32}`

<a id='NeuroAnalyzer.add_locs!' href='#NeuroAnalyzer.add_locs!'>#</a>
**`NeuroAnalyzer.add_locs!`** &mdash; *Function*.



```julia
add_locs!(obj; locs)
```

Load electrode positions from `locs` and return `NeuroAnalyzer.NEURO` object attached with channel locations data.

Electrode locations:

  * `labels`: channel label
  * `loc_theta`: polar angle
  * `loc_radius`: polar radius
  * `loc_x`: Cartesian X
  * `loc_y`: Cartesian Y
  * `loc_z`: Cartesian Z
  * `loc_radius_sph`: spherical radius
  * `loc_theta_sph`: spherical horizontal angle
  * `loc_phi_sph`: spherical azimuth angle

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `locs::DataFrame`

<a id='NeuroAnalyzer.cart2pol' href='#NeuroAnalyzer.cart2pol'>#</a>
**`NeuroAnalyzer.cart2pol`** &mdash; *Function*.



```julia
cart2pol(x, y)
```

Convert Cartesian coordinates to polar.

**Arguments**

  * `x::Real`
  * `y::Real`

**Returns**

  * `radius::Float64`
  * `theta::Float64`

<a id='NeuroAnalyzer.cart2sph' href='#NeuroAnalyzer.cart2sph'>#</a>
**`NeuroAnalyzer.cart2sph`** &mdash; *Function*.



```julia
cart2sph(x, y, z)
```

Convert spherical coordinates to Cartesian.

**Arguments**

  * `x::Real`
  * `y::Real`
  * `z::Real`

**Returns**

  * `radius::Float64`: spherical radius, the distance from the origin to the point
  * `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
  * `phi::Float64`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

<a id='NeuroAnalyzer.edit_locs' href='#NeuroAnalyzer.edit_locs'>#</a>
**`NeuroAnalyzer.edit_locs`** &mdash; *Function*.



```julia
edit_locs(obj; <keyword arguments>)
```

Edit electrode.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{String, Int64}`: channel number or name
  * `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
  * `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
  * `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
  * `theta::Union{Real, Nothing}=nothing`: polar angle
  * `radius::Union{Real, Nothing}=nothing`: polar radius
  * `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
  * `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
  * `name::String=""`: channel name
  * `type::String=""`: channel type

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.edit_locs!' href='#NeuroAnalyzer.edit_locs!'>#</a>
**`NeuroAnalyzer.edit_locs!`** &mdash; *Function*.



```julia
edit_locs!(obj; <keyword arguments>)
```

Edit electrode.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{String, Int64}`: channel number or name
  * `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
  * `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
  * `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
  * `theta::Union{Real, Nothing}=nothing`: polar angle
  * `radius::Union{Real, Nothing}=nothing`: polar radius
  * `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
  * `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
  * `name::String=""`: channel name
  * `type::String=""`: channel type

<a id='NeuroAnalyzer.locs_cart2pol' href='#NeuroAnalyzer.locs_cart2pol'>#</a>
**`NeuroAnalyzer.locs_cart2pol`** &mdash; *Function*.



```julia
locs_cart2pol(locs)
```

Convert Cartesian coordinates to polar.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_cart2pol!' href='#NeuroAnalyzer.locs_cart2pol!'>#</a>
**`NeuroAnalyzer.locs_cart2pol!`** &mdash; *Function*.



```julia
locs_cart2pol!(locs)
```

Convert Cartesian coordinates to polar.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_cart2sph' href='#NeuroAnalyzer.locs_cart2sph'>#</a>
**`NeuroAnalyzer.locs_cart2sph`** &mdash; *Function*.



```julia
locs_cart2sph(locs)
```

Convert Cartesian coordinates to spherical.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_cart2sph!' href='#NeuroAnalyzer.locs_cart2sph!'>#</a>
**`NeuroAnalyzer.locs_cart2sph!`** &mdash; *Function*.



```julia
locs_cart2sph!(locs)
```

Convert Cartesian coordinates to spherical.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_center' href='#NeuroAnalyzer.locs_center'>#</a>
**`NeuroAnalyzer.locs_center`** &mdash; *Function*.



```julia
locs_center(locs; polar, cart, spherical)
```

Center locs at (0, 0).

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.locs_center!' href='#NeuroAnalyzer.locs_center!'>#</a>
**`NeuroAnalyzer.locs_center!`** &mdash; *Function*.



```julia
locs_center!(locs; polar, cart, spherical)
```

Center locs at X=0.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_details' href='#NeuroAnalyzer.locs_details'>#</a>
**`NeuroAnalyzer.locs_details`** &mdash; *Function*.



```julia
locs_details(obj; ch, output)
```

Return locations of OBJ ch electrode.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, String}`: channel number or name
  * `out::Bool=true`: print output if true

**Returns**

Named tuple containing:

  * `ch::Int64`: channel number
  * `label::String`: location label
  * `theta::Float64`: polar angle
  * `radius::Float64`: polar radius
  * `x::Float64`: Cartesian X spherical coordinate
  * `y::Float64`: Cartesian Y spherical coordinate
  * `z::Float64`: Cartesian Z spherical coordinate
  * `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
  * `radius_sph::Float64`: spherical radius, the distance from the origin to the point
  * `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

<a id='NeuroAnalyzer.locs_flipx' href='#NeuroAnalyzer.locs_flipx'>#</a>
**`NeuroAnalyzer.locs_flipx`** &mdash; *Function*.



```julia
locs_flipx(locs; polar, cart, spherical)
```

Flip channel locations along x axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_flipx!' href='#NeuroAnalyzer.locs_flipx!'>#</a>
**`NeuroAnalyzer.locs_flipx!`** &mdash; *Function*.



```julia
locs_flipx!(locs; polar, cart, spherical)
```

Flip channel locations along x axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_flipy' href='#NeuroAnalyzer.locs_flipy'>#</a>
**`NeuroAnalyzer.locs_flipy`** &mdash; *Function*.



```julia
locs_flipy(locs; polar, cart, spherical)
```

Flip channel locations along y axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_flipy!' href='#NeuroAnalyzer.locs_flipy!'>#</a>
**`NeuroAnalyzer.locs_flipy!`** &mdash; *Function*.



```julia
locs_flipy!(locs; polar, cart, spherical)
```

Flip channel locations along y axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_flipz' href='#NeuroAnalyzer.locs_flipz'>#</a>
**`NeuroAnalyzer.locs_flipz`** &mdash; *Function*.



```julia
locs_flipz(locs; polar, cart, spherical)
```

Flip channel locations along z axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_flipz!' href='#NeuroAnalyzer.locs_flipz!'>#</a>
**`NeuroAnalyzer.locs_flipz!`** &mdash; *Function*.



```julia
locs_flipz!(locs; polar, cart, spherical)
```

Flip channel locations along z axis.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_generate' href='#NeuroAnalyzer.locs_generate'>#</a>
**`NeuroAnalyzer.locs_generate`** &mdash; *Function*.



```julia
locs_generate(locs)
```

Generate spherical coordinates according to 10/10 system.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`


```
locs_generate(obj)
```

Generate spherical coordinates according to 10/5 system.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `obj_new::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.locs_generate!' href='#NeuroAnalyzer.locs_generate!'>#</a>
**`NeuroAnalyzer.locs_generate!`** &mdash; *Function*.



```julia
locs_generate!(locs)
```

Generate spherical coordinates according to 10/5 system.

**Arguments**

  * `locs::DataFrame`


```
locs_generate!(obj)
```

Generate spherical coordinates according to 10/5 system.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.locs_normalize' href='#NeuroAnalyzer.locs_normalize'>#</a>
**`NeuroAnalyzer.locs_normalize`** &mdash; *Function*.



```julia
locs_normalize(locs; polar, cart, spherical)
```

Normalize channel locations to fit the unit sphere.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_normalize!' href='#NeuroAnalyzer.locs_normalize!'>#</a>
**`NeuroAnalyzer.locs_normalize!`** &mdash; *Function*.



```julia
locs_normalize!(locs; polar, cart, spherical)
```

Normalize channel locations to fit the unit sphere.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_pol2cart' href='#NeuroAnalyzer.locs_pol2cart'>#</a>
**`NeuroAnalyzer.locs_pol2cart`** &mdash; *Function*.



```julia
locs_pol2cart(locs)
```

Convert polar coordinates to Cartesian.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_pol2cart!' href='#NeuroAnalyzer.locs_pol2cart!'>#</a>
**`NeuroAnalyzer.locs_pol2cart!`** &mdash; *Function*.



```julia
locs_pol2cart!(locs)
```

Convert polar coordinates to Cartesian.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_pol2sph' href='#NeuroAnalyzer.locs_pol2sph'>#</a>
**`NeuroAnalyzer.locs_pol2sph`** &mdash; *Function*.



```julia
locs_pol2sph(locs)
```

Convert polar coordinates to spherical.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_pol2sph!' href='#NeuroAnalyzer.locs_pol2sph!'>#</a>
**`NeuroAnalyzer.locs_pol2sph!`** &mdash; *Function*.



```julia
locs_pol2sph!(locs)
```

Convert polar coordinates to spherical.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_rotx' href='#NeuroAnalyzer.locs_rotx'>#</a>
**`NeuroAnalyzer.locs_rotx`** &mdash; *Function*.



```julia
locs_rotx(locs; a, polar, cart, spherical)
```

Rotate channel locations around the X axis (in the YZ-plane).

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_rotx!' href='#NeuroAnalyzer.locs_rotx!'>#</a>
**`NeuroAnalyzer.locs_rotx!`** &mdash; *Function*.



```julia
locs_rotx!(locs; a, polar, cart, spherical)
```

Rotate channel locations around the X axis (in the YZ-plane).

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_roty' href='#NeuroAnalyzer.locs_roty'>#</a>
**`NeuroAnalyzer.locs_roty`** &mdash; *Function*.



```julia
locs_roty(locs; a, polar, cart, spherical)
```

Rotate channel locations around the Y axis (in the XZ-plane).

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_roty!' href='#NeuroAnalyzer.locs_roty!'>#</a>
**`NeuroAnalyzer.locs_roty!`** &mdash; *Function*.



```julia
locs_roty!(locs; a, polar, cart, spherical)
```

Rotate channel locations around the Y axis (in the XZ-plane).

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_rotz' href='#NeuroAnalyzer.locs_rotz'>#</a>
**`NeuroAnalyzer.locs_rotz`** &mdash; *Function*.



```julia
locs_rotz(locs; a, polar, cart, spherical)
```

Rotate channel locations around the Z axis.

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_rotz!' href='#NeuroAnalyzer.locs_rotz!'>#</a>
**`NeuroAnalyzer.locs_rotz!`** &mdash; *Function*.



```julia
locs_rotz!(locs; a, polar, cart, spherical)
```

Rotate channel locations in the xy-plane.

**Arguments**

  * `locs::DataFrame`
  * `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_scale' href='#NeuroAnalyzer.locs_scale'>#</a>
**`NeuroAnalyzer.locs_scale`** &mdash; *Function*.



```julia
locs_scale(locs; r, polar, cart, spherical)
```

Scale channel locations.

**Arguments**

  * `locs::DataFrame`
  * `r::Real`: scaling factor
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_scale!' href='#NeuroAnalyzer.locs_scale!'>#</a>
**`NeuroAnalyzer.locs_scale!`** &mdash; *Function*.



```julia
locs_scale!(locs, r, polar, cart, spherical)
```

Scale channel locations.

**Arguments**

  * `locs::DataFrame`
  * `r::Real`: scaling factor
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.locs_sph2cart' href='#NeuroAnalyzer.locs_sph2cart'>#</a>
**`NeuroAnalyzer.locs_sph2cart`** &mdash; *Function*.



```julia
locs_sph2cart(locs)
```

Convert spherical coordinates to Cartesian.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_sph2cart!' href='#NeuroAnalyzer.locs_sph2cart!'>#</a>
**`NeuroAnalyzer.locs_sph2cart!`** &mdash; *Function*.



```julia
locs_sph2cart!(locs)
```

Convert spherical coordinates to Cartesian.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_sph2pol' href='#NeuroAnalyzer.locs_sph2pol'>#</a>
**`NeuroAnalyzer.locs_sph2pol`** &mdash; *Function*.



```julia
locs_sph2pol(locs)
```

Convert spherical coordinates to polar.

**Arguments**

  * `locs::DataFrame`

**Returns**

  * `locs_new::DataFrame`

<a id='NeuroAnalyzer.locs_sph2pol!' href='#NeuroAnalyzer.locs_sph2pol!'>#</a>
**`NeuroAnalyzer.locs_sph2pol!`** &mdash; *Function*.



```julia
locs_sph2pol!(locs)
```

Convert Cartesian coordinates to polar.

**Arguments**

  * `locs::DataFrame`

<a id='NeuroAnalyzer.locs_swapxy' href='#NeuroAnalyzer.locs_swapxy'>#</a>
**`NeuroAnalyzer.locs_swapxy`** &mdash; *Function*.



```julia
locs_swapxy(locs; polar, cart, spherical)
```

Swap channel locations x and y axes.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

**Returns**

  * `obj::NeuroAnalyzer.NEURO`

<a id='NeuroAnalyzer.locs_swapxy!' href='#NeuroAnalyzer.locs_swapxy!'>#</a>
**`NeuroAnalyzer.locs_swapxy!`** &mdash; *Function*.



```julia
locs_swapxy!(locs; polar, cart, spherical)
```

Swap channel locations x and y axes.

**Arguments**

  * `locs::DataFrame`
  * `polar::Bool=true`: modify polar coordinates
  * `cart::Bool=true`: modify Cartesian coordinates
  * `spherical::Bool=true`: modify spherical coordinates

<a id='NeuroAnalyzer.pol2cart' href='#NeuroAnalyzer.pol2cart'>#</a>
**`NeuroAnalyzer.pol2cart`** &mdash; *Function*.



```julia
pol2cart(radius, theta)
```

Convert polar coordinates to Cartesian.

**Arguments**

  * `radius::Real`: polar radius, the distance from the origin to the point, in degrees
  * `theta::Real`: polar angle

**Returns**

  * `x::Float64`
  * `y::Float64`

<a id='NeuroAnalyzer.pol2sph' href='#NeuroAnalyzer.pol2sph'>#</a>
**`NeuroAnalyzer.pol2sph`** &mdash; *Function*.



```julia
pol2sph(radius, theta)
```

Convert polar coordinates to spherical.

**Arguments**

  * `radius::Real`: polar radius, the distance from the origin to the point, in degrees
  * `theta::Real`: polar angle

**Returns**

  * `radius::Float64`: spherical radius, the distance from the origin to the point
  * `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
  * `phi::Float64`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

<a id='NeuroAnalyzer.sph2cart' href='#NeuroAnalyzer.sph2cart'>#</a>
**`NeuroAnalyzer.sph2cart`** &mdash; *Function*.



```julia
sph2cart(radius, theta, phi)
```

Convert spherical coordinates to Cartesian.

**Arguments**

  * `radius::Real`: spherical radius, the distance from the origin to the point
  * `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
  * `phi::Real`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

**Returns**

  * `x::Float64`
  * `y::Float64`
  * `z::Float64`

<a id='NeuroAnalyzer.sph2pol' href='#NeuroAnalyzer.sph2pol'>#</a>
**`NeuroAnalyzer.sph2pol`** &mdash; *Function*.



```julia
sph2pol(radius, theta, phi)
```

Convert spherical coordinates to polar.

**Arguments**

  * `radius::Real`: spherical radius, the distance from the origin to the point
  * `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
  * `phi::Real`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

**Returns**

  * `radius::Real`: polar radius, the distance from the origin to the point
  * `theta::Real`: polar horizontal angle, the angle in the xy plane with respect to the x axis, in degrees


<a id='Analyze'></a>

<a id='Analyze-1'></a>

## Analyze


<a id='Analyze:-amplitude'></a>

<a id='Analyze:-amplitude-1'></a>

### Analyze: amplitude

<a id='NeuroAnalyzer.amp' href='#NeuroAnalyzer.amp'>#</a>
**`NeuroAnalyzer.amp`** &mdash; *Function*.



```julia
amp(s)
```

Calculate amplitudes.

**Arguments**

  * `s::AbstractVector`

**Returns**

Named tuple containing:

  * `p::Float64`: peak amplitude
  * `r::Float64`: RMS amplitude
  * `p2p::Float64`: peak-to-peak amplitude
  * `semi_p2p::Float64`: half of the peak-to-peak amplitude
  * `msa::Float64`: mean square amplitude
  * `rmsa::Float64`: root mean square amplitude
  * `energy::Float64`: total signal energy
  * `rms::Float64`: root mean square


```
amp(s)
```

Calculate amplitudes.

**Arguments**

  * `s::AbstractArray`

**Returns**

Named tuple containing:

  * `p::Matrix{Float64}`: peak amplitude
  * `r::Matrix{Float64}`: RMS amplitude
  * `p2p::Matrix{Float64}`: peak-to-peak amplitude
  * `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
  * `msa::Matrix{Float64}`: mean square amplitude
  * `rmsa::Matrix{Float64}`: root mean square amplitude
  * `energy::Matrix{Float64}`: total signal energy


```
amp(obj; ch)
```

Calculate amplitudes.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one

**Returns**

Named tuple containing:

  * `p::Matrix{Float64}`: peak amplitude
  * `r::Matrix{Float64}`: RMS amplitude
  * `p2p::Matrix{Float64}`: peak-to-peak amplitude
  * `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
  * `msa::Matrix{Float64}`: mean square amplitude
  * `rmsa::Matrix{Float64}`: root mean square amplitude
  * `energy::Matrix{Float64}`: total signal energy

<a id='NeuroAnalyzer.ampdiff' href='#NeuroAnalyzer.ampdiff'>#</a>
**`NeuroAnalyzer.ampdiff`** &mdash; *Function*.



```julia
ampdiff(s; ch)
```

Calculate amplitude difference between each channel and mean amplitude of reference channels.

**Arguments**

  * `s::AbstractArray`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=size(s, 1)`: index of reference channels, default is all channels except the analyzed one

**Returns**

  * `ad::Array{Float64, 3}`


```
ampdiff(obj; ch)
```

Calculate amplitude difference between each channel and mean amplitude of reference channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one

**Returns**

  * `ad::Array{Float64, 3}`

<a id='NeuroAnalyzer.rmse' href='#NeuroAnalyzer.rmse'>#</a>
**`NeuroAnalyzer.rmse`** &mdash; *Function*.



```julia
rmse(s1, s2)
```

Calculate Root Mean Square Error (RMSE).

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

  * `rmse::Float64`: RMSE


```
rmse(s1, s2)
```

Calculate Root Mean Square Error (RMSE).

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`

**Returns**

  * `r::Array{Float64, 2}`: RMSE


```
rmse(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate Root Mean Square Error (RMSE).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `r::Array{Float64, 3}`: RMSE
  * `cps_ph::Array{Float64, 3}`: cross power spectrum phase (in radians)
  * `cps_fq::Vector{Float64, 3}`: cross power spectrum frequencies

<a id='NeuroAnalyzer.snr' href='#NeuroAnalyzer.snr'>#</a>
**`NeuroAnalyzer.snr`** &mdash; *Function*.



```julia
snr(s)
```

Calculate mean-based SNR.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `snr::Float64`: SNR

**Source**

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278


```
snr(s; t, type)
```

Calculate SNR.

**Arguments**

  * `s::AbstractArray`
  * `t::Vector{Float64}`: epoch time
  * `type::Symbol=:rms`: SNR type:

      * `:mean`: mean-based
      * `:rms`: RMS-based

**Returns**

Named tuple containing:

  * `s::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
  * `f::Vector(Float64)`: frequencies


```
snr(obj; ch, type)
```

Calculate SNR.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `type::Symbol=:rms`: SNR type:

      * `:mean`: mean-based
      * `:rms`: RMS-based

**Returns**

Named tuple containing:

  * `sn::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
  * `f::Vector(Float64)`: frequencies

<a id='NeuroAnalyzer.snr2' href='#NeuroAnalyzer.snr2'>#</a>
**`NeuroAnalyzer.snr2`** &mdash; *Function*.



```julia
snr2(s)
```

Calculate RMS-based SNR.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `snr2::Float64`: SNR


<a id='Analyze:-frequency'></a>

<a id='Analyze:-frequency-1'></a>

### Analyze: frequency

<a id='NeuroAnalyzer.band_asymmetry' href='#NeuroAnalyzer.band_asymmetry'>#</a>
**`NeuroAnalyzer.band_asymmetry`** &mdash; *Function*.



```julia
band_asymmetry(obj; ch1, ch2, f, method, nt, wlen, woverlap)
```

Calculate band asymmetry: ln(channel 1 band power) - ln(channel 2 band power).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels, e.g. left frontal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels, e.g. right frontal channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `ba::Float64`: band asymmetry
  * `ba_norm::Float64`: normalized band asymmetry

<a id='NeuroAnalyzer.band_mpower' href='#NeuroAnalyzer.band_mpower'>#</a>
**`NeuroAnalyzer.band_mpower`** &mdash; *Function*.



```julia
band_mpower(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate mean and maximum band power and its frequency.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength(0, fs / 2)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `mbp::Float64`: mean band power
  * `maxfrq::Float64`: frequency of maximum band power
  * `maxbp::Float64`: power at maximum band frequency


```
band_mpower(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate absolute band power between two frequencies.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength(0, fs / 2)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `mbp::Matrix{Float64}`: mean band power per channel per epoch
  * `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
  * `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch


```
band_mpower(obj; ch, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate mean and maximum band power and its frequency.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `mbp::Matrix{Float64}`: mean band power per channel per epoch
  * `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
  * `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch

<a id='NeuroAnalyzer.band_power' href='#NeuroAnalyzer.band_power'>#</a>
**`NeuroAnalyzer.band_power`** &mdash; *Function*.



```julia
band_power(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate absolute band power between two frequencies.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `bp::Float64`: band power


```
band_power(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate absolute band power between two frequencies.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `bp::Matrix{Float64}`: band power


```
band_power(obj; ch, f, method, nt, wlen, woverlap)
```

Calculate absolute band power between two frequencies.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `bp::Matrix{Float64}`: band power

<a id='NeuroAnalyzer.coherence' href='#NeuroAnalyzer.coherence'>#</a>
**`NeuroAnalyzer.coherence`** &mdash; *Function*.



coherence(s1, s2; method, fs, frq_lim, demean, nt, wlen, woverlap, w, norm)

Calculate coherence and MSC (magnitude-squared coherence).

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

  * `coh::Vector{Float64}`: coherence
  * `mscoh::Vector{Float64}`: magnitude-squared coherence
  * `p::Vector{Float64}`: frequencies


coherence(s1, s2; method, fs, frq_lim, demean, nt, wlen, woverlap, w, norm)

Calculate coherence and MSC (magnitude-squared coherence).

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

  * `coh::Array{Float64, 3}`: coherence
  * `mscoh::Array{Float64, 3}`: magnitude-squared coherence
  * `p::Vector{Float64}`: frequencies


```
coherence(obj1, obj2; ch1, ch2, ep1, ep2, method, fs, frq_lim, demean, nt, wlen, woverlap, w, norm)
```

Calculate coherence and MSC (magnitude-squared coherence).

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

  * `coh::Array{Float64, 3}`: coherence
  * `mscoh::Array{Float64, 3}`: magnitude-squared coherence
  * `p::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.cpsd' href='#NeuroAnalyzer.cpsd'>#</a>
**`NeuroAnalyzer.cpsd`** &mdash; *Function*.



cpsd(s1, s2; method, fs, frq_lim, demean, nt, wlen, woverlap, w, norm)

Calculate cross power spectral density (CPSD).

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `norm::Bool=false`: normalize do dB

**Returns**

  * `pxy::Vector{Float64}`: cross-power spectrum
  * `p::Vector{Float64}`: frequencies


```
cpsd(s1, s2; method, fs, frq_lim, demean, nt, wlen, woverlap, w, norm)
```

Calculate cross power spectral density (CPSD).

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channelwise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `pxy::Array{Float64, 3}`: cross-power spectrum
  * `f::Vector{Float64}`: frequencies


```
cpsd(obj1, obj2; ch1, ch2, ep1, ep2, method, frq_lim, demean, nt, wlen, woverlap, w, norm)
```

Calculate cross power spectral density (CPSD).

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
  * `method::Symbol=:mt`: method used to calculate CPSD:

      * `:mt`: multi-tapered cross-power spectra
      * `:fft`: fast Fourier transformation
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj1) / 2)`: frequency bounds
  * `demean::Bool=false`: if true, the channelwise mean will be subtracted from the input signals before the cross spectral powers are computed
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj1)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `pxy::Array{Float64, 3}`: cross-power spectrum
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.cwtspectrogram' href='#NeuroAnalyzer.cwtspectrogram'>#</a>
**`NeuroAnalyzer.cwtspectrogram`** &mdash; *Function*.



```julia
cwtspectrogram(s; wt, pad, norm, frq_lim, fs)
```

Calculate spectrogram using continuous wavelet transformation (CWT).

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds for the spectrogram
  * `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `fs::Int64`: sampling rate
  * `w::Bool=true`: if true, apply Hanning window
  * `norm::Bool=true`: normalize powers to dB

**Returns**

Named tuple containing:

  * `sp::Matrix{Float64}`: powers
  * `sf::Vector{Float64}`: frequency indices

<a id='NeuroAnalyzer.erop' href='#NeuroAnalyzer.erop'>#</a>
**`NeuroAnalyzer.erop`** &mdash; *Function*.



```julia
erop(obj; <keyword arguments>)
```

Calculate ERO (Event-Related Oscillations) power-spectrum. If `obj` is ERP, `ero()` returns two epochs: ERP power-spectrum (`ero_s[:, :, 1]`) and averaged power-spectra of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged power-spectra of all `obj` epochs (`ero_s[:, :, 1]`)

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`: channel to analyze
  * `method::Symbol=:welch`: method of calculating power-spectrum:

      * `:welch`: Welch periodogram
      * `:stft`: short time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:fft`: Fast Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `ero_p::Array{Float64, 3}`: powers
  * `ero_f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.eros' href='#NeuroAnalyzer.eros'>#</a>
**`NeuroAnalyzer.eros`** &mdash; *Function*.



```julia
eros(obj; <keyword arguments>)
```

Calculate ERO (Event-Related Oscillations) spectrogram. If `obj` is ERP, `ero()` returns two epochs: ERP spectrogram (`ero_s[:, :, 1]`) and averaged spectrograms of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged spectrograms of all `obj` epochs (`ero_s[:, :, 1]`)

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`: channel to analyze
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length, default is 4 seconds
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `pad::Int64=0`: number of zeros to add
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets

**Returns**

Named tuple containing:

  * `ero_s::Array{Float64, 3}`: spectrogram(s)
  * `ero_f::Vector{Float64}`: frequencies
  * `ero_t::Vector{Float64}`: time

<a id='NeuroAnalyzer.frqinst' href='#NeuroAnalyzer.frqinst'>#</a>
**`NeuroAnalyzer.frqinst`** &mdash; *Function*.



```julia
frqinst(s; fs)
```

Calculate instantaneous frequency.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `f::Vector{Float64}`


```
frqinst(s; fs)
```

Calculate instantaneous frequency.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `f::Array{Float64, 2}`


```
frqinst(obj; channel)
```

Calculate instantaneous frequency.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `f::Array{Float64, 3}`

<a id='NeuroAnalyzer.ghspectrogram' href='#NeuroAnalyzer.ghspectrogram'>#</a>
**`NeuroAnalyzer.ghspectrogram`** &mdash; *Function*.



```julia
ghspectrogram(s; fs, norm, frq_lim, frq_n, frq, gw)
```

Calculate spectrogram using Gaussian and Hilbert transform.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `norm::Bool=true`: normalize powers to dB
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds for the spectrogram
  * `frq_n::Int64_tlength(frq_lim)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `sp::Matrix{Float64}`: powers
  * `sph::Matrix{Float64}`: phases
  * `sf::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.hjorth' href='#NeuroAnalyzer.hjorth'>#</a>
**`NeuroAnalyzer.hjorth`** &mdash; *Function*.



hjorth(s)

Calculate Hjorths parameters.

**Arguments**

  * `s::AbstractVector`

**Returns**

Named tuple containing:

  * `h_act::Float64`: activity
  * `h_mob::Float64`: mobility
  * `h_comp::Float64`: complexity

**Notes:**

  * Activity: the total power of the signal
  * Mobility: an estimate of the mean frequency
  * Complexity: indicates the similarity of the shape of the signal to a pure sine wave


hjorth(s)

Calculate Hjorths parameters.

**Arguments**

  * `s::AbstractArray`

**Returns**

Named tuple containing:

  * `h_act::Matrix{Float64}`: activity
  * `h_mob::Matrix{Float64}`: mobility
  * `h_comp::Matrix{Float64}`: complexity

**Notes:**

  * Activity: the total power of the signal
  * Mobility: an estimate of the mean frequency
  * Complexity: indicates the similarity of the shape of the signal to a pure sine wave


hjorth(obj)

Calculate Hjorths parameters.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `h_act::Matrix{Float64}`: activity
  * `h_mob::Matrix{Float64}`: mobility
  * `h_comp::Matrix{Float64}`: complexity

**Notes:**

  * Activity: the total power of the signal
  * Mobility: an estimate of the mean frequency
  * Complexity: indicates the similarity of the shape of the signal to a pure sine wave

<a id='NeuroAnalyzer.hspectrum' href='#NeuroAnalyzer.hspectrum'>#</a>
**`NeuroAnalyzer.hspectrum`** &mdash; *Function*.



```julia
hspectrum(s; pad=0)
```

Calculate amplitudes, powers and phases using Hilbert transform.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64`: number of zeros to add at the end of the signal
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `hc::Vector(ComplexF64}`: Hilbert components
  * `sa::Vector{Float64}`: amplitudes
  * `sp::Vector{Float64}`: powers
  * `sph::Vector{Float64}`: phases


```
hspectrum(s; pad, norm)
```

Calculate amplitudes, powers and phases using Hilbert transform.

**Arguments**

  * `s::AbstractArray`
  * `pad::Int64`: number of zeros to add at the end of the signal
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `hc::Array(ComplexF64, 3}`: Hilbert components
  * `sa::Array{Float64, 3}`: amplitudes
  * `sp::Array{Float64, 3}`: powers
  * `sph::Array{Float64, 3}`: phases

<a id='NeuroAnalyzer.mwpsd' href='#NeuroAnalyzer.mwpsd'>#</a>
**`NeuroAnalyzer.mwpsd`** &mdash; *Function*.



```julia
mwpsd(s; pad, norm, frq_n, frq, fs, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64=0`: pad with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `fs::Int64`: sampling rate
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `pw::Matrix{Float64}`: powers
  * `pf::Vector{Float64}`: frequencies


```
mwpsd(s; pad, norm, frq_n, frq, fs, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `s::AbstractMatrix`
  * `pad::Int64=0`: pad with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `fs::Int64`: sampling rate
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 2}`: powers
  * `pf::Vector{Float64}`: frequencies


```
mwpsd(s; pad, norm, frq_n, frq, fs, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `s::AbstractArray`
  * `pad::Int64=0`: pad with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `fs::Int64`: sampling rate
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies


```
mwpsd(obj; ch, pad, norm, frq_n, frq, ncyc)
```

Calculate power spectrum using Morlet wavelet convolution.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all s channels
  * `pad::Int64=0`: pad with `pad` zeros
  * `norm::Bool`=true: normalize powers to dB
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.mwspectrogram' href='#NeuroAnalyzer.mwspectrogram'>#</a>
**`NeuroAnalyzer.mwspectrogram`** &mdash; *Function*.



```julia
mwspectrogram(s; pad, norm, fs, frq_lim, frq_n, frq, ncyc)
```

Calculate spectrogram using wavelet convolution.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64`: pad with `pad` zeros
  * `norm::Bool=true`: normalize powers to dB
  * `fs::Int64`: sampling rate
  * `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds for the spectrogram
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `cs::Matrix(ComplexF64}`: convoluted signal
  * `sp::Matrix{Float64}`: powers
  * `sph::Matrix{Float64}`: phases
  * `sf::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.peak_frq' href='#NeuroAnalyzer.peak_frq'>#</a>
**`NeuroAnalyzer.peak_frq`** &mdash; *Function*.



```julia
peak_frq(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate peak frequency in a band.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `pf::Float64`: peak frequency


```
peak_frq(s; fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate peak frequency in a band.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `pf::Matrix{Float64}`: peak frequency


```
peak_frq(obj; ch, f, method, nt, wlen, woverlap)
```

Calculate peak frequency in a band.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `f::Tuple{Real, Real}`: lower and upper frequency bounds
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `pf::Matrix{Float64}`: peak frequency

<a id='NeuroAnalyzer.psd' href='#NeuroAnalyzer.psd'>#</a>
**`NeuroAnalyzer.psd`** &mdash; *Function*.



```julia
psd(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::Vector{Float64}`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Vector{Float64}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractMatrix`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 2}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd(obj; ch, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.psd_rel' href='#NeuroAnalyzer.psd_rel'>#</a>
**`NeuroAnalyzer.psd_rel`** &mdash; *Function*.



```julia
psd_rel(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate relative power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Vector{Float64}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd_rel(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate relative power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractMatrix`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd_rel(s; fs, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate relative power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `norm::Bool=false`: normalize do dB
  * `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Vector{Float64}`: frequencies


```
psd_rel(obj; ch, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate relative power spectrum density. Default method is Welch periodogram.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `norm::Bool=false`: normalize do dB
  * `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `pw::Array{Float64, 3}`: powers
  * `pf::Array{Float64, 3}`: frequencies

<a id='NeuroAnalyzer.psd_slope' href='#NeuroAnalyzer.psd_slope'>#</a>
**`NeuroAnalyzer.psd_slope`** &mdash; *Function*.



```julia
psd_slope(s; fs, f, norm, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate PSD linear fit and slope. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `lf::Vector{Float64}`: linear fit
  * `ls::Float64`: slopes of linear fit
  * `pf::Vector{Float64}`: range of frequencies for the linear fit


```
psd_slope(s; fs, f, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate PSD linear fit and slope. Default method is Welch periodogram.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `lf::Matrix{Float64}`: linear fit
  * `s::Vector{Float64}`: slope of linear fit
  * `pf::Vector{Float64}`: range of frequencies for the linear fit


```
psd_slope(obj; ch, f, norm, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate PSD linear fit and slope. Default method is Welch periodogram.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `f::Tuple{Real, Real}=(0, sr(obj) / 2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
  * `norm::Bool=false`: normalize do dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `lf::Array{Float64, 3}`: linear fit
  * `ls::Array{Float64, 2}`: slope of linear fit
  * `pf::Vector{Float64}`: range of frequencies for the linear fit

<a id='NeuroAnalyzer.sef' href='#NeuroAnalyzer.sef'>#</a>
**`NeuroAnalyzer.sef`** &mdash; *Function*.



```julia
sef(s; x, fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate spectral edge frequency (SEF) – the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

**Arguments**

  * `s::AbstractVector`
  * `x::Float64=0.95`: threshold
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `sef_frq::Float64`: spectral edge frequency


```
sef(s; x, fs, f, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate spectral edge frequency (SEF) – the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

**Arguments**

  * `s::AbstractArray`
  * `x::Float64=0.95`: threshold
  * `fs::Int64`: sampling rate
  * `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `sef_frq::Matrix{Float64}`: spectral edge frequency


```
sef(obj; ch, x, f, method, nt, wlen, woverlap)
```

Calculate spectral edge frequency (SEF) – the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `x::Float64=0.95`: threshold
  * `f::Tuple{Real, Real}=(0, sr(obj) / 2)`: lower and upper frequency bounds, default is total power
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `sef_frq::Matrix{Float64}`: spectral edge frequency

<a id='NeuroAnalyzer.spectrogram' href='#NeuroAnalyzer.spectrogram'>#</a>
**`NeuroAnalyzer.spectrogram`** &mdash; *Function*.



```julia
spectrogram(s; fs, norm, mt, st, wlen, woverlap, w)
```

Calculate spectrogram. Default method is short time Fourier transform.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling frequency
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:stft`: method used to calculate PSD:

      * `:stft`: short time Fourier transform
      * `:mt`: multi-tapered periodogram
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length, default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `sp::Matrix{Float64}`: powers
  * `sf::Vector{Float64}`: frequencies
  * `st::Vector{Float64}`: time


```
spectrogram(obj; ch, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)
```

Calculate spectrogram. Default method is short time Fourier transform.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `pad::Int64=0`: number of zeros to add
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `nt::Int64=7`: number of Slepian tapers
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `sp::Array{Float64, 3}`: powers
  * `sf::Vector{Float64}`: frequencies (frequency indices for continuous wavelet transformation)
  * `st::Vector{Float64}`: time points

<a id='NeuroAnalyzer.spectrum' href='#NeuroAnalyzer.spectrum'>#</a>
**`NeuroAnalyzer.spectrum`** &mdash; *Function*.



```julia
spectrum(s; pad, norm)
```

Calculate FFT, amplitudes, powers and phases.

**Arguments**

  * `s::AbstractVector`
  * `pad::Int64=0`: number of zeros to add at the end of the signal
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `ft::Vector{ComplexF64}`: Fourier transforms
  * `sa::Vector{Float64}`: amplitudes
  * `sp::Vector{Float64}`: powers
  * `sph::Vector{Float64}`: phases


```
spectrum(s; pad, h, norm)
```

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

**Arguments**

  * `s::AbstractArray`
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `h::Bool=false`: use Hilbert transform for calculations instead of FFT
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
  * `sa::Array{Float64, 3}`: amplitudes
  * `sp::Array{Float64, 3}`: powers
  * `sph::Array{Float64, 3}: phase angles


```
spectrum(obj; ch, pad, h, norm)
```

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `pad::Int64=0`: number of zeros to add signal for FFT
  * `h::Bool=false`: use Hilbert transform for calculations instead of FFT
  * `norm::Bool=false`: normalize do dB

**Returns**

Named tuple containing:

  * `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
  * `sa::Array{Float64, 3}`: amplitudes
  * `sp::Array{Float64, 3}`: powers
  * `sph::Array{Float64, 3}: phase angles

<a id='NeuroAnalyzer.total_power' href='#NeuroAnalyzer.total_power'>#</a>
**`NeuroAnalyzer.total_power`** &mdash; *Function*.



```julia
total_power(s; fs, method, nt, wlen, woverlap, w, frq_n, frq, fs, ncyc)
```

Calculate total power.

**Arguments**

  * `s::AbstractVector`
  * `fs::Int64`: sampling rate
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `tp::Float64`: total power


```
total_power(s; fs, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate total power.

`# Arguments

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, fs / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `tp::Matrix{Float64}`: total power


```
total_power(obj, ch, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate total power.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

  * `tp::Matrix{Float64}`: total power


<a id='Analyze:-time'></a>

<a id='Analyze:-time-1'></a>

### Analyze: time

<a id='NeuroAnalyzer.acor' href='#NeuroAnalyzer.acor'>#</a>
**`NeuroAnalyzer.acor`** &mdash; *Function*.



acor(s; l, demean, biased, method)

Calculate auto-correlation.

**Arguments**

  * `s::AbstractVector`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating auto-correlation:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`, `biased` value is ignored
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acor(s; l, demean, biased, method)

Calculate auto-correlation.

**Arguments**

  * `s::AbstractMatrix`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating auto-correlation:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`, `biased` value is ignored
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acor(s; l, demean, biased, method)

Calculate auto-correlation.

**Arguments**

  * `s::AbstractArray`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating auto-correlation:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`, `biased` value is ignored
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acor(obj; ch, lag, demean, biased, method)

Calculate auto-correlation. For ERP return trial-averaged auto-correlation.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `l::Real=1`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating auto-correlation:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`, `biased` value is ignored
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

Named tuple containing:

  * `ac::Array{Float64, 3}`
  * `l::Vector{Float64}`: lags [s]

<a id='NeuroAnalyzer.acov' href='#NeuroAnalyzer.acov'>#</a>
**`NeuroAnalyzer.acov`** &mdash; *Function*.



acov(s; l, demean, biased, method)

Calculate autocovariance.

**Arguments**

  * `s::AbstractVector`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing autocovariance
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating autocovariance:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cov`: `acf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `autocov()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acov(s; l, demean, biased, method)

Calculate autocovariance.

**Arguments**

  * `s::AbstractMatrix`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing autocovariance
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating autocovariance:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acov(s; l, demean, biased, method)

Calculate autocovariance.

**Arguments**

  * `s::AbstractArray`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
  * `demean::Bool=true`: demean signal before computing autocovariance
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating autocovariance:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

  * `ac::Matrix{Float64}`


acov(obj; ch, l, demean, biased, method)

Calculate autocovariance. For ERP return trial-averaged autocovariance.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `l::Int64=1`: lags range is `0:lag` [samples]
  * `demean::Bool=true`: demean signal before computing autocovariance
  * `biased::Bool=true`: calculate biased or unbiased autocovariance
  * `method::Symbol=:sum`: method of calculating autocovariance:

      * `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `autocor()`, `biased` value is ignored

**Returns**

Named tuple containing:

  * `ac::Array{Float64, 3}`
  * `l::Vector{Float64}`: lags [s]

<a id='NeuroAnalyzer.axc2frq' href='#NeuroAnalyzer.axc2frq'>#</a>
**`NeuroAnalyzer.axc2frq`** &mdash; *Function*.



axc2frq(c, l)

Detect peaks in auto-/cross- correlation/covariance and transform them into frequencies.

**Arguments**

  * `c::AbstractVector`: auto-/cross- correlation/covariance values
  * `l::AbstractVector`: lags

**Returns**

  * `frq::Vector{Float64}`: list of frequencies dominating in the auto-/cross- correlation/covariance

<a id='NeuroAnalyzer.pacor' href='#NeuroAnalyzer.pacor'>#</a>
**`NeuroAnalyzer.pacor`** &mdash; *Function*.



pacor(s; l, demean, method)

Calculate partial auto-correlation.

**Arguments**

  * `s::AbstractVector`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `method::Symbol=:yw`: method of calculating auto-correlation:

      * `:yw`: computes the partial autocorrelations using the Yule-Walker equations
      * `:reg`: computes the partial autocorrelations via successive regression models

**Returns**

  * `pac::Matrix{Float64}`

**Notes**

If you get `ERROR: PosDefException: matrix is not positive definite; Cholesky factorization failed.`, try lowering `l` value or change method to `:yw`.


pacor(s; l, demean, method)

Calculate partial auto-correlation.

**Arguments**

  * `s::AbstractMatrix`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `method::Symbol=:yw`: method of calculating auto-correlation:

      * `:yw`: computes the partial autocorrelations using the Yule-Walker equations
      * `:reg`: computes the partial autocorrelations via successive regression models

**Returns**

  * `pac::Matrix{Float64}`


pacor(s; l, demean, method)

Calculate partial auto-correlation.

**Arguments**

  * `s::AbstractArray`
  * `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `method::Symbol=:yw`: method of calculating auto-correlation:

      * `:yw`: computes the partial autocorrelations using the Yule-Walker equations
      * `:reg`: computes the partial autocorrelations via successive regression models

**Returns**

  * `pac::Matrix{Float64}`


pacor(obj; ch, lag, demean, method)

Calculate partial auto-correlation. For ERP return trial-averaged auto-correlation.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `l::Real=1`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing auto-correlation
  * `method::Symbol=:yw`: method of calculating auto-correlation:

      * `:yw`: computes the partial autocorrelations using the Yule-Walker equations
      * `:reg`: computes the partial autocorrelations via successive regression models

**Returns**

Named tuple containing:

  * `pac::Array{Float64, 3}`
  * `l::Vector{Float64}`: lags [s]

<a id='NeuroAnalyzer.xcor' href='#NeuroAnalyzer.xcor'>#</a>
**`NeuroAnalyzer.xcor`** &mdash; *Function*.



xcor(s1, s2; l, demean, biased, method)

Calculate cross-correlation.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `l::Int64=round(Int64, min(length(s1) - 1, 10 * log10(length(s1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-correlation
  * `biased::Bool=true`: calculate biased or unbiased cross-correlation
  * `method::Symbol=:sum`: method of calculating cross-correlation:

      * `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ (std(s1) × std(s2))`
      * `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`, `biased` value is ignored
      * `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

**Returns**

  * `xc::Matrix{Float64}`


xcor(s1, s2; l, demean, biased, method)

Calculate cross-correlation.

**Arguments**

  * `s1::AbstractMatrix`
  * `s2::AbstractMatrix`
  * `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-correlation
  * `biased::Bool=true`: calculate biased or unbiased cross-correlation
  * `method::Symbol=:sum`: method of calculating cross-correlation:

      * `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
      * `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

**Returns**

  * `xc::Array{Float64, 3}`


xcor(s1, s2; l, demean, biased, method)

Calculate cross-correlation.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `l::Int64=round(Int64, min(size(s1[1, :, 1], 1) - 1, 10 * log10(size(s1[1, :, 1], 1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-correlation
  * `biased::Bool=true`: calculate biased or unbiased cross-correlation
  * `method::Symbol=:sum`: method of calculating cross-correlation:

      * `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
      * `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

**Returns**

  * `xc::Array{Float64, 3}`


```
xcor(obj1, obj2; ch1, ch2, ep1, ep2, l, demean, biased, method)
```

Calculate cross-correlation. For ERP return trial-averaged cross-correlation.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
  * `l::Real=1`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-correlation
  * `biased::Bool=true`: calculate biased or unbiased cross-correlation
  * `method::Symbol=:sum`: method of calculating cross-correlation:

      * `:sum`: `acf = Σ(s1[1:end - l] .* s1[1+l:end]) ./ var(s)`
      * `:cor`: `acf = cor(s1[1:end - l], s2[1+l:end])`
      * `:stat`: use StatsBase `crosscor()`, `biased` value is ignored

**Returns**

Named tuple containing:

  * `xc::Array{Float64, 3}`: cross-correlation
  * `l::Vector{Float64}`: lags [s]

<a id='NeuroAnalyzer.xcov' href='#NeuroAnalyzer.xcov'>#</a>
**`NeuroAnalyzer.xcov`** &mdash; *Function*.



xcov(s1, s2; l, demean, biased, method)

Calculate cross-covariance.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `l::Int64=round(Int64, min(length(s1) - 1, 10 * log10(length(s1))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-covariance
  * `biased::Bool=true`: calculate biased or unbiased cross-covariance
  * `method::Symbol=:sum`: method of calculating cross-covariance:

      * `:sum`: `xcf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cov`: `xcf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `crosscov()`, `biased` value is ignored

**Returns**

  * `xc::Matrix{Float64}`


xcov(s1, s2; l, demean, biased, method)

Calculate cross-covariance.

**Arguments**

  * `s1::AbstractMatrix`
  * `s2::AbstractMatrix`
  * `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-covariance
  * `biased::Bool=true`: calculate biased or unbiased cross-covariance
  * `method::Symbol=:sum`: method of calculating cross-covariance:

      * `:sum`: `xcf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cov`: `xcf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `crosscov()`, `biased` value is ignored

**Returns**

  * `xc::Array{Float64, 3}`


xcov(s1, s2; l, demean, biased, method)

Calculate cross-covariance.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `l::Int64=round(Int64, min(size(s1, 2), 10 * log10(size(s1, 2))))`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-covariance
  * `biased::Bool=true`: calculate biased or unbiased cross-covariance
  * `method::Symbol=:sum`: method of calculating cross-covariance:

      * `:sum`: `xcf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cov`: `xcf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `crosscov()`, `biased` value is ignored

**Returns**

  * `xc::Array{Float64, 3}`


```
xcov(obj1, obj2; ch1, ch2, ep1, ep2, l, demean, biased, method)
```

Calculate cross-covariance. For ERP return trial-averaged cross-covariance.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
  * `l::Real=1`: lags range is `-l:l`
  * `demean::Bool=true`: demean signal before computing cross-covariance
  * `biased::Bool=true`: calculate biased or unbiased cross-covariance
  * `method::Symbol=:sum`: method of calculating cross-covariance:

      * `:sum`: `xcf = Σ(s[1:end - l] .* s[1+l:end])`
      * `:cov`: `xcf = cov(s[1:end - l], s[1+l:end])`
      * `:stat`: use StatsBase `crosscov()`, `biased` value is ignored

**Returns**

Named tuple containing:

  * `xc::Array{Float64, 3}`: cross-covariance
  * `l::Vector{Float64}`: lags [s]


<a id='Analyze:-phase'></a>

<a id='Analyze:-phase-1'></a>

### Analyze: phase

<a id='NeuroAnalyzer.cph' href='#NeuroAnalyzer.cph'>#</a>
**`NeuroAnalyzer.cph`** &mdash; *Function*.



```julia
cph(s1, s2; fs)
```

Calculate cross-phases.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Vector{Float64}`: cross-power spectrum phase (in radians)
  * `f::Vector{Float64}`: cross-power spectrum frequencies


```
cph(s; fs)
```

Calculate cross-phases.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
  * `f::Vector{Float64}`: cross-power spectrum frequencies


```
cph(s1, s2; fs)
```

Calculate cross-phases.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
  * `f::Vector{Float64}`: cross-power spectrum frequencies


```
cph(obj; ch)
```

Calculate cross-phases.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
  * `f::Vector{Float64, 4}`: cross-power spectrum frequencies


```
cph(obj1, obj2; ch1, ch2, ep1, ep2, norm)
```

Calculate cross-phases.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 3}`: cross-power spectrum phase (in radians)
  * `f::Vector{Float64, 3}`: cross-power spectrum frequencies

<a id='NeuroAnalyzer.ispc' href='#NeuroAnalyzer.ispc'>#</a>
**`NeuroAnalyzer.ispc`** &mdash; *Function*.



```julia
ispc(s1, s2)
```

Calculate ISPC (Inter-Site-Phase Clustering) between `s1` and `s2`.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

Named tuple containing:

  * `ispc_value::Float64`: ISPC value
  * `ispc_angle::Float64`: ISPC angle
  * `s_diff::Vector{Float64}`: signal difference (s2 - s1)
  * `ph_diff::Vector{Float64}`: phase difference (s2 - s1)
  * `s1_phase::Vector{Float64}`: signal 1 phase
  * `s2_phase::Vector{Float64}`: signal 2 phase


```
ispc(obj; ch)
```

Calculate ISPCs (Inter-Site-Phase Clustering).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `ispc_value::Array{Float64, 3}`: ISPC value matrices over epochs
  * `ispc_angle::Array{Float64, 3}`: ISPC angle matrices over epochs


```
ispc(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate ISPC (Inter-Site-Phase Clustering).

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `ispc_value::Array{Float64, 2}`: ISPC value
  * `ispc_angle::Array{Float64, 2}`: ISPC angle
  * `s_diff::Array{Float64, 3}`: signal difference (s2 - s1)
  * `ph_diff::Array{Float64, 3}`: phase difference (s2 - s1)
  * `s1_phase::Array{Float64, 3}`: signal 1 phase
  * `s2_phase::Array{Float64, 3}`: signal 2 phase

<a id='NeuroAnalyzer.itpc' href='#NeuroAnalyzer.itpc'>#</a>
**`NeuroAnalyzer.itpc`** &mdash; *Function*.



```julia
itpc(s; t, w)
```

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

**Arguments**

  * `s::AbstractArray`: one channel over epochs
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_value::Float64`: ITPC value
  * `itpcz_value::Float64`: Rayleigh's ITPC Z value
  * `itpc_angle::Float64`: ITPC angle
  * `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs


```
itpc(obj; <keyword arguments>)
```

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `t::Int64`: time point (sample number) at which ITPC is calculated
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_value::Vector{Float64}`: ITPC or wITPC value
  * `itpcz_value::Vector{Float64}`: Rayleigh's ITPC Z value
  * `itpc_angle::Vector{Float64}`: ITPC angle
  * `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)

<a id='NeuroAnalyzer.itpc_spec' href='#NeuroAnalyzer.itpc_spec'>#</a>
**`NeuroAnalyzer.itpc_spec`** &mdash; *Function*.



```julia
itpc_spec(s; w)
```

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

**Arguments**

  * `s::AbstractArray`: one channel over epochs
  * `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_values::Vector{Float64}`: ITPC values
  * `itpcz_values::Vector{Float64}`: Rayleigh's ITPC Z values
  * `itpc_angles::Vector{Float64}`: ITPC angles
  * `itpc_phases::Matrix{Float64}`: phases at time `t` averaged across trials/epochs


```
itpc_spec(obj; <keyword arguments>)
```

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Int64`
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for the spectrogram
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

**Returns**

Named tuple containing:

  * `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
  * `itpcz_s::Array{Float64, 3}`: spectrogram itpcz_value values
  * `itpc_f::Vector{Float64}`: frequencies list

<a id='NeuroAnalyzer.phdiff' href='#NeuroAnalyzer.phdiff'>#</a>
**`NeuroAnalyzer.phdiff`** &mdash; *Function*.



```julia
phdiff(s1, s2; pad, h)
```

Calculate phase difference between signals.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `pad::Int64=0`: number of zeros to add
  * `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

**Returns**

Named tuple containing:

  * `phd::Vector{Float64}`: phase differences in radians


```
phdiff(s; ch, pad, h)
```

Calculate phase difference between channels and mean phase of reference `ch`.

**Arguments**

  * `s::AbstractArray`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=size(s, 1)`: index of reference channels, default is all  channels except the analyzed one
  * `avg::Symbol=:phase`: method of averaging:

      * `:phase`: phase is calculated for each reference channel separately and then averaged
      * `:signal`: signals are averaged prior to phase calculation
  * `pad::Int64=0`: pad signals with 0s
  * `h::Bool=false`: use FFT or Hilbert transformation

**Returns**

  * `phd::Array{Float64, 3}`


```
phdiff(obj; ch, pad, h)
```

Calculate phase difference between channels and mean phase of reference `ch`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one
  * `avg::Symbol=:phase`: method of averaging:

      * `:phase`: phase is calculated for each reference channel separately and then averaged
      * `:signal`: signals are averaged prior to phase calculation
  * `pad::Int64=0`: pad signals with 0s
  * `h::Bool=false`: use FFT or Hilbert transformation

**Returns**

  * `phd::Array{Float64, 3}`

<a id='NeuroAnalyzer.phsd' href='#NeuroAnalyzer.phsd'>#</a>
**`NeuroAnalyzer.phsd`** &mdash; *Function*.



```julia
phsd(s; fs)
```

Calculate phase spectral density.

**Arguments**

  * `s::Vector{Float64}`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Vector{Float64}`: phases
  * `f::Vector{Float64}`: frequencies


```
phsd(s; fs)
```

Calculate phase spectral density.

**Arguments**

  * `s::AbstractMatrix`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 2}`: phases
  * `f::Vector{Float64}`: frequencies


```
phsd(s; fs)
```

Calculate phase spectral density.

**Arguments**

  * `s::AbstractArray`
  * `fs::Int64`: sampling rate

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 3}`: phases
  * `f::Vector{Float64}`: frequencies


```
phsd(obj; ch)
```

Calculate phase spectral density.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `ph::Array{Float64, 3}`: phases
  * `f::Vector{Float64}`: frequencies

<a id='NeuroAnalyzer.pli' href='#NeuroAnalyzer.pli'>#</a>
**`NeuroAnalyzer.pli`** &mdash; *Function*.



```julia
pli(s1, s2)
```

Calculate PLI (Phase-Lag Index).

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

Named tuple containing:

  * `pv::Float64`: PLI value
  * `sd::Vector{Float64}`: signal difference (s2 - s1)
  * `phd::Vector{Float64}`: phase difference (s2 - s1)
  * `s1ph::Vector{Float64}`: signal 1 phase
  * `s2ph::Vector{Float64}`: signal 2 phase


```
pli(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate PLI (Phase Lag Index).

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `pv::Array{Float64, 2}`: PLI value
  * `sd::Array{Float64, 3}`: signal difference (s2 - s1)
  * `phd::Array{Float64, 3}`: phase difference (s2 - s1)
  * `s1ph::Array{Float64, 3}`: signal 1 phase
  * `s2ph::Array{Float64, 3}`: signal 2 phase


```
pli(obj; ch)
```

Calculate PLIs (Phase Lag Index).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `pv::Array{Float64, 3}`: PLI value matrices over epochs


<a id='Analyze:-ERP'></a>

<a id='Analyze:-ERP-1'></a>

### Analyze: ERP

<a id='NeuroAnalyzer.amp_at' href='#NeuroAnalyzer.amp_at'>#</a>
**`NeuroAnalyzer.amp_at`** &mdash; *Function*.



```julia
amp_at(obj; t)
```

Calculate amplitude at given time.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `t::Real`: time in seconds

**Returns**

  * `p::Matrix{Float64, 2}`: amplitude for each channel per epoch

<a id='NeuroAnalyzer.avgamp_at' href='#NeuroAnalyzer.avgamp_at'>#</a>
**`NeuroAnalyzer.avgamp_at`** &mdash; *Function*.



```julia
avgamp_at(obj; t)
```

Calculate average amplitude at given time segment.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `t::Tuple{Real, Real}`: time segment in seconds

**Returns**

  * `p::Matrix{Float64, 2}`: mean amplitude for each channel per epoch

<a id='NeuroAnalyzer.erp_peaks' href='#NeuroAnalyzer.erp_peaks'>#</a>
**`NeuroAnalyzer.erp_peaks`** &mdash; *Function*.



```julia
erp_peaks(obj)
```

Detect a pair of positive and negative peaks of ERP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position

<a id='NeuroAnalyzer.maxamp_at' href='#NeuroAnalyzer.maxamp_at'>#</a>
**`NeuroAnalyzer.maxamp_at`** &mdash; *Function*.



```julia
maxamp_at(obj; t)
```

Calculate maximum amplitude at given time segment.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `t::Tuple{Real, Real}`: time segment in seconds

**Returns**

  * `p::Matrix{Float64, 2}`: maximum amplitude for each channel per epoch

<a id='NeuroAnalyzer.minamp_at' href='#NeuroAnalyzer.minamp_at'>#</a>
**`NeuroAnalyzer.minamp_at`** &mdash; *Function*.



```julia
minamp_at(obj; t)
```

Calculate minimum amplitude at given time segment.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `t::Tuple{Real, Real}`: time segment in seconds

**Returns**

  * `p::Matrix{Float64, 2}`: minimum amplitude for each channel per epoch


<a id='Analyze:-MEP'></a>

<a id='Analyze:-MEP-1'></a>

### Analyze: MEP

<a id='NeuroAnalyzer.mep_peaks' href='#NeuroAnalyzer.mep_peaks'>#</a>
**`NeuroAnalyzer.mep_peaks`** &mdash; *Function*.



```julia
mep_peaks(obj)
```

Detect a pair of positive and negative peaks of MEP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position


<a id='Analyze:-misc'></a>

<a id='Analyze:-misc-1'></a>

### Analyze: misc

<a id='NeuroAnalyzer.corm' href='#NeuroAnalyzer.corm'>#</a>
**`NeuroAnalyzer.corm`** &mdash; *Function*.



corm(s; norm=true)

Calculate correlation matrix of `s * s'`.

**Arguments**

  * `s::AbstractVector`
  * `norm::Bool`: normalize correlation matrix

**Returns**

  * `cm::Matrix{Float64}`


corm(s1, s2; norm=true)

Calculate correlation matrix of `s1 * s2'`.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `norm::Bool`: normalize correlation matrix

**Returns**

  * `cm::Matrix{Float64}`


corm(s; norm=true)

Calculate correlation matrix.

**Arguments**

  * `s::AbstractArray`
  * `norm::Bool=false`: normalize covariance

**Returns**

  * `cm::Array{Float64, 4}`


```
corm(obj; ch, norm)
```

Calculate correlation matrix.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `norm::Bool=true`: normalize matrix

**Returns**

  * `cm::Array{Float64, 3}`: correlation matrix for each epoch

<a id='NeuroAnalyzer.covm' href='#NeuroAnalyzer.covm'>#</a>
**`NeuroAnalyzer.covm`** &mdash; *Function*.



covm(s; norm=true)

Calculate covariance matrix of `s * s'`.

**Arguments**

  * `s::AbstractVector`
  * `norm::Bool=false`: normalize covariance

**Returns**

  * `cm::Matrix{Float64}`


covm(s1, s2; norm=true)

Calculate covariance matrix of `s1 * s2'`.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `norm::Bool=false`: normalize covariance

**Returns**

  * `cm::Matrix{Float64}`


covm(s; norm=true)

Calculate covariance matrix.

**Arguments**

  * `s::AbstractArray`
  * `norm::Bool=false`: normalize covariance

**Returns**

  * `cm::Matrix{Float64}`


```
covm(obj; ch, norm)
```

Calculate covariance matrix of `signal * signal'`.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `norm::Bool=false`: normalize matrix

**Returns**

  * `cm::Array{Float64, 3}`: covariance matrix for each epoch

<a id='NeuroAnalyzer.diss' href='#NeuroAnalyzer.diss'>#</a>
**`NeuroAnalyzer.diss`** &mdash; *Function*.



```julia
diss(s1, s2)
```

Calculate DISS (global dissimilarity) and spatial correlation between `s1` and `s2`.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

Named tuple containing:

  * `gd::Float64`: global dissimilarity
  * `sc::Float64`: spatial correlation


```
diss(s)
```

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

**Arguments**

  * `s::AbstractArray`

**Returns**

Named tuple containing:

  * `gd::Array{Float64, 3}`: global dissimilarity
  * `sc::Array{Float64, 3}`: spatial correlation


```
diss(s1, s2)
```

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`

**Returns**

Named tuple containing:

  * `gd::Array{Float64, 3}`: global dissimilarity
  * `sc::Array{Float64, 3}`: spatial correlation


```
diss(obj; ch)
```

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `gd::Array{Float64, 3}`: global dissimilarity
  * `sc::Array{Float64, 3}`: spatial correlation


```
diss(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate DISS (global dissimilarity) and spatial correlation (`ch1` of `obj1` vs `ch2` of `obj2`)..

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `gd::Array{Float64, 3}`: global dissimilarity
  * `sc::Array{Float64, 3}`: spatial correlation

<a id='NeuroAnalyzer.entropy' href='#NeuroAnalyzer.entropy'>#</a>
**`NeuroAnalyzer.entropy`** &mdash; *Function*.



```julia
entropy(s)
```

Calculate entropy.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `ent::Float64`
  * `sent::Float64`: Shanon entropy
  * `leent::Float64`: log energy entropy


```
entropy(s)
```

Calculate entropy.

**Arguments**

  * `s::AbstractArray`

**Returns**

Named tuple containing:

  * `ent::Array{Float64, 2}`
  * `sent::Array{Float64, 2}`: Shanon entropy
  * `leent::Array{Float64, 2}`: log energy entropy


```
entropy(obj; ch)
```

Calculate entropy.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `ent::Array{Float64, 2}`
  * `sent::Array{Float64, 2}`: Shanon entropy
  * `leent::Array{Float64, 2}`: log energy entropy

<a id='NeuroAnalyzer.env_cor' href='#NeuroAnalyzer.env_cor'>#</a>
**`NeuroAnalyzer.env_cor`** &mdash; *Function*.



```julia
env_cor(env1, env2)
```

Calculate envelope correlation.

**Arguments**

  * `env1::Array{Float64, 3}`
  * `env2::Array{Float64, 3}`

**Returns**

Named tuple containing:

  * `ec::Vector{Float64}`: power correlation value
  * `p::Vector{Float64}`: p-value

<a id='NeuroAnalyzer.env_lo' href='#NeuroAnalyzer.env_lo'>#</a>
**`NeuroAnalyzer.env_lo`** &mdash; *Function*.



```julia
env_lo(s, x; d)
```

Calculate lower envelope.

**Arguments**

  * `s::AbstractVector`: signal
  * `x::AbstractVector`: x-axis points
  * `d::Int64=32`: distance between peeks in points, lower values get better envelope fit

**Returns**

  * `e::Vector{Float64}`: envelope

<a id='NeuroAnalyzer.env_up' href='#NeuroAnalyzer.env_up'>#</a>
**`NeuroAnalyzer.env_up`** &mdash; *Function*.



```julia
env_up(s, x; d)
```

Calculate upper envelope.

**Arguments**

  * `s::AbstractVector`: signal
  * `x::AbstractVector`: x-axis points
  * `d::Int64=32`: distance between peeks in points, lower values get better envelope fit

**Returns**

  * `e::Vector{Float64}`: envelope

<a id='NeuroAnalyzer.ged' href='#NeuroAnalyzer.ged'>#</a>
**`NeuroAnalyzer.ged`** &mdash; *Function*.



```julia
ged(s1, s2)
```

Perform generalized eigendecomposition.

**Arguments**

  * `s1::AbstractArray`: signal to be analyzed
  * `s2::AbstractArray`: original signal

**Returns**

Named tuple containing:

  * `sged::Matrix{Float64}`
  * `ress::Vector{Float64}`
  * `ress_norm::Vector{Float64}`: RESS normalized to -1..1


```
ged(obj1, obj2; ch1, ch2, ep1, ep2)
```

Perform generalized eigendecomposition.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: signal data to be analyzed
  * `obj2::NeuroAnalyzer.NEURO`: original signal data
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

  * `sged::Array{Float64, 3}`
  * `ress::Matrix{Float64}`
  * `ress_norm::Matrix{Float64}`

<a id='NeuroAnalyzer.gfp' href='#NeuroAnalyzer.gfp'>#</a>
**`NeuroAnalyzer.gfp`** &mdash; *Function*.



```julia
gfp(s)
```

Calculate GFP (Global Field Power).

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `gfp::Float64`

<a id='NeuroAnalyzer.gfp_norm' href='#NeuroAnalyzer.gfp_norm'>#</a>
**`NeuroAnalyzer.gfp_norm`** &mdash; *Function*.



```julia
gfp_norm(s)
```

Calculate signal normalized for GFP (Global Field Power).

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `gfp_norm::Float64`

<a id='NeuroAnalyzer.henv' href='#NeuroAnalyzer.henv'>#</a>
**`NeuroAnalyzer.henv`** &mdash; *Function*.



```julia
henv(obj; ch, d)
```

Calculate Hilbert spectrum amplitude envelope.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.henv_lo' href='#NeuroAnalyzer.henv_lo'>#</a>
**`NeuroAnalyzer.henv_lo`** &mdash; *Function*.



```julia
henv_lo(s)
```

Calculate lower envelope using Hilbert transform.

**Arguments**

  * `s::AbstractVector`: signal

**Returns**

  * `e::Vector{Float64}`: envelope

<a id='NeuroAnalyzer.henv_mean' href='#NeuroAnalyzer.henv_mean'>#</a>
**`NeuroAnalyzer.henv_mean`** &mdash; *Function*.



```julia
henv_mean(obj; ch, dims, d)
```

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
  * `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  * `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.henv_median' href='#NeuroAnalyzer.henv_median'>#</a>
**`NeuroAnalyzer.henv_median`** &mdash; *Function*.



```julia
henv_median(obj; ch, dims, d)
```

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
  * `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  * `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.henv_up' href='#NeuroAnalyzer.henv_up'>#</a>
**`NeuroAnalyzer.henv_up`** &mdash; *Function*.



```julia
henv_up(s)
```

Calculate upper envelope using Hilbert transform.

**Arguments**

  * `s::AbstractVector`: signal

**Returns**

  * `e::Vector{Float64}`: envelope

<a id='NeuroAnalyzer.hrv_analyze' href='#NeuroAnalyzer.hrv_analyze'>#</a>
**`NeuroAnalyzer.hrv_analyze`** &mdash; *Function*.



```julia
hrv_analyze(nn_seg)
```

Analyze heart rate variability (HRV).

**Arguments**

  * `nn_seg::Vector{Float64}`: list of NN segments [msec]

**Returns**

Named tuple containing:

  * `menn::Float64`: the mean of NN segments
  * `mdnn::Float64`: the median of NN segments
  * `vnn::Float64`: the variance of NN segments
  * `sdnn::Float64`: the standard deviation of NN segments
  * `rmssd::Float64`: ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent NNs
  * `sdsd::Float64`: ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent NNs
  * `nn50::Float64`: the number of pairs of successive NNs that differ by more than 50 ms
  * `pnn50::Float64`, the proportion of NN50 divided by total number of NNs
  * `nn20::Float64`, the number of pairs of successive NNs that differ by more than 20 ms
  * `pnn20::Float64`, the proportion of NN20 divided by total number of NNs

<a id='NeuroAnalyzer.hrv_detect' href='#NeuroAnalyzer.hrv_detect'>#</a>
**`NeuroAnalyzer.hrv_detect`** &mdash; *Function*.



```julia
hrv_detect(obj)
```

Detect heart rate variability (HRV). Requires ECG channel (which will be automatically detected based on `channel_type` field).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

**Returns**

  * `nn_seg::Vector{Float64}`: list of NN segments [msec]
  * `r_idx::Vector{Float64}`: index of R peaks

<a id='NeuroAnalyzer.mutual_information' href='#NeuroAnalyzer.mutual_information'>#</a>
**`NeuroAnalyzer.mutual_information`** &mdash; *Function*.



```julia
mutual_information(s1, s2)
```

Calculate mutual information.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

  * `mutual_information::Float64`


```
mutual_information(s1, s2)
```

Calculate mutual information (channels of `s1` vs channels of `s2`).

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`

**Returns**

  * `mutual_information::Array{Float64}`


```
mutual_information(s)
```

Calculate mutual information (channels vs channels).

**Arguments**

  * `s::AbstractArray`

**Returns**


```
mutual_information(obj; channel)
```

Calculate mutual information between channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `mutual_information::Array{Float64, 3}`


```
mutual_information(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate mutual information between two channels.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

  * `m::Array{Float64, 3}`

<a id='NeuroAnalyzer.negentropy' href='#NeuroAnalyzer.negentropy'>#</a>
**`NeuroAnalyzer.negentropy`** &mdash; *Function*.



```julia
negentropy(signal)
```

Calculate negentropy.

**Arguments**

  * `signal::AbstractVector`

**Returns**

  * `negent::Float64`


```
negentropy(s)
```

Calculate negentropy.

**Arguments**

  * `s::AbstractArray`

**Returns**

  * `ne::Array{Float64, 2}`


```
negentropy(obj; ch)
```

Calculate negentropy.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

  * `ne::Array{Float64, 2}`

<a id='NeuroAnalyzer.penv' href='#NeuroAnalyzer.penv'>#</a>
**`NeuroAnalyzer.penv`** &mdash; *Function*.



```julia
penv(obj; ch, d, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate power (in dB) envelope.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `p_env::Array{Float64, 3}`: power spectrum envelope
  * `p_env_frq::Vector{Float64}`: frequencies for each envelope

<a id='NeuroAnalyzer.penv_mean' href='#NeuroAnalyzer.penv_mean'>#</a>
**`NeuroAnalyzer.penv_mean`** &mdash; *Function*.



```julia
penv_mean(obj; ch, dims, d, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate power (in dB) envelope: mean and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroAnalyzer.penv_median' href='#NeuroAnalyzer.penv_median'>#</a>
**`NeuroAnalyzer.penv_median`** &mdash; *Function*.



```julia
penv_median(obj; ch, dims, d, method, nt, wlen, woverlap, w, frq_n, frq, ncyc)
```

Calculate power (in dB) envelope: median and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
  * `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_n::Int64=_tlength((0, sr(obj) / 2))`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

**Returns**

Named tuple containing:

  * `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
  * `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  * `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  * `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)

<a id='NeuroAnalyzer.senv' href='#NeuroAnalyzer.senv'>#</a>
**`NeuroAnalyzer.senv`** &mdash; *Function*.



```julia
senv(obj; ch, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)
```

Calculate spectral envelope.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `pad::Int64=0`: number of zeros to add
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `nt::Int64=7`: number of Slepian tapers
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `s_env::Array{Float64, 3}`: spectral envelope
  * `s_env_t::Vector{Float64}`: spectrogram time

<a id='NeuroAnalyzer.senv_mean' href='#NeuroAnalyzer.senv_mean'>#</a>
**`NeuroAnalyzer.senv_mean`** &mdash; *Function*.



```julia
senv_mean(obj; ch, dims, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)
```

Calculate spectral envelope: mean and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `pad::Int64=0`: number of zeros to add
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `nt::Int64=7`: number of Slepian tapers
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: mean
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroAnalyzer.senv_median' href='#NeuroAnalyzer.senv_median'>#</a>
**`NeuroAnalyzer.senv_median`** &mdash; *Function*.



```julia
senv_median(obj; ch, dims, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)
```

Calculate spectral envelope: median and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  * `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `pad::Int64=0`: number of zeros to add
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `norm::Bool=true`: normalize powers to dB
  * `nt::Int64=7`: number of Slepian tapers
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window

**Returns**

Named tuple containing:

  * `s_env_m::Array{Float64, 3}`: spectral envelope: median
  * `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  * `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  * `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)

<a id='NeuroAnalyzer.stationarity' href='#NeuroAnalyzer.stationarity'>#</a>
**`NeuroAnalyzer.stationarity`** &mdash; *Function*.



```julia
stationarity(obj; ch, window, method)
```

Calculate stationarity.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `window::Int64=10`: time window in samples
  * `method::Symbol=:euclid`: stationarity method:

      * `:mean`: mean across `window`-long windows
      * `:var`: variance across `window`-long windows
      * `:cov`: covariance stationarity based on Euclidean distance between covariance matrix of adjacent time windows
      * `:hilbert`: phase stationarity using Hilbert transformation
      * `:adf`: Augmented Dickey–Fuller test; returns ADF-test value and p-value (H0: signal is non-stationary; p-value < alpha means that signal is stationary)

**Returns**

  * `stationarity::Array{Float64, 3}`

<a id='NeuroAnalyzer.stationarity_hilbert' href='#NeuroAnalyzer.stationarity_hilbert'>#</a>
**`NeuroAnalyzer.stationarity_hilbert`** &mdash; *Function*.



```julia
stationarity_hilbert(s)
```

Calculate phase stationarity using Hilbert transformation.

**Arguments**

  * `s::AbstractVector`

**Returns**

  * `stph::Vector{Float64}`

<a id='NeuroAnalyzer.stationarity_mean' href='#NeuroAnalyzer.stationarity_mean'>#</a>
**`NeuroAnalyzer.stationarity_mean`** &mdash; *Function*.



```julia
stationarity_mean(s; window)
```

Calculate mean stationarity. Signal is split into `window`-long windows and averaged across windows.

**Arguments**

  * `s::AbstractVector`
  * `window::Int64`: time window in samples

**Returns**

  * `stm::Vector{Float64}`

<a id='NeuroAnalyzer.stationarity_var' href='#NeuroAnalyzer.stationarity_var'>#</a>
**`NeuroAnalyzer.stationarity_var`** &mdash; *Function*.



```julia
stationarity_var(s; window)
```

Calculate variance stationarity. Signal is split into `window`-long windows and variance is calculated across windows.

**Arguments**

  * `s::AbstractVector`
  * `window::Int64`: time window in samples

**Returns**

  * `stv::Vector{Float64}`

<a id='NeuroAnalyzer.tenv' href='#NeuroAnalyzer.tenv'>#</a>
**`NeuroAnalyzer.tenv`** &mdash; *Function*.



```julia
tenv(obj; ch, d)
```

Calculate temporal envelope (amplitude).

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env::Array{Float64, 3}`: temporal envelope
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.tenv_mean' href='#NeuroAnalyzer.tenv_mean'>#</a>
**`NeuroAnalyzer.tenv_mean`** &mdash; *Function*.



```julia
tenv_mean(obj; ch, dims, d)
```

Calculate temporal envelope (amplitude): mean and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.tenv_median' href='#NeuroAnalyzer.tenv_median'>#</a>
**`NeuroAnalyzer.tenv_median`** &mdash; *Function*.



```julia
tenv_median(obj; ch, dims, d)
```

Calculate temporal envelope (amplitude): median and 95% CI.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  * `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

**Returns**

Named tuple containing:

  * `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
  * `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  * `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  * `s_t::Vector{Float64}`: signal time

<a id='NeuroAnalyzer.tkeo' href='#NeuroAnalyzer.tkeo'>#</a>
**`NeuroAnalyzer.tkeo`** &mdash; *Function*.



```julia
tkeo(s, t; method)
```

Calculate Teager-Kaiser energy-tracking operator.

**Arguments**

  * `s::AbstractVector`: signal
  * `t::AbstractVector=collect(1:length(s))`: time points
  * `method::Symbol=:pow`:

      * `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
      * `:der`: TKEO = f'(t) - f(t) × f''(t)
      * `:amp`: TKEO = envelope(amplitude)^2

**Returns**

  * `tk::Vector{Float64}`


```
tkeo(s, t; method)
```

Calculate Teager-Kaiser energy-tracking operator

**Arguments**

  * `s::AbstractArray`: signal
  * `t::AbstractArray=collect(1:length(s))`: time points
  * `method::Symbol=:pow`:

      * `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
      * `:der`: TKEO = f'(t) - f(t) × f''(t)
      * `:amp`: TKEO = envelope(amplitude)^2

**Returns**

  * `tk::Array{Float64, 3}`


```
tkeo(obj; channel, method)
```

Calculate Teager-Kaiser energy-tracking operator.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `method::Symbol=:pow`:

      * `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
      * `:der`: TKEO = f'(t) - f(t) × f''(t)
      * `:amp`: TKEO = envelope(amplitude)^2

**Returns**

  * `tk::Array{Float64, 3}`


<a id='Plot'></a>

<a id='Plot-1'></a>

## Plot


!!! warning "Missing docstring."
    Missing docstring for `NeuroAnalyzer.add_locs`. Check Documenter's build log for details.


<a id='NeuroAnalyzer.add_to_canvas' href='#NeuroAnalyzer.add_to_canvas'>#</a>
**`NeuroAnalyzer.add_to_canvas`** &mdash; *Function*.



```julia
add_to_canvas(c1, c2; x, y, title, view, file_name)
```

Place CairoSurfaceBase at another canvas at `x, y`. If `file_name` is provided, the plot is saved as PNG file.

**Arguments**

  * `c1::Cairo.CairoSurfaceBase{UInt32}`
  * `c2::Cairo.CairoSurfaceBase{UInt32}`
  * `x::Int64`
  * `y::Int64`
  * `title::String=""`: title of the subplot
  * `view::Bool=true`: view the output image
  * `file_name::String=""`: output image filename

**Returns**

  * `c::Cairo.CairoSurfaceBase{UInt32}`

<a id='NeuroAnalyzer.add_topmargin_canvas' href='#NeuroAnalyzer.add_topmargin_canvas'>#</a>
**`NeuroAnalyzer.add_topmargin_canvas`** &mdash; *Function*.



```julia
add_topmargin_canvas(c1, c2)
```

Resize CairoSurfaceBase to make space for another canvas

**Arguments**

  * `c1::Cairo.CairoSurfaceBase{UInt32}`
  * `c2::Cairo.CairoSurfaceBase{UInt32}`

**Returns**

  * `c::Cairo.CairoSurfaceBase{UInt32}`

<a id='NeuroAnalyzer.plot' href='#NeuroAnalyzer.plot'>#</a>
**`NeuroAnalyzer.plot`** &mdash; *Function*.



```julia
plot(obj; <keyword arguments>)
```

Plot signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ep::Union{Int64, AbstractRange}=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is no label
  * `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or gray palette
  * `emarkers::Bool`: draw epoch markers if available
  * `markers::Bool`: draw markers if available
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `type::Symbol=:normal`: plot type:

      * `:normal`
      * `:mean`: mean ± 95%CI
      * `:butterfly`: butterfly plot
  * `norm::Bool=false`: normalize signal for butterfly and averaged plots
  * `bad::Union{Bool, Matrix{Bool}}=false`: list of bad channels; if not empty – plot bad channels using this list
  * `s_pos::Tuple{Real, Real}=(0, 0)`: draw segment borders if different than (0, 0), used by `iedit()`
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot(obj, c; <keyword arguments>)
```

Plot embedded or external component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `ep::Union{Int64, AbstractRange}=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is no label
  * `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or gray palette
  * `emarkers::Bool`: draw epoch markers if available
  * `markers::Bool`: draw markers if available
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `type::Symbol=:normal`: plot type:

      * `:normal`
      * `:mean`: mean ± 95%CI
      * `:butterfly`: butterfly plot
  * `norm::Bool=false`: normalize signal for butterfly and averaged plots
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot(obj1, obj2; <keyword arguments>)
```

Plot two signals. This function is used to compare two signals, e.g. before and after `ica_recovery()`. Both signals must have the same data type, dimensions and sampling rate.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
  * `ep::Union{Int64, AbstractRange}=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is no label
  * `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
  * `emarkers::Bool`: draw epoch markers if available
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot2canvas' href='#NeuroAnalyzer.plot2canvas'>#</a>
**`NeuroAnalyzer.plot2canvas`** &mdash; *Function*.



```julia
plot2canvas(c)
```

Convert Plots.Plot to CairoSurfaceBase.

**Arguments**

  * `p::Plots.Plot{Plots.GRBackend}`

**Returns**

  * `c::Cairo.CairoSurfaceBase{UInt32}`

<a id='NeuroAnalyzer.plot_2signals' href='#NeuroAnalyzer.plot_2signals'>#</a>
**`NeuroAnalyzer.plot_2signals`** &mdash; *Function*.



```julia
plot_2signals(t, s1, s2; <keyword arguments>)
```

Plot amplitude of single- or multi-channel `s1` and `s2`.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s1::Union{AbstractVector, AbstractArray}`: data to plot (before) - drawn in black
  * `s2::Union{AbstractVector, AbstractArray}`: data to plot (after) - drawn in red
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_bar' href='#NeuroAnalyzer.plot_bar'>#</a>
**`NeuroAnalyzer.plot_bar`** &mdash; *Function*.



```julia
plot_bar(s; <keyword arguments>)
```

Bar plot.

**Arguments**

  * `s::AbstractVector`
  * `xlabels::Vector{String}`: x-ticks labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_box' href='#NeuroAnalyzer.plot_box'>#</a>
**`NeuroAnalyzer.plot_box`** &mdash; *Function*.



```julia
plot_box(s; <keyword arguments>)
```

Box plot.

**Arguments**

  * `s::AbstractArray`
  * `glabels::Vector{String}`: group labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_ci' href='#NeuroAnalyzer.plot_ci'>#</a>
**`NeuroAnalyzer.plot_ci`** &mdash; *Function*.



```julia
plot_ci(s, s_ci_l, s_ci_h; <keyword arguments>)
```

Dots plot.

**Arguments**

  * `s::AbstractVector`: signal
  * `s_l::AbstractVector`: CI lower bound
  * `s_u::AbstractVector`: CI upper bound
  * `t::AbstractVector`: time points
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_coherence' href='#NeuroAnalyzer.plot_coherence'>#</a>
**`NeuroAnalyzer.plot_coherence`** &mdash; *Function*.



```julia
plot_coherence(coh, f; <keyword arguments>)
```

Plot coherence.

**Arguments**

  * `coh::Vector{Float64}`: coherence
  * `f::Vector{Float64}`: frequencies
  * `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Coherence"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_coherence(coh, f; <keyword arguments>)
```

Plot multi-channel coherence.

**Arguments**

  * `coh::Matrix{Float64}`: coherence
  * `f::Vector{Float64}`: frequencies
  * `clabels::Vector{String}=[""]`: channel pairs labels vector
  * `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_coherence_avg' href='#NeuroAnalyzer.plot_coherence_avg'>#</a>
**`NeuroAnalyzer.plot_coherence_avg`** &mdash; *Function*.



```julia
plot_coherence_avg(coh, f; <keyword arguments>)
```

Plot coherence mean and ±95% CI of averaged channels.

**Arguments**

  * `coh::Matrix{Float64}`: coherence
  * `f::Vector{Float64}`: frequencies
  * `clabels::Vector{String}=[""]`: channel pairs labels vector
  * `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the X-axis
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Coherence"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_coherence_butterfly' href='#NeuroAnalyzer.plot_coherence_butterfly'>#</a>
**`NeuroAnalyzer.plot_coherence_butterfly`** &mdash; *Function*.



```julia
plot_coherence_butterfly(coh, f; <keyword arguments>)
```

Butterfly PSD plot.

**Arguments**

  * `coh::Array{Float64, 3}`: coherence
  * `f::Vector{Float64}`: frequencies
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `frq_lim::Tuple{Real, Real}=(f[1], f[end]): frequency limit for the x-axis
  * `xlabel::String="Frequency [Hz]"`: x-axis label
  * `ylabel::String="Coherence"`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_compose' href='#NeuroAnalyzer.plot_compose'>#</a>
**`NeuroAnalyzer.plot_compose`** &mdash; *Function*.



```julia
plot_compose(p; <keyword arguments>)
```

Compose a complex plot of various plots contained in vector `p` using layout `layout`. Layout scheme is:

  * `(2, 2)`: 2 × 2 plots, regular layout
  * `grid(4, 1, heights=[0.6, 0.1, 0.1, 0.1]`: 4 × 1 plots, irregular layout
  * `@layout [a{0.2w} b{0.8w};_ c{0.6}]`: complex layout using Plots.jl `@layout` macro

**Arguments**

  * `p::Vector{Plots.Plot{Plots.GRBackend}}`: vector of plots
  * `layout::Union(Matrix{Any}, Tuple{Int64, Int64}, Plots.GridLayout}`: layout
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for `p` vector plots

**Returns**

  * `pc::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_connections' href='#NeuroAnalyzer.plot_connections'>#</a>
**`NeuroAnalyzer.plot_connections`** &mdash; *Function*.



```julia
plot_connections(obj; <keyword arguments>)
```

Plot connections between channels.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `connections::Matrix{<:Real}`: matrix of connections weights (channels by channels)
  * `threshold::Real=0`: threshold for plotting, see below
  * `threshold_type::Symbol=:neq`: rule for thresholding:

      * `:eq`: plot if connection weight is equal to threshold
      * `:neq`: plot if connection weight is not equal to threshold
      * `:geq`: plot if connection weight is ≥ to threshold
      * `:leq`: plot if connection weight is ≤ to threshold
      * `:g`: plot if connection weight is > to threshold
      * `:l`: plot if connection weight is < to threshold
  * `weights::Bool=true`: weight line widths and alpha based on connection value, if false connections values will be drawn
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `ch_labels::Bool=false`: plot channel labels
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)
  * `title::String=""`: plot title

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_connections(obj; <keyword arguments>)
```

Plot weights at electrode positions.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `connections::Matrix{<:Real}`: matrix of connections weights
  * `threshold::Real=0`: threshold value
  * `threshold_type::Symbol=:neq`: rule for thresholding:

      * `:eq`: plot if connection weight is equal to threshold
      * `:neq`: plot if connection weight is not equal to threshold
      * `:geq`: plot if connection weight is ≥ to threshold
      * `:leq`: plot if connection weight is ≤ to threshold
      * `:g`: plot if connection weight is > to threshold
      * `:l`: plot if connection weight is < to threshold
  * `weights::Bool=true`: weight line widths and alpha based on connection value
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `ch_labels::Bool=false`: plot ch_labels
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)
  * `title::String=""`: plot title

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_dipole2d' href='#NeuroAnalyzer.plot_dipole2d'>#</a>
**`NeuroAnalyzer.plot_dipole2d`** &mdash; *Function*.



```julia
plot_dipole2d(d; <keyword arguments>)
```

Plot dipole in 2D.

**Arguments**

  * `d::NeuroAnalyzer.DIPOLE`

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

**Notes**

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)

<a id='NeuroAnalyzer.plot_dipole3d' href='#NeuroAnalyzer.plot_dipole3d'>#</a>
**`NeuroAnalyzer.plot_dipole3d`** &mdash; *Function*.



```julia
plot_dipole3d(d; <keyword arguments>)
```

Plot dipole in 3D.

**Arguments**

  * `d::NeuroAnalyzer.DIPOLE`
  * `project::Bool=true`: plot lines projected onto X, Y and Z axes

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

**Notes**

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)

<a id='NeuroAnalyzer.plot_dots' href='#NeuroAnalyzer.plot_dots'>#</a>
**`NeuroAnalyzer.plot_dots`** &mdash; *Function*.



```julia
plot_dots(s; <keyword arguments>)
```

Dots plot.

**Arguments**

  * `s::Vector{Vector{Float64}}`
  * `glabels::Vector{String}`: group labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_empty' href='#NeuroAnalyzer.plot_empty'>#</a>
**`NeuroAnalyzer.plot_empty`** &mdash; *Function*.



```julia
plot_empty()
```

Return an empty plot, useful for filling matrices of plots. 

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erop' href='#NeuroAnalyzer.plot_erop'>#</a>
**`NeuroAnalyzer.plot_erop`** &mdash; *Function*.



```julia
plot_erop(p, f; <keyword arguments>)
```

Plot ERO (Event-Related Oscillations) power-spectrum.

**Arguments**

  * `p::AbstractArray`: ERO powers
  * `f::AbstractVector`: ERO frequencies
  * `xlabel::String="default"`
  * `ylabel::String="default"`
  * `title::String="default"`
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_eros' href='#NeuroAnalyzer.plot_eros'>#</a>
**`NeuroAnalyzer.plot_eros`** &mdash; *Function*.



```julia
plot_eros(s, f, t; <keyword arguments>)
```

Plot ERO (Event-Related Oscillations) spectrogram.

**Arguments**

  * `s::AbstractArray`: ERO spectrogram
  * `f::AbstractVector`: ERO frequencies
  * `t::AbstractVector`: ERO time
  * `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
  * `xlabel::String="default"`
  * `ylabel::String="default"`
  * `title::String="default"`
  * `cb::Bool=true`: draw color bar
  * `cb_title::String="Power [dB]"`: color bar title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erp' href='#NeuroAnalyzer.plot_erp'>#</a>
**`NeuroAnalyzer.plot_erp`** &mdash; *Function*.



```julia
plot_erp(t, s, bad; <keyword arguments>)
```

Plot ERP.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractVector`: data to plot
  * `rt::Union{Nothing, Real}=nothing`:: response time value
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_erp(obj; <keyword arguments>)
```

Plot ERP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
  * `xlabel::String="default"`: x-axis label, default is Time [ms]
  * `ylabel::String="default"`: y-axis label, default is Amplitude [units]
  * `title::String="default"`: plot title, default is ERP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
  * `cb::Bool=true`: plot color bar
  * `cb_title::String="default"`: color bar title, default is Amplitude [units]
  * `mono::Bool=false`: use color or gray palette
  * `peaks::Bool=true`: draw peaks
  * `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
  * `type::Symbol=:normal`: plot type: `:normal`, butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
  * `yrev::Bool=false`: reverse Y axis
  * `avg::Bool=false`: plot average ERP for `:butterfly` plot
  * `smooth::Bool=false`: smooth the image using Gaussian blur
  * `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  * `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
  * `sort_epochs::Bool=false`:: sort epochs by rt vector
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erp_avg' href='#NeuroAnalyzer.plot_erp_avg'>#</a>
**`NeuroAnalyzer.plot_erp_avg`** &mdash; *Function*.



```julia
plot_erp_avg(t, s; <keyword arguments>)
```

Plot ERP amplitude mean and ±95% CI.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `rt::Union{Nothing, Real}=nothing`:: response time value
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erp_butterfly' href='#NeuroAnalyzer.plot_erp_butterfly'>#</a>
**`NeuroAnalyzer.plot_erp_butterfly`** &mdash; *Function*.



```julia
plot_erp_butterfly(t, s; <keyword arguments>)
```

Butterfly plot of ERP.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `rt::Union{Nothing, Real}=nothing`:: response time value
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `avg::Bool=false`: plot average ERP
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erp_stack' href='#NeuroAnalyzer.plot_erp_stack'>#</a>
**`NeuroAnalyzer.plot_erp_stack`** &mdash; *Function*.



```julia
plot_erp_stack(s; <keyword arguments>)
```

Plot EPRs stacked by channels or by epochs.

**Arguments**

  * `t::AbstractVector`: x-axis values
  * `s::AbstractArray`
  * `rt::Union{Nothing, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `cb::Bool=true`: plot color bar
  * `cb_title::String=""`: color bar title
  * `mono::Bool=false`: use color or gray palette
  * `smooth::Bool=false`: smooth the image using Gaussian blur
  * `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_erp_topo' href='#NeuroAnalyzer.plot_erp_topo'>#</a>
**`NeuroAnalyzer.plot_erp_topo`** &mdash; *Function*.



```julia
plot_erp_topo(locs, t, s; <keyword arguments>)
```

Plot topographical map ERPs.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `t::Vector{Float64}`: time vector
  * `s::Array{Float64, 2}`: ERPs
  * `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `yrev::Bool=false`: reverse Y axis
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_filter_response' href='#NeuroAnalyzer.plot_filter_response'>#</a>
**`NeuroAnalyzer.plot_filter_response`** &mdash; *Function*.



```julia
plot_filter_response(<keyword arguments>)
```

Plot filter response.

**Arguments**

  * `fs::Int64`: sampling rate
  * `fprototype::Symbol`: filter prototype:

      * `:butterworth`
      * `:chebyshev1`
      * `:chebyshev2`
      * `:elliptic`
      * `:fir`
      * `:iirnotch`: second-order IIR notch filter
      * `:remez`: Remez FIR filter
  * `ftype::Union{Symbol, Nothing}=nothing`: filter type:

      * `:lp`: low pass
      * `:hp`: high pass
      * `:bp`: band pass
      * `:bs`: band stop
  * `cutoff::Union{Real, Tuple{Real, Real}}`: filter cutoff in Hz (tuple for `:bp` and `:bs`)
  * `n::Int64=2560`: signal length in samples
  * `fs::Int64`: sampling rate
  * `order::Int64=8`: filter order (6 dB/octave), number of taps for `:remez`, attenuation (× 4 dB) for `:fir` filters
  * `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for `:elliptic`, 2 dB for others
  * `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for `:elliptic`, 20 dB for others
  * `bw::Real=-1`: bandwidth for `:iirnotch` and :remez filters
  * `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using fred harris' rule-of-thumb
  * `mono::Bool=false`: use color or gray palette
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_histogram' href='#NeuroAnalyzer.plot_histogram'>#</a>
**`NeuroAnalyzer.plot_histogram`** &mdash; *Function*.



```julia
plot_histogram(s; <keyword arguments>)
```

Plot histogram.

**Arguments**

  * `s::AbstractVector`
  * `x::Union{Nothing, Real}=nothing`: value to plot against the histogram
  * `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
  * `bins::Union{Int64, Symbol, AbstractVector}=(length(s) ÷ 10)`: histogram bins: number of bins, range or `:sturges`, `:sqrt`, `:rice`, `:scott` or `:fd`)
  * `label::String=""`: channel label
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_icatopo' href='#NeuroAnalyzer.plot_icatopo'>#</a>
**`NeuroAnalyzer.plot_icatopo`** &mdash; *Function*.



```julia
plot_icatopo(obj; <keyword arguments>)
```

Topographical plot of embedded ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `cb::Bool=false`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_icatopo(obj, ic, ic_mw; <keyword arguments>)
```

Topographical plot of external ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `cb::Bool=false`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_line' href='#NeuroAnalyzer.plot_line'>#</a>
**`NeuroAnalyzer.plot_line`** &mdash; *Function*.



```julia
plot_line(s; <keyword arguments>)
```

Line plot.

**Arguments**

  * `s::AbstractVector`
  * `xlabels::Vector{String}`: x-ticks labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_line(s; <keyword arguments>)
```

Line plot.

**Arguments**

  * `s::AbstractArray`
  * `rlabels::Vector{String}`: signal rows labels
  * `xlabels::Vector{String}`: x-ticks labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_locs' href='#NeuroAnalyzer.plot_locs'>#</a>
**`NeuroAnalyzer.plot_locs`** &mdash; *Function*.



```julia
plot_locs(locs; <keyword arguments>)
```

Preview channel locations.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
  * `ch_labels::Bool=true`: plot channel labels
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `grid::Bool=false`: draw grid, useful for locating positions
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_locs(obj; <keyword arguments>)
```

Preview of channel locations.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: which channel should be highlighted
  * `ch_labels::Bool=true`: plot channel labels
  * `src_labels::Bool=false`: plot source labels
  * `det_labels::Bool=false`: plot detector labels
  * `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `threed::Bool=false`: 3-dimensional plot
  * `mono::Bool=false`: use color or gray palette
  * `grid::Bool=false`: draw grid, useful for locating positions
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)
  * `interactive::Bool=true`: if true, use interactive 3-dimensional plot
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_locs3d' href='#NeuroAnalyzer.plot_locs3d'>#</a>
**`NeuroAnalyzer.plot_locs3d`** &mdash; *Function*.



```julia
plot_locs3d(locs; <keyword arguments>)
```

3D preview of channel locations.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
  * `ch_labels::Bool=true`: plot channel labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
  * `camera::Tuple{Real, Real}=(20, 45)`: camera position – (XY plane angle, XZ plane angle)

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_locs_nirs' href='#NeuroAnalyzer.plot_locs_nirs'>#</a>
**`NeuroAnalyzer.plot_locs_nirs`** &mdash; *Function*.



```julia
plot_locs_nirs(locs; <keyword arguments>)
```

Preview of NIRS optodes and channel locations. It uses Cartesian `:loc_x` and `:loc_y` locations.

**Arguments**

  * `locs::DataFrame`: columns: labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `opt_pairs::Matrix{Int64}`: pairs of source and detector
  * `src_n::Int64`: number of sources
  * `det_n::Int64`: number of detectors
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
  * `src_labels::Bool=false`: plot source labels
  * `det_labels::Bool=false`: plot detector labels
  * `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `grid::Bool=false`: draw grid, useful for locating positions
  * `plot_size::Int64=400`: plot dimensions in pixels (size × size)

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_matrix' href='#NeuroAnalyzer.plot_matrix'>#</a>
**`NeuroAnalyzer.plot_matrix`** &mdash; *Function*.



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
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_mep' href='#NeuroAnalyzer.plot_mep'>#</a>
**`NeuroAnalyzer.plot_mep`** &mdash; *Function*.



```julia
plot_mep(t, s, bad; <keyword arguments>)
```

Plot MEP.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractVector`: data to plot
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_mep(obj; <keyword arguments>)
```

Plot MEP.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
  * `xlabel::String="default"`: x-axis label, default is Time [ms]
  * `ylabel::String="default"`: y-axis label, default is Amplitude [units]
  * `title::String="default"`: plot title, default is MEP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
  * `cb::Bool=true`: plot color bar
  * `cb_title::String="default"`: color bar title, default is Amplitude [units]
  * `mono::Bool=false`: use color or gray palette
  * `peaks::Bool=true`: draw peaks
  * `peaks_detect::Bool=true`: if true, detect MEP peaks, otherwise use embedded
  * `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
  * `type::Symbol=:normal`: plot type: `:normal`, butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
  * `yrev::Bool=false`: reverse Y axis
  * `avg::Bool=false`: plot average MEP for `:butterfly` plot
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_mep_avg' href='#NeuroAnalyzer.plot_mep_avg'>#</a>
**`NeuroAnalyzer.plot_mep_avg`** &mdash; *Function*.



```julia
plot_mep_avg(t, s; <keyword arguments>)
```

Plot MEP amplitude mean and ±95% CI.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_mep_butterfly' href='#NeuroAnalyzer.plot_mep_butterfly'>#</a>
**`NeuroAnalyzer.plot_mep_butterfly`** &mdash; *Function*.



```julia
plot_mep_butterfly(t, s; <keyword arguments>)
```

Butterfly plot of MEP.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `avg::Bool=false`: plot average MEP
  * `yrev::Bool=false`: reverse Y axis
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_mep_stack' href='#NeuroAnalyzer.plot_mep_stack'>#</a>
**`NeuroAnalyzer.plot_mep_stack`** &mdash; *Function*.



```julia
plot_mep_stack(s; <keyword arguments>)
```

Plot EPRs stacked by channels or by epochs.

**Arguments**

  * `t::AbstractVector`: x-axis values
  * `s::AbstractArray`
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `cb::Bool=true`: plot color bar
  * `cb_title::String=""`: color bar title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_paired' href='#NeuroAnalyzer.plot_paired'>#</a>
**`NeuroAnalyzer.plot_paired`** &mdash; *Function*.



```julia
plot_paired(signal; <keyword arguments>)
```

Plot paired data.

**Arguments**

  * `signal::Vector{Vector{Float64}}`
  * `glabels::Vector{String}`: group labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_phsd' href='#NeuroAnalyzer.plot_phsd'>#</a>
**`NeuroAnalyzer.plot_phsd`** &mdash; *Function*.



```julia
plot_phsd(sf, sp; <keyword arguments>)
```

Plot PHSD (phase spectral density).

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Vector{Float64}`:phases
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_phsd(sf, sp; <keyword arguments>)
```

Plot multi-channel PHSD (phase spectral density).

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Matrix{Float64}`:phases
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_phsd(obj; <keyword arguments>)
```

Plot phase spectral density.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Phase [rad]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Phase [rad]
  * `title::String="default"`: plot title, default is PHSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or gray palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_phsd(obj; <keyword arguments>)
```

Plot phase spectral density of embedded or external component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Phase [rad]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Phase [rad]
  * `title::String="default"`: plot title, default is PHSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or gray palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `units::String=""`
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_phsd_3d' href='#NeuroAnalyzer.plot_phsd_3d'>#</a>
**`NeuroAnalyzer.plot_phsd_3d`** &mdash; *Function*.



```julia
plot_phsd_w3d(sf, sp; <keyword arguments>)
```

Plot 3-d waterfall PHSD plot.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`:phases
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_phsd_avg' href='#NeuroAnalyzer.plot_phsd_avg'>#</a>
**`NeuroAnalyzer.plot_phsd_avg`** &mdash; *Function*.



```julia
plot_phsd_avg(sf, sp; <keyword arguments>)
```

Plot PHSD mean and ±95% CI of averaged channels.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`:phases
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_phsd_butterfly' href='#NeuroAnalyzer.plot_phsd_butterfly'>#</a>
**`NeuroAnalyzer.plot_phsd_butterfly`** &mdash; *Function*.



```julia
plot_phsd_butterfly(sf, sp; <keyword arguments>)
```

Butterfly PHSD plot.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`:phases
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_phsd_topo' href='#NeuroAnalyzer.plot_phsd_topo'>#</a>
**`NeuroAnalyzer.plot_phsd_topo`** &mdash; *Function*.



```julia
plot_phsd_topo(locs, sf, sp; <keyword arguments>)
```

Plot topographical map PHSDs.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`:phases
  * `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling: linear-linear (`:linlin`), log10-linear (`:loglin`)
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_polar' href='#NeuroAnalyzer.plot_polar'>#</a>
**`NeuroAnalyzer.plot_polar`** &mdash; *Function*.



```julia
plot_polar(s; <keyword arguments>)
```

Polar plot.

**Arguments**

  * `s::Union{AbstractVector, AbstractArray}`
  * `m::Tuple{Real, Real}=(0, 0)`: major value to plot
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd' href='#NeuroAnalyzer.plot_psd'>#</a>
**`NeuroAnalyzer.plot_psd`** &mdash; *Function*.



```julia
plot_psd(sf, sp; <keyword arguments>)
```

Plot PSD (power spectrum density).

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Vector{Float64}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_psd(sf, sp; <keyword arguments>)
```

Plot multi-channel PSD (power spectrum density).

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Matrix{Float64}`: powers
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end])`: frequency limit for the X-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_psd(obj; <keyword arguments>)
```

Plot power spectrum density.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=fs`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [units^2/Hz]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [units^2/Hz]
  * `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or gray palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_psd(obj; <keyword arguments>)
```

Plot power spectrum density of embedded or external component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:welch`: method used to calculate PSD:

      * `:welch`: Welch periodogram
      * `:fft`: fast Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:stft`: short time Fourier transform
      * `:mw`: Morlet wavelet convolution
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
  * `ref::Symbol=:abs`: type of PSD reference: absolute power (no reference) (`:abs`) or relative to: total power (`:total`), `:delta`, `:theta`, `:alpha`, `:beta`, `:beta_high`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower` or `:gamma_higher`
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `xlabel::String="default"`: x-axis label, default is Frequency [Hz]
  * `ylabel::String="default"`: y-axis label, default is Power [dB] or Power [units^2/Hz]
  * `zlabel::String="default"`: z-axis label for 3-d plots, default is Power [dB] or Power [units^2/Hz]
  * `title::String="default"`: plot title, default is PSD [frequency limit: 0-128 Hz] [channel: 1, epoch: 1, time window: 0 ms:10 s]
  * `mono::Bool=false`: use color or gray palette
  * `type::Symbol=:normal`: plot type: `:normal`, `:butterfly`, `:mean`, 3-d waterfall (`:w3d`), 3-d surface (`:s3d`), topographical (`:topo`)
  * `units::String=""`
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_3d' href='#NeuroAnalyzer.plot_psd_3d'>#</a>
**`NeuroAnalyzer.plot_psd_3d`** &mdash; *Function*.



```julia
plot_psd_w3d(sf, sp; <keyword arguments>)
```

Plot 3-d waterfall PSD plot.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`: powers
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `zlabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `variant::Symbol`: waterfall (`:w`) or surface (`:s`)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_avg' href='#NeuroAnalyzer.plot_psd_avg'>#</a>
**`NeuroAnalyzer.plot_psd_avg`** &mdash; *Function*.



```julia
plot_psd_avg(sf, sp; <keyword arguments>)
```

Plot PSD mean and ±95% CI of averaged channels.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Matrix{Float64}`: powers
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_butterfly' href='#NeuroAnalyzer.plot_psd_butterfly'>#</a>
**`NeuroAnalyzer.plot_psd_butterfly`** &mdash; *Function*.



```julia
plot_psd_butterfly(sf, sp; <keyword arguments>)
```

Butterfly PSD plot.

**Arguments**

  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`: powers
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_psd_topo' href='#NeuroAnalyzer.plot_psd_topo'>#</a>
**`NeuroAnalyzer.plot_psd_topo`** &mdash; *Function*.



```julia
plot_psd_topo(locs, sf, sp; <keyword arguments>)
```

Plot topographical map PSDs.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `sf::Vector{Float64}`: frequencies
  * `sp::Array{Float64, 3}`: powers
  * `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(sf[1], sf[end]): frequency limit for the x-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `ax::Symbol=:linlin`: type of axes scaling:

      * `:linlin`: linear-linear
      * `:loglin`: log10-linear
      * `:linlog`: linear-log10
      * `:loglog`: log10-log10
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_save' href='#NeuroAnalyzer.plot_save'>#</a>
**`NeuroAnalyzer.plot_save`** &mdash; *Function*.



```julia
plot_save(p; file_name::String)
```

Saves plot as file (PNG/PDF). File format is determined using `file_name` extension.

**Arguments**

  * `p::Plots.Plot{Plots.GRBackend}`
  * `file_name::String`

<a id='NeuroAnalyzer.plot_signal' href='#NeuroAnalyzer.plot_signal'>#</a>
**`NeuroAnalyzer.plot_signal`** &mdash; *Function*.



```julia
plot_signal(t, s; <keyword arguments>)
```

Plot amplitude of single- or multi-channel `s`.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::Union{AbstractVector, AbstractArray}`: data to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_signal(t, s, bad; <keyword arguments>)
```

Plot amplitude of single- or multi-channel `s`.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::Union{AbstractVector, AbstractArray}`: data to plot
  * `bad::Vector{Bool}`: list of bad channels
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_signal_avg' href='#NeuroAnalyzer.plot_signal_avg'>#</a>
**`NeuroAnalyzer.plot_signal_avg`** &mdash; *Function*.



```julia
plot_signal_avg(t, signal; <keyword arguments>)
```

Plot amplitude mean and ±95% CI of averaged `signal` channels.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `norm::Bool=false`: normalize to -1 .. +1
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_signal_butterfly' href='#NeuroAnalyzer.plot_signal_butterfly'>#</a>
**`NeuroAnalyzer.plot_signal_butterfly`** &mdash; *Function*.



```julia
plot_signal_butterfly(t, s; <keyword arguments>)
```

Butterfly plot of `s` channels.

**Arguments**

  * `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  * `s::AbstractArray`: data to plot
  * `clabels::Vector{String}=[""]`: signal channel labels vector
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `scale::Bool=true`: draw scale
  * `units::String=""`: units of the scale
  * `mono::Bool=false`: use color or gray palette
  * `norm::Bool=false`: normalize to -1 .. +1
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_spectrogram' href='#NeuroAnalyzer.plot_spectrogram'>#</a>
**`NeuroAnalyzer.plot_spectrogram`** &mdash; *Function*.



```julia
plot_spectrogram(st, sf, sp; <keyword arguments>)
```

Plot single-channel spectrogram.

**Arguments**

  * `st::Vector{Float64}`: time
  * `sf::Vector{<:Real}`: frequencies
  * `sp::Array{Float64, 2}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limit for the Y-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `units::String=""`
  * `smooth::Bool=false`: smooth the image using Gaussian blur
  * `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_spectrogram(sch, sf, sp; <keyword arguments>)
```

Plot multiple-channel spectrogram.

**Arguments**

  * `sch::Vector{String}`: channel labels
  * `sf::Vector{<:Real}`: frequencies
  * `sp::Array{Float64, 2}`: powers
  * `norm::Bool=true`: whether powers are normalized to dB
  * `frq_lim::Tuple{Real, Real}=(0, 0): frequency limit for the Y-axis
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `units::String=""`
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_spectrogram(obj; <keyword arguments>)
```

Plots spectrogram.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
  * `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]

[channel: 1, epoch: 1, time window: 0 ms:10 s]

  * `mono::Bool=false`: use color or gray palette
  * `markers::Bool`: draw markers if available
  * `smooth::Bool=false`: smooth the image using Gaussian blur
  * `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_spectrogram(obj, c; <keyword arguments>)
```

Plots spectrogram of embedded or external component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `ep::Int64=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `norm::Bool=true`: normalize powers to dB
  * `method::Symbol=:stft`: method of calculating spectrogram:

      * `:stft`: short-time Fourier transform
      * `:mt`: multi-tapered periodogram
      * `:mw`: Morlet wavelet convolution
      * `:gh`: Gaussian and Hilbert transform
      * `:cwt`: continuous wavelet transformation
  * `nt::Int64=7`: number of Slepian tapers
  * `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  * `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
  * `w::Bool=true`: if true, apply Hanning window
  * `gw::Real=5`: Gaussian width in Hz
  * `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet
  * `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
  * `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  * `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
  * `frq_lim::Tuple{Real, Real}=(0, 0)`: y-axis limits
  * `xlabel::String="default"`: x-axis label, default is Time [s]
  * `ylabel::String="default"`: y-axis label, default is Frequency [Hz]
  * `title::String="default"`: plot title, default is Spectrogram [frequency limit: 0-128 Hz]

[component: 1, epoch: 1, time window: 0 ms:10 s]

  * `mono::Bool=false`: use color or gray palette
  * `markers::Bool`: draw markers if available
  * `units::String=""`
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_topo' href='#NeuroAnalyzer.plot_topo'>#</a>
**`NeuroAnalyzer.plot_topo`** &mdash; *Function*.



```julia
plot_topo(c; <keyword arguments>)
```

Plot topographical view.

**Arguments**

  * `s::Vector{<:Real}`: values to plot (one value per channel)
  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `head::Bool=true`: draw head
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_topo(obj; <keyword arguments>)
```

Topographical plot.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ep::Union{Int64, AbstractRange}=0`: epoch to display
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or gray palette
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `head::Bool=true`: draw head
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_topo(obj; <keyword arguments>)
```

Topographical plot of embedded or external component.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `ep::Union{Int64, AbstractRange}=0`: epoch to display
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
  * `mono::Bool=false`: use color or gray palette
  * `cb::Bool=true`: plot color bar
  * `cb_label::String="[A.U.]"`: color bar label
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  * `plot_contours::Bools=true`: plot contours over topo plot
  * `plot_electrodes::Bools=true`: plot electrodes over topo plot
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `head::Bool=true`: draw head
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_violin' href='#NeuroAnalyzer.plot_violin'>#</a>
**`NeuroAnalyzer.plot_violin`** &mdash; *Function*.



```julia
plot_violin(s; <keyword arguments>)
```

Violin plot.

**Arguments**

  * `s::AbstractArray`
  * `glabels::Vector{String}`: group labels
  * `xlabel::String=""`: x-axis label
  * `ylabel::String=""`: y-axis label
  * `title::String=""`: plot title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_weights' href='#NeuroAnalyzer.plot_weights'>#</a>
**`NeuroAnalyzer.plot_weights`** &mdash; *Function*.



```julia
plot_weights(locs; <keyword arguments>)
```

Plot weights at electrode positions.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `weights::Vector{<:Real}=[]`: weights vector
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `ch_labels::Bool=true`: plot channel labels
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)
  * `title::String=""`: plot title

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
plot_weights(obj; <keyword arguments>)
```

Plot weights at electrode positions.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `weights::Matrix{<:Real}`: matrix of weights
  * `ch_labels::Bool=false`: plot ch_labels
  * `head::Bool=true`: draw head
  * `head_labels::Bool=false`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  * `plane::Symbol=:xy`: which plane to plot:

      * `:xy`: horizontal (top)
      * `:xz`: coronary (front)
      * `:yz`: sagittal (side)
  * `title::String=""`: plot title

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.plot_xac' href='#NeuroAnalyzer.plot_xac'>#</a>
**`NeuroAnalyzer.plot_xac`** &mdash; *Function*.



```julia
plot_xac(m, lags; <keyword arguments>)
```

Plot cross/auto-covariance/correlation.

**Arguments**

  * `m::Abstractvector`: covariance matrix
  * `lags::AbstractVector`: covariance lags
  * `xlabel::String="lag"`
  * `ylabel::String=""`
  * `title::String=""`
  * `cb_title::String=""`: color bar title
  * `mono::Bool=false`: use color or gray palette
  * `kwargs`: optional arguments for plot() function

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.resize_canvas' href='#NeuroAnalyzer.resize_canvas'>#</a>
**`NeuroAnalyzer.resize_canvas`** &mdash; *Function*.



```julia
resize_canvas(c; r)
```

Resize CairoSurfaceBase by a factor.

**Arguments**

  * `c::Cairo.CairoSurfaceBase{UInt32}`
  * `r::Real`: resizing factor

**Returns**

  * `c::Cairo.CairoSurfaceBase{UInt32}`


<a id='GUI'></a>

<a id='GUI-1'></a>

## GUI

<a id='NeuroAnalyzer.iedit' href='#NeuroAnalyzer.iedit'>#</a>
**`NeuroAnalyzer.iedit`** &mdash; *Function*.



```julia
iedit(obj, ch, mono, zoom, snap)
```

Interactive edit of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `mono::Bool=true`: use color or gray palette
  * `zoom::Real=5`: how many seconds are displayed in one segment
  * `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

<a id='NeuroAnalyzer.iedit_ch' href='#NeuroAnalyzer.iedit_ch'>#</a>
**`NeuroAnalyzer.iedit_ch`** &mdash; *Function*.



```julia
iedit_ch(obj)
```

Interactive edit signal channels properties and locations.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object

<a id='NeuroAnalyzer.iedit_cont' href='#NeuroAnalyzer.iedit_cont'>#</a>
**`NeuroAnalyzer.iedit_cont`** &mdash; *Function*.



```julia
iedit_cont(obj, ch, mono, zoom)
```

Interactive edit of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `mono::Bool=true`: use color or gray palette
  * `zoom::Real=5`: how many seconds are displayed in one segment
  * `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

<a id='NeuroAnalyzer.iedit_ep' href='#NeuroAnalyzer.iedit_ep'>#</a>
**`NeuroAnalyzer.iedit_ep`** &mdash; *Function*.



```julia
iedit_ep(obj, ch, mono)
```

Interactive edit of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `mono::Bool=true`: use color or gray palette

<a id='NeuroAnalyzer.iplot' href='#NeuroAnalyzer.iplot'>#</a>
**`NeuroAnalyzer.iplot`** &mdash; *Function*.



```julia
iplot(obj, ch, zoom)
```

Interactive plot of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iplot(obj1, obj2, ch, zoom)
```

Interactive plot of two continuous or epoched signals.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.iplot_cont' href='#NeuroAnalyzer.iplot_cont'>#</a>
**`NeuroAnalyzer.iplot_cont`** &mdash; *Function*.



```julia
iplot_cont(obj, ch, zoom)
```

Interactive plot of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iplot_cont(obj1, obj2, ch, zoom)
```

Interactive plot of two continuous signals.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.iplot_ep' href='#NeuroAnalyzer.iplot_ep'>#</a>
**`NeuroAnalyzer.iplot_ep`** &mdash; *Function*.



```julia
iplot_ep(obj, ch)
```

Interactive plot of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels


```
iplot_ep(obj1, obj2, ch)
```

Interactive plot of two epoched signal.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels

<a id='NeuroAnalyzer.iplot_icatopo' href='#NeuroAnalyzer.iplot_icatopo'>#</a>
**`NeuroAnalyzer.iplot_icatopo`** &mdash; *Function*.



```julia
iplot_icatopo(obj; <keyword arguments>)
```

Interactive topographical plot of embedded ("ic" and "ic_mw") ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`


```
iplot_icatopo(obj, ic, ic_mw; <keyword arguments>)
```

Interactive topographical plot of external ICA components.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ic::Matrix{Float64}`: components IC(1)..IC(n)
  * `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  * `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
  * `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
  * `amethod::Symbol=:mean`: averaging method:

      * `:mean`
      * `:median`
  * `imethod::Symbol=:sh`: interpolation method:

      * `:sh`: Shepard
      * `:mq`: Multiquadratic
      * `:imq`: InverseMultiquadratic
      * `:tp`: ThinPlate
      * `:nn`: NearestNeighbour
      * `:ga`: Gaussian
  * `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

<a id='NeuroAnalyzer.iplot_locs3d' href='#NeuroAnalyzer.iplot_locs3d'>#</a>
**`NeuroAnalyzer.iplot_locs3d`** &mdash; *Function*.



```julia
iplot_locs3d(locs; <keyword arguments>)
```

3D interactive preview of channel locations.

**Arguments**

  * `locs::DataFrame`: columns: channel, labels, loc*radius, loc*theta, loc*x, loc*y, loc*z, loc*radius*sph, loc*theta*sph, loc*phi_sph
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
  * `ch_labels::Bool=true`: plot channel labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
  * `camera::Tuple{Real, Real}=(20, 45)`: camera position – (XY plane angle, XZ plane angle)


```
iplot_locs3d(obj; <keyword arguments>)
```

3D interactive preview of channel locations.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
  * `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
  * `ch_labels::Bool=true`: plot channel labels
  * `head_labels::Bool=true`: plot head labels
  * `mono::Bool=false`: use color or gray palette
  * `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
  * `camera::Tuple{Real, Real}=(20, 45)`: camera position – (XY plane angle, XZ plane angle)

<a id='NeuroAnalyzer.ipsd' href='#NeuroAnalyzer.ipsd'>#</a>
**`NeuroAnalyzer.ipsd`** &mdash; *Function*.



```julia
ipsd(obj; <keyword arguments>)
```

Interactive PSD of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.ipsd_cont' href='#NeuroAnalyzer.ipsd_cont'>#</a>
**`NeuroAnalyzer.ipsd_cont`** &mdash; *Function*.



```julia
ipsd_cont(obj, ch, zoom)
```

Interactive PSD of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.ipsd_ep' href='#NeuroAnalyzer.ipsd_ep'>#</a>
**`NeuroAnalyzer.ipsd_ep`** &mdash; *Function*.



```julia
ipsd_ep(obj, ch)
```

Interactive PSD of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel

<a id='NeuroAnalyzer.iselect_seg' href='#NeuroAnalyzer.iselect_seg'>#</a>
**`NeuroAnalyzer.iselect_seg`** &mdash; *Function*.



```julia
iselect_seg(m; c, extract, v)
```

Interactive selection of matrix area.

**Arguments**

  * `m::AbstractMatrix`
  * `c::Bool=false`: if true, select circular segment
  * `extract::Bool=false`
  * `v::Bool=false`

**Returns**

  * `r1::Int64`: upper-left corner
  * `r2::Int64`: bottom-right corner
  * `c1::Int64`: upper-left corner
  * `c2::Int64`: bottom-right corner

or

  * `seg::Union{AbstractMatrix, AbstractVector}`: extracted segment

<a id='NeuroAnalyzer.iselect_ts' href='#NeuroAnalyzer.iselect_ts'>#</a>
**`NeuroAnalyzer.iselect_ts`** &mdash; *Function*.



```julia
iselect_ts(obj, ch, mono, zoom, snap)
```

Select time segment.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type=datatype(obj))`: channel(s) to plot, default is EEG/MEG/ERP channels
  * `mono::Bool=true`: use color or gray palette
  * `zoom::Real=5`: how many seconds are displayed in one segment
  * `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

**Returns**

  * `seg::Tuple{Float64, Float64}`

<a id='NeuroAnalyzer.ispectrogram' href='#NeuroAnalyzer.ispectrogram'>#</a>
**`NeuroAnalyzer.ispectrogram`** &mdash; *Function*.



```julia
ispectrogram(obj, ch, zoom)
```

Interactive spectrogram of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.ispectrogram_cont' href='#NeuroAnalyzer.ispectrogram_cont'>#</a>
**`NeuroAnalyzer.ispectrogram_cont`** &mdash; *Function*.



```julia
ispectrogram_cont(obj, ch, zoom)
```

Interactive spectrogram of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.ispectrogram_ep' href='#NeuroAnalyzer.ispectrogram_ep'>#</a>
**`NeuroAnalyzer.ispectrogram_ep`** &mdash; *Function*.



```julia
ispectrogram_ep(obj, ch)
```

Interactive spectrogram of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel

<a id='NeuroAnalyzer.itopo' href='#NeuroAnalyzer.itopo'>#</a>
**`NeuroAnalyzer.itopo`** &mdash; *Function*.



```julia
itopo(obj, ch, seg)
```

Interactive topographical map.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type=datatype(obj))`: channel(s) to plot, default is EEG/MEG/ERP channels
  * `seg::Tuple{Real, Real}=(0, 0)`: segment (from, to) in seconds to display

<a id='NeuroAnalyzer.iview' href='#NeuroAnalyzer.iview'>#</a>
**`NeuroAnalyzer.iview`** &mdash; *Function*.



```julia
iview(obj; ch, zoom)
```

Interactive view of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iview(obj1, obj2; ch, zoom)
```

Interactive view of continuous or epoched signal.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iview(obj, c; c_idx, zoom)
```

Interactive view of embedded or external component of continuous or epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.iview_cont' href='#NeuroAnalyzer.iview_cont'>#</a>
**`NeuroAnalyzer.iview_cont`** &mdash; *Function*.



```julia
iview_cont(obj; ch, zoom)
```

Interactive view of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iview_cont(obj1, obj2; ch, zoom)
```

Interactive view of continuous signal.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
  * `zoom::Real=5`: how many seconds are displayed in one segment


```
iview_cont(obj, c; c_idx, zoom)
```

Interactive view of embedded or external component of continuous signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
  * `zoom::Real=5`: how many seconds are displayed in one segment

<a id='NeuroAnalyzer.iview_ep' href='#NeuroAnalyzer.iview_ep'>#</a>
**`NeuroAnalyzer.iview_ep`** &mdash; *Function*.



```julia
iview_ep(obj, ch, mono)
```

Interactive view of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels


```
iview_ep(obj1, obj2, ch, mono)
```

Interactive view of epoched signal.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels


```
iview_ep(obj, c; c_idx, mono)
```

Interactive view of embedded or external component of epoched signal.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  * `c::Union{Symbol, AbstractArray}`: component to plot
  * `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels

<a id='NeuroAnalyzer.iview_plot' href='#NeuroAnalyzer.iview_plot'>#</a>
**`NeuroAnalyzer.iview_plot`** &mdash; *Function*.



```julia
iview_plot(p)
```

View plot object.

**Arguments**

  * `p::Plots.Plot{Plots.GRBackend}`


```
iview_plot(file_name)
```

View PNG image.

**Arguments**

  * `file_name::String`


```
iview_plot(c)
```

View Cairo surface object.

**Arguments**

  * `c::Cairo.CairoSurfaceBase{UInt32}`


<a id='Statistics'></a>

<a id='Statistics-1'></a>

## Statistics

<a id='NeuroAnalyzer.binom_prob' href='#NeuroAnalyzer.binom_prob'>#</a>
**`NeuroAnalyzer.binom_prob`** &mdash; *Function*.



```julia
binom_prob(p, r, n)
```

Calculate probability of exactly `r` successes in `n` trials.

**Arguments**

  * `p::Float64`: proportion of successes
  * `r::Int64`: number of successes
  * `n::Int64`: number of trials

**Returns**

  * `binom_prob::Float64`: probability

<a id='NeuroAnalyzer.binom_stat' href='#NeuroAnalyzer.binom_stat'>#</a>
**`NeuroAnalyzer.binom_stat`** &mdash; *Function*.



```julia
binom_stat(p, n)
```

Calculate mean and standard deviation for probability `p`.

**Arguments**

  * `p::Float64`: proportion of successes
  * `n::Int64`: number of trials

**Returns**

Named tuple containing:

  * `m::Float64`: mean
  * `s::Float64`: standard deviation

<a id='NeuroAnalyzer.bootstrap_ci' href='#NeuroAnalyzer.bootstrap_ci'>#</a>
**`NeuroAnalyzer.bootstrap_ci`** &mdash; *Function*.



```julia
bootstrap_ci(s; n1, n2, ci)
```

Calculate Confidence Interval using bootstrapping.

**Arguments**

  * `s::AbstractMatrix`: signal (time points × epochs)
  * `n1::Int64=3000`: number of stage 1 resamplings – number of samples in the resampled signal
  * `n2::Int64=1000`: number of stage 2 resamplings – number of samples sampled from the signal
  * `ci::Float64=0.95`: confidence interval

**Returns**

Named tuple containing:

  * `s_avg::Vector{Float64}`: averaged signal
  * `s_ci_l::Vector{Float64}`: lower bound of the confidence interval
  * `s_ci_h::Vector{Float64}`: upper bound of the confidence interval

<a id='NeuroAnalyzer.bootstrap_stat' href='#NeuroAnalyzer.bootstrap_stat'>#</a>
**`NeuroAnalyzer.bootstrap_stat`** &mdash; *Function*.



```julia
bootstrap_stat(s; n1, n2, f)
```

Calculate signal statistic using bootstrapping.

**Arguments**

  * `s::AbstractMatrix`: signal (time points × epochs)
  * `n1::Int64=3000`: number of stage 1 resamplings – number of samples in the resampled signal
  * `n2::Int64=1000`: number of stage 2 resamplings – number of samples sampled from the signal
  * `f::String`: statistic function to be applied, e.g. `f="abs(maximum(OBJ))"; signal is given using variable`OBJ` here.

**Returns**

  * `out::AbstractVector`

<a id='NeuroAnalyzer.channel_stats' href='#NeuroAnalyzer.channel_stats'>#</a>
**`NeuroAnalyzer.channel_stats`** &mdash; *Function*.



```julia
channel_stats(obj)
```

Calculate channels statistics per epoch.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

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

<a id='NeuroAnalyzer.ci2z' href='#NeuroAnalyzer.ci2z'>#</a>
**`NeuroAnalyzer.ci2z`** &mdash; *Function*.



```julia
ci2z(ci_level)
```

Convert Confidence Interval level to z score.

**Arguments**

  * `ci_level::Float64=0.95`: confidence level

**Returns**

  * `z_score::Float64`

<a id='NeuroAnalyzer.ci_median' href='#NeuroAnalyzer.ci_median'>#</a>
**`NeuroAnalyzer.ci_median`** &mdash; *Function*.



```julia
ci_median(x; ci_level)
```

Calculate confidence interval for a median.

**Arguments**

  * `x::AbstractVector`
  * `ci_level::Float64=0.95`: confidence level

**Returns**

  * `ci_median::Tuple(Float64, Float64)`


```
ci_median(x; ci_level)
```

Calculate confidence interval for a median.

**Arguments**

  * `x::AbstractArray`
  * `ci_level::Float64=0.95`: confidence level

**Returns**

  * `ci_median::Tuple(Float64, Float64)`

<a id='NeuroAnalyzer.ci_r' href='#NeuroAnalyzer.ci_r'>#</a>
**`NeuroAnalyzer.ci_r`** &mdash; *Function*.



```julia
ci_r(x, y; ci_level)
```

Calculate confidence interval for a correlation coefficient.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`
  * `ci_level::Float64=0.95`: confidence level

**Returns**

  * `ci_r::Tuple(Float64, Float64)`


```
ci_r(; r, n, ci_level)
```

Calculate confidence interval for a correlation coefficient.

**Arguments**

  * `r::Float64`
  * `n::Int64`: number of observations
  * `ci_level::Float64=0.95`: confidence level

**Returns**

  * `ci_r::Tuple(Float64, Float64)`

<a id='NeuroAnalyzer.cmp_stat' href='#NeuroAnalyzer.cmp_stat'>#</a>
**`NeuroAnalyzer.cmp_stat`** &mdash; *Function*.



```julia
cmp_stat(stat_dist, stat_value)
```

Calculate proportion of elements below or above a given statistic value.

**Arguments**

  * `stat_dist::AbstractVector`: statistic values distribution
  * `stat_value::Real`: statistic value
  * `type::Symbol=:g`: calculation proportion of elements greater (`:g`) or lesser (`:l`) than `stat_value`

**Returns**

  * `p::Float64`

<a id='NeuroAnalyzer.cmp_test' href='#NeuroAnalyzer.cmp_test'>#</a>
**`NeuroAnalyzer.cmp_test`** &mdash; *Function*.



```julia
cmp_test(seg1, seg2, paired, alpha, type, exact)
```

Compare two vectors; Kruskall-Wallis test is used first, next t-test (paired on non-paired) or non-parametric test (paired: Wilcoxon signed rank, non-paired: Mann-Whitney U test) is applied.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`
  * `paired::Bool`
  * `alpha::Float64=0.05`: confidence level
  * `type::Symbol=:auto`: choose test automatically (`:auto`), permutation-based (`:perm`), parametric (`:p`) or non-parametric (`:np`)
  * `exact::Bool=false`: if true, use exact Wilcoxon test
  * `nperm::Int64=1000`: number of permutation for `:perm` method

**Returns**

Named tuple containing for type !== `:perm`:

  * `t`: test results
  * `ts::Tuple{Float64, String}`: test statistics
  * `tc::Tuple{Float64, Float64}`: test statistics confidence interval
  * `df::Int64`: degrees of freedom
  * `p::Float64`: p-value

Named tuple containing for type === `:perm`:

  * `t::Tuple{perm_diff::Vector{Float64}, obs_diff::Float64}`: test results: (permutation difference, observed difference)
  * `p1::Float64`: one-sided p-value
  * `p2::Float64`: two-sided p-value

<a id='NeuroAnalyzer.cor_test' href='#NeuroAnalyzer.cor_test'>#</a>
**`NeuroAnalyzer.cor_test`** &mdash; *Function*.



```julia
cor_test(seg1, seg2)
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

<a id='NeuroAnalyzer.count_thresh' href='#NeuroAnalyzer.count_thresh'>#</a>
**`NeuroAnalyzer.count_thresh`** &mdash; *Function*.



```julia
count_thresh(x; t, t_type)
```

Collect thresholded elements, e.g. in a topographical map.

**Arguments**

  * `x::AbstractMatrix`
  * `t::Real`: threshold value
  * `t_type::Symbol=:g`: rule for thresholding: = (`:eq`), ≥ (`:geq`), ≤ (`:leq`), > (`:g`), < (`:l`)

**Returns**

Named tuple containing:

  * `x_t::Int64`: thresholded matrix
  * `n::Int64`: number of elements

<a id='NeuroAnalyzer.crit_z' href='#NeuroAnalyzer.crit_z'>#</a>
**`NeuroAnalyzer.crit_z`** &mdash; *Function*.



```julia
crit_z(ci_level)
```

Calculate critical Z value.

**Arguments**

  * `c::Float64=0.95`: confidence level

**Returns**

  * `z_score::Float64`

<a id='NeuroAnalyzer.cvar' href='#NeuroAnalyzer.cvar'>#</a>
**`NeuroAnalyzer.cvar`** &mdash; *Function*.



```julia
cvar(se, s)
```

Calculate coefficient of variation for statistic `s`.

**Arguments**

  * `se::Real`: standard error
  * `s::Real`: statistics, e.g. mean value

**Returns**

  * `cvar::Float64`

<a id='NeuroAnalyzer.cvar_mean' href='#NeuroAnalyzer.cvar_mean'>#</a>
**`NeuroAnalyzer.cvar_mean`** &mdash; *Function*.



```julia
cvar_mean(x)
```

Calculate coefficient of variation for a mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `cvar_mean::Float64`

<a id='NeuroAnalyzer.cvar_median' href='#NeuroAnalyzer.cvar_median'>#</a>
**`NeuroAnalyzer.cvar_median`** &mdash; *Function*.



```julia
cvar_median(x)
```

Calculate coefficient of variation for a median.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `cvar_median::Float64`

<a id='NeuroAnalyzer.dprime' href='#NeuroAnalyzer.dprime'>#</a>
**`NeuroAnalyzer.dprime`** &mdash; *Function*.



```julia
dprime(p1::Real, p2::Real)
```

Calculate d' and response bias for two proportions.

**Arguments**

  * `p1::Real`
  * `p2::Real`

**Returns**

Named tuple containing:

  * `dprime::Float64`
  * `rb::Float64`: response bias

<a id='NeuroAnalyzer.dranks' href='#NeuroAnalyzer.dranks'>#</a>
**`NeuroAnalyzer.dranks`** &mdash; *Function*.



```julia
dranks(x, nbins)
```

Calculate ranks scaled in 0..nbins.

**Arguments**

  * `x::AbstractArray`: some continuous variable such as reaction time (the time it takes to indicate the response)
  * `nbins::Int64`: number of bins, default is Sturges' formula

**Returns**

  * `caf::Array{Float64}`

<a id='NeuroAnalyzer.effsize' href='#NeuroAnalyzer.effsize'>#</a>
**`NeuroAnalyzer.effsize`** &mdash; *Function*.



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


```
effsize(p1, p2)
```

Calculate effect size for two proportions `p1` and `p2`.

**Arguments**

  * `p1::Float64`: 1st proportion, e.g. 0.7
  * `p2::Float64`: 2nd proportion, e.g. 0.3

**Returns**

  * `e::Float64`

<a id='NeuroAnalyzer.epoch_stats' href='#NeuroAnalyzer.epoch_stats'>#</a>
**`NeuroAnalyzer.epoch_stats`** &mdash; *Function*.



```julia
epoch_stats(obj)
```

Calculate epochs statistics.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`

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

<a id='NeuroAnalyzer.f1' href='#NeuroAnalyzer.f1'>#</a>
**`NeuroAnalyzer.f1`** &mdash; *Function*.



```julia
f1(; tp, tn, fp, fn)
```

Assess performance of the classification model using F1-score. F1-score value ranges from 0 to 1.

**Arguments**

  * `tp::Int64`: number of true positives
  * `tn::Int64`: number of true negatives
  * `fp::Int64`: number of false positives
  * `fn::Int64`: number of false negatives

**Returns**

Named tuple containing:

  * `f1::Float64`: F1-score
  * `p::Float64`: precision
  * `r::Float64`: recall

**Source**

https://www.statology.org/what-is-a-good-f1-score/

<a id='NeuroAnalyzer.friedman' href='#NeuroAnalyzer.friedman'>#</a>
**`NeuroAnalyzer.friedman`** &mdash; *Function*.



```julia
friedman(m)
```

Estimate Friedman's nonparametric two-way analysis of variance (and Kendall's coefficient of concordance).

**Arguments**

  * `m::AbstractArray`: values × groups

**Returns**

Named tuple containing:

  * `f::Float64`: Friedman
  * `k::Float64`: Kendall
  * `p::Float64`: P-value

**Notes**

  * H0 (Friedman) is that the treatments are equal
  * H0 (Kendall) is that there is agreement between rankings or test results
  * Kendall's coefficient of concordance ranges from 0 to 1, with 0 meaning no agreement across raters (judges)

<a id='NeuroAnalyzer.grubbs' href='#NeuroAnalyzer.grubbs'>#</a>
**`NeuroAnalyzer.grubbs`** &mdash; *Function*.



```julia
grubbs(x; alpha, t)
```

Perform Grubbs test for outlier.

**Arguments**

  * `x::AbstractVector`
  * `alpha::Float64=0.95`
  * `t::Int64=0`: test type:

      * `-1`: test whether the minimum value is an outlier
      * `0`: two-sided test
      * `1`: test whether the maximum value is an outlier

**Returns**

  * `g::Bool`: true: outlier exists, false: there is no outlier

<a id='NeuroAnalyzer.hildebrand_rule' href='#NeuroAnalyzer.hildebrand_rule'>#</a>
**`NeuroAnalyzer.hildebrand_rule`** &mdash; *Function*.



```julia
hildebrand_rule(x)
```

Calculate Hildebrand rule for vector `x`. If H < 0.2 then the vector `x` is symmetrical.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `h::Float64`

<a id='NeuroAnalyzer.infcrit' href='#NeuroAnalyzer.infcrit'>#</a>
**`NeuroAnalyzer.infcrit`** &mdash; *Function*.



```julia
infcrit(m)
```

Calculate Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression model `m`.

**Arguments**

  * `m::StatsModels.TableRegressionModel`: linear regression model

**Returns**

Named tuple containing:

  * `aic::Float64`
  * `bic::Float64`

<a id='NeuroAnalyzer.jaccard_similarity' href='#NeuroAnalyzer.jaccard_similarity'>#</a>
**`NeuroAnalyzer.jaccard_similarity`** &mdash; *Function*.



```julia
jaccard_similarity(x, y)
```

Calculate Jaccard similarity between two vectors.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `j::Float64`

<a id='NeuroAnalyzer.k_categories' href='#NeuroAnalyzer.k_categories'>#</a>
**`NeuroAnalyzer.k_categories`** &mdash; *Function*.



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

<a id='NeuroAnalyzer.linreg' href='#NeuroAnalyzer.linreg'>#</a>
**`NeuroAnalyzer.linreg`** &mdash; *Function*.



```julia
linreg(x, y)
```

Linear regression between two vectors.

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

<a id='NeuroAnalyzer.mcc' href='#NeuroAnalyzer.mcc'>#</a>
**`NeuroAnalyzer.mcc`** &mdash; *Function*.



```julia
mcc(; tp, tn, fp, fn)
```

Assess performance of the classification model using Matthews correlation coefficient (MCC).

MCC’s value ranges from -1 to 1, depending on:

  * a score of -1 denotes a complete discrepancy between expected and actual classes
  * 0 is equivalent to making an entirely arbitrary guess
  * total agreement between expected and actual classes is indicated by a score of 1

**Arguments**

  * `tp::Int64`: number of true positives
  * `tn::Int64`: number of true negatives
  * `fp::Int64`: number of false positives
  * `fn::Int64`: number of false negatives

**Returns**

  * `mcc::Float64`

**Source**

https://finnstats.com/index.php/2022/09/06/assess-performance-of-the-classification-model/

<a id='NeuroAnalyzer.mdiff' href='#NeuroAnalyzer.mdiff'>#</a>
**`NeuroAnalyzer.mdiff`** &mdash; *Function*.



```julia
mdiff(s1, s2; n, method)
```

Calculate mean difference and 95% confidence interval for 2 signals.

**Arguments**

  * `s1::AbstractMatrix`
  * `s2::AbstractMatrix`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:absdiff`:

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `st::Vector{Float64}`
  * `sts::Float64`
  * `p::Float64`


```
mdiff(s1, s2; n, method)
```

Calculate mean difference and its 95% CI between channels.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:absdiff`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `st::Matrix{Float64}`
  * `sts::Vector{Float64}`
  * `p::Vector{Float64}`


```
mdiff(obj1, obj2; ch1, ch2, ep1, ep2, n, method)
```

Calculates mean difference and 95% confidence interval for two channels.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2:NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs
  * `n::Int64`: number of bootstraps
  * `method::Symbol[:absdiff, :diff2int]`

      * `:absdiff`: maximum difference
      * `:diff2int`: integrated area of the squared difference

**Returns**

Named tuple containing:

  * `st::Matrix{Float64}`
  * `sts::Vector{Float64}`
  * `p::Vector{Float64}`

<a id='NeuroAnalyzer.meanc' href='#NeuroAnalyzer.meanc'>#</a>
**`NeuroAnalyzer.meanc`** &mdash; *Function*.



```julia
meanc(x; rad)
```

Calculate circular mean.

**Arguments**

  * `x::AbstractVector`: angles
  * `rad::Bool=false`: angles in radians (`rad=true`) or degrees (`rad=false`)

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.meang' href='#NeuroAnalyzer.meang'>#</a>
**`NeuroAnalyzer.meang`** &mdash; *Function*.



```julia
meang(x)
```

Calculate geometric mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.meanh' href='#NeuroAnalyzer.meanh'>#</a>
**`NeuroAnalyzer.meanh`** &mdash; *Function*.



```julia
meanh(x)
```

Calculate harmonic mean.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.meanw' href='#NeuroAnalyzer.meanw'>#</a>
**`NeuroAnalyzer.meanw`** &mdash; *Function*.



```julia
meanw(x, w)
```

Calculate weighted mean.

**Arguments**

  * `x::AbstractVector`
  * `w::AbstractVector`: weights

**Returns**

  * `m::Float64`

<a id='NeuroAnalyzer.moe' href='#NeuroAnalyzer.moe'>#</a>
**`NeuroAnalyzer.moe`** &mdash; *Function*.



```julia
moe(n)
```

Calculate margin of error for given sample size `n`.

**Arguments**

  * `n::Int64`

**Returns**

  * `moe::Float64`

<a id='NeuroAnalyzer.msci95' href='#NeuroAnalyzer.msci95'>#</a>
**`NeuroAnalyzer.msci95`** &mdash; *Function*.



```julia
msci95(s; n, method)
```

Calculate mean, standard deviation and 95% confidence interval.

**Arguments**

  * `s::AbstractVector`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

**Returns**

Named tuple containing:

  * `sm::Float64`: mean
  * `ss::Float64`: standard deviation
  * `su::Float64`: upper 95% CI
  * `sl::Float64`: lower 95% CI


```
msci95(s; n, method)
```

Calculate mean, standard deviation and 95% confidence interval.

**Arguments**

  * `s::AbstractMatrix`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

**Returns**

Named tuple containing:

  * `sm::Vector{Float64}`: mean
  * `ss::Vector{Float64}`: standard deviation
  * `su::Vector{Float64}`: upper 95% CI
  * `sl::Vector{Float64}`: lower 95% CI


```
msci95(s; n, method)
```

Calculate mean, standard deviation and 95% confidence interval.

**Arguments**

  * `s::AbstractArray`
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal method (`:normal`) or `n`-times boostrapping (`:boot`)

**Returns**

Named tuple containing:

  * `sm::Array{Float64}`: mean
  * `ss::Array{Float64}`: standard deviation
  * `su::Array{Float64}`: upper 95% CI
  * `sl::Array{Float64}`: lower 95% CI


```
msci95(obj; ch, n, method)
```

Calculate mean, standard deviation and 95% confidence interval.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
  * `n::Int64=3`: number of bootstraps
  * `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

**Returns**

Named tuple containing:

  * `sm::Matrix{Float64}`: mean
  * `ss::Matrix{Float64}`: standard deviation
  * `su::Matrix{Float64}`: upper 95% CI
  * `sl::Matrix{Float64}`: lower 95% CI


```
msci95(s1, s2)
```

Calculate mean difference, standard deviation and 95% confidence interval.

**Arguments**

  * `s1::AbstractVector`
  * `s2::AbstractVector`

**Returns**

Named tuple containing:

  * `sm::Float64`: mean
  * `ss::Float64`: standard deviation
  * `su::Float64`: upper 95% CI
  * `sl::Float64`: lower 95% CI


```
msci95(s1, s2)
```

Calculate mean difference, standard deviation and 95% confidence interval.

**Arguments**

  * `s1::AbstractArray`
  * `s2::AbstractArray`

**Returns**

Named tuple containing:

  * `sm::Array{Float64}`: mean
  * `ss::Array{Float64}`: standard deviation
  * `su::Array{Float64}`: upper 95% CI
  * `sl::Array{Float64}`: lower 95% CI


```
msci95(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate mean difference, standard deviation and 95% confidence interval.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2:NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `sm::Matrix{Float64}`: mean
  * `ss::Matrix{Float64}`: standard deviation
  * `su::Matrix{Float64}`: upper 95% CI bound
  * `sl::Matrix{Float64}`: lower 95% CI bound

<a id='NeuroAnalyzer.mscr' href='#NeuroAnalyzer.mscr'>#</a>
**`NeuroAnalyzer.mscr`** &mdash; *Function*.



```julia
mscr(; tp, tn, fp, fn)
```

Assess performance of the classification model using misclassification rate.

**Arguments**

  * `tp::Int64`: number of true positives
  * `tn::Int64`: number of true negatives
  * `fp::Int64`: number of false positives
  * `fn::Int64`: number of false negatives

**Returns**

Named tuple containing:

  * `mr::Float64`: misclassification rate
  * `acc::Float64`: accuracy

**Source**

https://www.statology.org/misclassification-rate/

<a id='NeuroAnalyzer.norminv' href='#NeuroAnalyzer.norminv'>#</a>
**`NeuroAnalyzer.norminv`** &mdash; *Function*.



```julia
norminv(x::Real)
```

Convert probability to a normal distribution with a peak at 0.5.

**Arguments**

  * `x::Real`

**Returns**

  * `norminv::Float64`

<a id='NeuroAnalyzer.outlier_detect' href='#NeuroAnalyzer.outlier_detect'>#</a>
**`NeuroAnalyzer.outlier_detect`** &mdash; *Function*.



```julia
outlier_detect(x; method)
```

Detect outliers.

**Arguments**

  * `x::AbstractVector`
  * `method::Symbol=iqr`: detecting methods:

      * `:iqr`: interquartile range
      * `:z`: z-score
      * `:g`: Grubbs test

**Returns**

  * `o::Vector{Bool}`: index of outliers

<a id='NeuroAnalyzer.power_c1g' href='#NeuroAnalyzer.power_c1g'>#</a>
**`NeuroAnalyzer.power_c1g`** &mdash; *Function*.



```julia
power_c1g(x)
```

Calculate study power for a continuous variable (group 1 vs population).

**Arguments**

  * `m0::Real`: population mean
  * `s0::Real`: population standard deviation
  * `m1::Real`: group 1 mean
  * `n1::Int64`: group 1 sample size
  * `alpha::Float64=0.05`: the probability of type I error

**Returns**

  * `p::Float64`: study power

<a id='NeuroAnalyzer.power_c2g' href='#NeuroAnalyzer.power_c2g'>#</a>
**`NeuroAnalyzer.power_c2g`** &mdash; *Function*.



```julia
power_c2g(x)
```

Calculate study power for a continuous variable (group 1 vs group 2).

**Arguments**

  * `m1::Real`: group 1 mean
  * `s1::Real`: group 1 standard deviation
  * `n1::Int64`: group 1 sample size
  * `m2::Real`: group 2 mean
  * `s2::Real`: group 2 standard deviation
  * `n2::Int64`: group 2 sample size
  * `alpha::Float64=0.05`: the probability of type I error

**Returns**

  * `p::Float64`: study power

<a id='NeuroAnalyzer.power_p1g' href='#NeuroAnalyzer.power_p1g'>#</a>
**`NeuroAnalyzer.power_p1g`** &mdash; *Function*.



```julia
power_p1g(x)
```

Calculate required sample size for a proportion (group 1 vs population).

**Arguments**

  * `p0::Float64`: population incidence
  * `p1::Float64`: group 1 anticipated incidence
  * `n1::Int64`: group 1 sample size
  * `alpha::Float64=0.05`: the probability of type I error

**Returns**

  * `p::Float64`: study power

<a id='NeuroAnalyzer.power_p2g' href='#NeuroAnalyzer.power_p2g'>#</a>
**`NeuroAnalyzer.power_p2g`** &mdash; *Function*.



```julia
power_p2g(x)
```

Calculate required sample size for a proportion (group 1 vs group 2).

**Arguments**

  * `p1::Float64`: group 1 incidence
  * `p2::Float64`: group 2 incidence
  * `n1::Int64`: group 1 sample size
  * `n2::Int64`: group 2 sample size
  * `alpha::Float64=0.05`: the probability of type I error

**Returns**

  * `p::Float64`: study power

<a id='NeuroAnalyzer.prank' href='#NeuroAnalyzer.prank'>#</a>
**`NeuroAnalyzer.prank`** &mdash; *Function*.



```julia
prank(x)
```

Calculate percentile rank.

**Arguments**

  * `x::AbstractVector`: the vector to analyze

**Returns**

  * `p::Vector{Float64}`: percentile ranks

<a id='NeuroAnalyzer.pred_int' href='#NeuroAnalyzer.pred_int'>#</a>
**`NeuroAnalyzer.pred_int`** &mdash; *Function*.



```julia
pred_int(n)
```

Calculates the prediction interval (95% CI adjusted for sample size)

**Arguments**

  * `n::Int64`: sample size

**Returns**

  * `pred_int::Float64`

<a id='NeuroAnalyzer.res_norm' href='#NeuroAnalyzer.res_norm'>#</a>
**`NeuroAnalyzer.res_norm`** &mdash; *Function*.



```julia
res_norm(x, g)
```

Test normal distribution of residuals.

**Arguments**

  * `x::AbstractVector`: data values
  * `g::Vector{Int64}`: group(s) to which each data value belongs

**Returns**

Named tuple containing:

  * `adt_p::Vector{Float64}`: p-values for k-sample Anderson–Darling test vs normal distribution
  * `ks_p::Vector{Float64}`: p-values for one-sample exact Kolmogorov–Smirnov test vs normal distribution

**Notes**

p-values are reported for each group and for the whole sample. If there is only one group, p-values are returned only for the whole sample p-values are reported.

<a id='NeuroAnalyzer.rng' href='#NeuroAnalyzer.rng'>#</a>
**`NeuroAnalyzer.rng`** &mdash; *Function*.



```julia
rng(x)
```

Calculate range.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `rng::Float64`


```
rng(x)
```

Calculate range.

**Arguments**

  * `x::AbstractArray`

**Returns**

  * `rng::Float64`

<a id='NeuroAnalyzer.r_test' href='#NeuroAnalyzer.r_test'>#</a>
**`NeuroAnalyzer.r_test`** &mdash; *Function*.



```julia
r_test(; r1, r2, n1, n2)
```

Test if two correlation coefficients are significantly different.

**Arguments**

  * `r1::Float64`: correlation coefficient, group 1
  * `r2::Float64`: correlation coefficient, group 2
  * `n1::Int64`: number of observations, group 1
  * `n2::Int64`: number of observations, group 2

**Returns**

  * `z_r1r2::Float64`

<a id='NeuroAnalyzer.se' href='#NeuroAnalyzer.se'>#</a>
**`NeuroAnalyzer.se`** &mdash; *Function*.



```julia
se(x)
```

Calculate standard error.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `se::Float64`

<a id='NeuroAnalyzer.seg_extract' href='#NeuroAnalyzer.seg_extract'>#</a>
**`NeuroAnalyzer.seg_extract`** &mdash; *Function*.



```julia
seg_extract(m, r1, c1, r2, c2; v)
```

Extract segment from a matrix.

**Arguments**

  * `m::AbstractMatrix`
  * `rc::NTuple{4, Int64}`: upper-left corner row and column, bottom-right corner row and column
  * `c::Bool=false`: if true, use circular segment; for circular segment the segment is always returned as vector
  * `v::Bool=false`: if true, return as vector (matrix m by rows over columns)

**Returns**

  * `seg::Union{AbstractMatrix, AbstractVector}`

<a id='NeuroAnalyzer.seg_mean' href='#NeuroAnalyzer.seg_mean'>#</a>
**`NeuroAnalyzer.seg_mean`** &mdash; *Function*.



```julia
seg_mean(seg)
```

Calculate mean of a segment (e.g. spectrogram).

**Arguments**

  * `seg::AbstractArray`

**Returns**

  * `seg_mean::Vector{Float64}`: averaged segment


```
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

<a id='NeuroAnalyzer.sem_diff' href='#NeuroAnalyzer.sem_diff'>#</a>
**`NeuroAnalyzer.sem_diff`** &mdash; *Function*.



```julia
sem_diff(x, y)
```

Calculate SEM (standard error of the mean) for the difference of two means.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`

**Returns**

  * `sem_diff::Float64`

<a id='NeuroAnalyzer.size_c1diff' href='#NeuroAnalyzer.size_c1diff'>#</a>
**`NeuroAnalyzer.size_c1diff`** &mdash; *Function*.



```julia
size_c1diff(s1, s2, alpha, power)
```

Calculate required sample size for detecting a difference in a continuous variable (group 1 vs population).

**Arguments**

  * `s0::Real`: population standard deviation
  * `s1::Real`: study standard deviation that we want to detect
  * `two_sided::Bool=true`: if true, the estimation is for two-sided difference
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

  * `n::Int64`: study sample size

<a id='NeuroAnalyzer.size_c1g' href='#NeuroAnalyzer.size_c1g'>#</a>
**`NeuroAnalyzer.size_c1g`** &mdash; *Function*.



```julia
size_c1g(m0, s0, m1, alpha, power)
```

Calculate required sample size for a continuous variable (group 1 vs population).

**Arguments**

  * `m0::Real`: population mean
  * `s0::Real`: population standard deviation
  * `m1::Real`: group 1 mean (expected)
  * `alpha::Float64=0.05`: the probability of type I error
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

Named tuple containing:

  * `n::Int64`: group 1 sample size

<a id='NeuroAnalyzer.size_c2g' href='#NeuroAnalyzer.size_c2g'>#</a>
**`NeuroAnalyzer.size_c2g`** &mdash; *Function*.



```julia
size_c2g(m1, s1, m2, r, alpha, power)
```

Calculate required sample size for a continuous variable (group 1 vs group 2).

**Arguments**

  * `m1::Real`: group 1 mean
  * `s1::Real`: group 1 standard deviation
  * `m2::Real`: group 2 mean (expected)
  * `r::Int64=1`: enrollment ratio – the ratio of group 2 to group 1 enrollment
  * `alpha::Float64=0.05`: the probability of type I error
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

Named tuple containing:

  * `n1::Int64`: group 1 sample size
  * `n2::Int64`: group 2 sample size

<a id='NeuroAnalyzer.size_p1diff' href='#NeuroAnalyzer.size_p1diff'>#</a>
**`NeuroAnalyzer.size_p1diff`** &mdash; *Function*.



```julia
size_p1diff(p0, p1, power)
```

Calculate required sample size for detecting a difference in a proportion (group 1 vs population).

**Arguments**

  * `p0::Real`: population proportion
  * `p1::Real`: study proportion that we want to detect
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

  * `n::Int64`: study sample size (for both study groups)

<a id='NeuroAnalyzer.size_p1g' href='#NeuroAnalyzer.size_p1g'>#</a>
**`NeuroAnalyzer.size_p1g`** &mdash; *Function*.



```julia
size_p1g(x)
```

Calculate required sample size for a proportion (group 1 vs population).

**Arguments**

  * `p0::Float64`: population incidence
  * `p1::Float64`: group anticipated incidence
  * `alpha::Float64=0.05`: the probability of type I error
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

  * `n::Int64`: group 1 sample size

<a id='NeuroAnalyzer.size_p2g' href='#NeuroAnalyzer.size_p2g'>#</a>
**`NeuroAnalyzer.size_p2g`** &mdash; *Function*.



```julia
size_p2g(x)
```

Calculate required sample size for a proportion (group 1 vs group 2).

**Arguments**

  * `p1::Float64`: group 1 anticipated incidence
  * `p2::Float64`: group 2 anticipated incidence
  * `r::Int64=1`: enrollment ratio – the ratio of group 2 to group 1 enrollment
  * `alpha::Float64=0.05`: the probability of type I error
  * `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

**Returns**

Named tuple containing:

  * `n1::Int64`: group 1 sample size
  * `n2::Int64`: group 2 sample size

<a id='NeuroAnalyzer.spec_seg' href='#NeuroAnalyzer.spec_seg'>#</a>
**`NeuroAnalyzer.spec_seg`** &mdash; *Function*.



```julia
spec_seg(sp, st, sf; t, f)
```

Return spectrogram segment.

**Arguments**

  * `sp::Matrix{Float64}`: spectrogram powers
  * `sf::Vector{Float64}`: spectrogram frequencies
  * `st::Vector{Float64}`: spectrogram time
  * `t::Tuple{Real, Real}`: time bounds
  * `f::Tuple{Real, Real}`: frequency bounds

**Returns**

Named tuple containing:

  * `segp::Matrix{Float64}`: powers
  * `segs::Vector{Tuple{Float64, Float64}}`: segment coordinates, for plotting should be converted by `Plots.Shape(segs)`
  * `tidx::Tuple{Real, Real}`: time indices
  * `fidx::Tuple{Real, Real}`: frequency indices


```
spec_seg(sp, sf, st; ch, t, f)
```

Return spectrogram segment.

**Arguments**

  * `sp::AbstractArray`: spectrogram powers
  * `sf::AbstractVector`: spectrogram frequencies
  * `st::AbstractVector`: spectrogram time
  * `ch::Int64`: channel
  * `t::Tuple{Real, Real}`: time bounds
  * `f::Tuple{Real, Real}`: frequency bounds

**Returns**

Named tuple containing:

  * `segp::Array{Float64, 3}`: segment of powers
  * `segs::Vector{Tuple{Float64, Float64}}`: segment coordinates, for plotting should be converted by `Plots.Shape(segs)`
  * `tidx::Tuple{Real, Real}`: time indices
  * `fidx::Tuple{Real, Real}`: frequency indices

<a id='NeuroAnalyzer.summary' href='#NeuroAnalyzer.summary'>#</a>
**`NeuroAnalyzer.summary`** &mdash; *Function*.



```julia
summary(x)
```

Return summary statistics.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `mm::Float64`: mean
  * `s::Float64`: standard deviation
  * `v::Float64`: variance
  * `me::Float64`: median
  * `mo::Float64`: mode


```
summary(x, y)
```

Return summary statistics.

**Arguments**

  * `x::AbstractVector`
  * `y::AbstractVector`
  * `g1::String="1"`: group 1 name
  * `g2::String="2"`: group 2 name

**Returns**

  * `mm1::Float64`: mean
  * `mm2::Float64`: mean
  * `s1::Float64`: standard deviation
  * `s2::Float64`: standard deviation
  * `v1::Float64`: variance
  * `v2::Float64`: variance
  * `me1::Float64`: median
  * `me1::Float64`: median
  * `mo1::Float64`: mode
  * `mo2::Float64`: mode

<a id='NeuroAnalyzer.vartest' href='#NeuroAnalyzer.vartest'>#</a>
**`NeuroAnalyzer.vartest`** &mdash; *Function*.



```julia
vartest(obj; ch)
```

Calculate variance F-test.

**Arguments**

  * `obj::NeuroAnalyzer.NEURO`
  * `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

**Returns**

Named tuple containing:

  * `f::Array{Float64, 3}`
  * `p::Array{Float64, 3}`


```
vartest(obj1, obj2; ch1, ch2, ep1, ep2)
```

Calculate variance F-test.

**Arguments**

  * `obj1::NeuroAnalyzer.NEURO`
  * `obj2::NeuroAnalyzer.NEURO`
  * `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
  * `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
  * `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  * `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

**Returns**

Named tuple containing:

  * `f::Array{Float64, 3}`
  * `p::Array{Float64, 3}`

<a id='NeuroAnalyzer.z2pow' href='#NeuroAnalyzer.z2pow'>#</a>
**`NeuroAnalyzer.z2pow`** &mdash; *Function*.



```julia
z2pow(z)
```

Calculate power for a given Z value.

**Arguments**

  * `z::Real`: Z value

**Returns**

  * `power::Float64`

<a id='NeuroAnalyzer.z_score' href='#NeuroAnalyzer.z_score'>#</a>
**`NeuroAnalyzer.z_score`** &mdash; *Function*.



```julia
z_score(x)
```

Calculate Z-scores for each value of the vector `x`.

**Arguments**

  * `x::AbstractVector`

**Returns**

  * `z_score::Vector{Float64}`


<a id='Study'></a>

<a id='Study-1'></a>

## Study

<a id='NeuroAnalyzer.create_study' href='#NeuroAnalyzer.create_study'>#</a>
**`NeuroAnalyzer.create_study`** &mdash; *Function*.



```julia
study(obj, group)
```

Create NeuroAnalyzer STUDY object.

**Arguments**

  * `obj::Vector{NeuroAnalyzer.NEURO}`
  * `group::Vector{Symbol}`

**Returns**

  * `study::NeuroAnalyzer.STUDY`


!!! warning "Missing docstring."
    Missing docstring for `NeuroAnalyzer.epoch_len`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `NeuroAnalyzer.nchannels`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `NeuroAnalyzer.nepochs`. Check Documenter's build log for details.


<a id='NeuroAnalyzer.obj_n' href='#NeuroAnalyzer.obj_n'>#</a>
**`NeuroAnalyzer.obj_n`** &mdash; *Function*.



```julia
obj_n(study)
```

Return number of NeuroAnalyzer NEURO objects in the study.

**Arguments**

  * `study::NeuroAnalyzer.STUDY`

**Returns**

  * `n::Int64`


!!! warning "Missing docstring."
    Missing docstring for `NeuroAnalyzer.sr`. Check Documenter's build log for details.



<a id='NeuroRecorder'></a>

<a id='NeuroRecorder-1'></a>

## NeuroRecorder

<a id='NeuroAnalyzer.ftt' href='#NeuroAnalyzer.ftt'>#</a>
**`NeuroAnalyzer.ftt`** &mdash; *Function*.



```julia
ftt(; duration, trials, interval, gpio)
```

Perform Finger Tapping Test (FTT) in CLI mode. Use computer keyboard (SPACEBAR key) or switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded. When using computer keyboard, only the number of taps and their time points are recorded; tap durations are set to -1.

**Arguments**

  * `duration::Int64=5`: single trial duration in seconds
  * `trials::Int64=2`: number of trials
  * `interval::Int64=2`: interval between trials in seconds
  * `gpio::Int64=23`: Raspberry Pi GPIO to which the switch is connected (default is GPIO 23 = BOARD 16 pin); set to -1 to disable using GPIO

**Returns**

Named tuple containing:

  * `taps::Vector{Int64}`: number of taps per trial
  * `tap_t::Vector{Vector{Float64}}`: taps time point [s]
  * `tap_d::Vector{Vector{Float64}}`: taps duration [s]
  * `taps_int::Vector{Int64}`: number of taps per trial during intervals
  * `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
  * `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals

<a id='NeuroAnalyzer.iftt' href='#NeuroAnalyzer.iftt'>#</a>
**`NeuroAnalyzer.iftt`** &mdash; *Function*.



```julia
iftt(; duration, trials, interval)
```

Perform Finger Tapping Test (FTT) in GUI mode. Use computer keyboard (SPACEBAR key) or switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded.

**Arguments**

  * `duration::Int64=5`: single trial duration in seconds
  * `trials::Int64=2`: number of trials
  * `interval::Int64=2`: interval between trials in seconds

**Returns**

Named tuple containing:

  * `taps::Vector{Int64}`: number of taps per trial
  * `tap_t::Vector{Vector{Float64}}`: taps time point [s]
  * `tap_d::Vector{Vector{Float64}}`: taps duration [s]
  * `taps_int::Vector{Int64}`: number of taps per trial during intervals
  * `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
  * `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals


<a id='NeuroStim'></a>

<a id='NeuroStim-1'></a>

## NeuroStim

<a id='NeuroAnalyzer.ect_charge' href='#NeuroAnalyzer.ect_charge'>#</a>
**`NeuroAnalyzer.ect_charge`** &mdash; *Function*.



```julia
ect_charge(; pw, pint, pf, duration)
```

Calculate charge administered during ECT.

**Arguments**

  * `pw::Real`: pulse width [ms]
  * `pint::Real`: pulse intensity [mA]
  * `pf::Real`: pulse frequency [Hz]
  * `duration::Real`: stimulation duration [s]

**Returns**

  * `charge::Float64`: charge [mC]

<a id='NeuroAnalyzer.tes_dose' href='#NeuroAnalyzer.tes_dose'>#</a>
**`NeuroAnalyzer.tes_dose`** &mdash; *Function*.



```julia
tes_dose(; current, pad_area, duration)
```

Convert `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

**Arguments**

  * `current::Real`: stimulation current [mA]
  * `pad_area::Real`: electrode pad area [cm²]
  * `duration::Int64`: stimulation duration [s]

**Returns**

Named tuple containing:

  * `charge::Float64`: charge [C]
  * `current_density::Float64`: current density [A/m²]
  * `charge_density::Float64`: delivered charge density [kC/m²]

**Source**

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.

<a id='NeuroAnalyzer.tes_model' href='#NeuroAnalyzer.tes_model'>#</a>
**`NeuroAnalyzer.tes_model`** &mdash; *Function*.



```julia
tes_model(; anode, cathode, anode_curr, cathode_curr)
```

Create model of TES stimulation.

**Arguments**

  * `anode::String`: anode location
  * `cathode::String`: cathode location
  * `anode_curr::Real=2.0`: anode current [mA]
  * `cathode_curr::Real=-2.0`: cathode current [mA]

**Returns**

  * `p::Plots.Plot{Plots.GRBackend}`

**Notes**

This is a very initial version, simplified model – just superficial spread of the electric field

**To do**

Model spread of electric field at the cortical surface – reduce charge for skull resistance

<a id='NeuroAnalyzer.tes_protocol' href='#NeuroAnalyzer.tes_protocol'>#</a>
**`NeuroAnalyzer.tes_protocol`** &mdash; *Function*.



```julia
tes_protocol(; <keyword arguments>)
```

Create TES (tDCS/tACS/tRNS/tPCS) protocol.

**Arguments**

  * `type::Symbol`: stimulation type (`:tDCS`, `:tACS`, `:tRNS`, `:tPCS`)
  * `hd::Bool`: protocol includes HD electrodes
  * `current::Real`: stimulation current [mA]
  * `frequency::Real=0`: stimulation frequency [mA]
  * `anode_size::Tuple{Int64, Int64}`: anode dimensions [mm]
  * `cathode_size::Tuple{Int64, Int64}`: cathode dimensions [mm]
  * `anode_loc::Symbol`: anode location (according to 10-20 Positioning System)
  * `cathode_loc::Symbol`: cathode location (according to 10-20 Positioning System)
  * `duration::Real`: stimulation duration [s]
  * `ramp_in::Real`: stimulation duration [s]
  * `ramp_out::Real`: stimulation duration [s]
  * `sham::Bool`: protocol includes sham stimulations

**Returns**

  * `protocol::Dict`


<a id='NeuroTester'></a>

<a id='NeuroTester-1'></a>

## NeuroTester


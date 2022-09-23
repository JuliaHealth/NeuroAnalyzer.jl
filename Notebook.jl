### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d82bea74-724e-4212-bbf3-e42e241572bb
using Pkg

# ╔═╡ ede1280a-1de5-4c77-ad6d-2a26603b6de6
Pkg.activate(".")

# ╔═╡ 5bcaa3f7-abd5-404a-85b3-b2ea3d7004a7
Pkg.add(url="https://codeberg.org/AdamWysokinski/Simpson.jl")

# ╔═╡ 6c53f862-ace6-4a18-ad8f-4cfd516d8027
using NeuroAnalyzer

# ╔═╡ b91df3a9-fb1d-4f3e-a462-2f7820ec6854
using PlutoUI

# ╔═╡ 53358ed5-e464-4557-85bb-a609051b9963
cd(expanduser("~/Documents/Code/NeuroAnalyzer.jl"))

# ╔═╡ 81ca33a6-5311-48c3-8c91-32482b887760
Pkg.update()

# ╔═╡ ed3bb59b-f1cb-47ce-91fa-5adbd1db1fd5
Pkg.resolve()

# ╔═╡ 52051ec7-ff32-4414-928e-af9572420e94
Pkg.instantiate()

# ╔═╡ 04d52f0f-e04c-4263-8c26-4bcd4a3361c5
na_info()

# ╔═╡ b0aebf11-9d67-4f74-a008-ef03162cc39c
edf = eeg_import_edf("test/eeg-test-edf.edf");

# ╔═╡ eeddafd4-7018-41c6-999a-91d85863c511
eeg_delete_channel!(edf, channel=[17, 18, 22, 23, 24])

# ╔═╡ 938f2be0-819a-4c80-806e-e88f691d77d9
eeg_filter!(edf, fprototype=:iirnotch, cutoff=50, bw=2)

# ╔═╡ 0a73c09f-d93d-4422-a090-4225f0819c58
eeg_filter!(edf, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)

# ╔═╡ 0ad7086d-d0cf-49bd-af89-059d976222dc
eeg_filter!(edf, fprototype=:butterworth, ftype=:lp, cutoff=45, order=8)

# ╔═╡ 80c58c02-37e1-4839-8091-3ddb563ebf37
eeg_reference_car!(edf)

# ╔═╡ 1028224e-38fd-4dde-861e-98d5eb2c85c9
eeg_epochs!(edf, epoch_len=10*eeg_sr(edf))

# ╔═╡ e1db2528-f2a9-4b61-80ad-43554ec46961
eeg_epoch_n(edf)

# ╔═╡ 74a64a34-5021-41a2-aa40-b7d7b058f10e
@bind epoch Slider(1:eeg_epoch_n(edf))

# ╔═╡ 858b526b-6441-474a-bbed-f92d487b596e
eeg_plot_signal(edf, scaled=true, epoch=epoch)

# ╔═╡ Cell order:
# ╠═53358ed5-e464-4557-85bb-a609051b9963
# ╠═d82bea74-724e-4212-bbf3-e42e241572bb
# ╠═81ca33a6-5311-48c3-8c91-32482b887760
# ╠═5bcaa3f7-abd5-404a-85b3-b2ea3d7004a7
# ╠═ede1280a-1de5-4c77-ad6d-2a26603b6de6
# ╠═ed3bb59b-f1cb-47ce-91fa-5adbd1db1fd5
# ╠═52051ec7-ff32-4414-928e-af9572420e94
# ╠═6c53f862-ace6-4a18-ad8f-4cfd516d8027
# ╠═b91df3a9-fb1d-4f3e-a462-2f7820ec6854
# ╠═04d52f0f-e04c-4263-8c26-4bcd4a3361c5
# ╠═b0aebf11-9d67-4f74-a008-ef03162cc39c
# ╠═eeddafd4-7018-41c6-999a-91d85863c511
# ╠═938f2be0-819a-4c80-806e-e88f691d77d9
# ╠═0a73c09f-d93d-4422-a090-4225f0819c58
# ╠═0ad7086d-d0cf-49bd-af89-059d976222dc
# ╠═80c58c02-37e1-4839-8091-3ddb563ebf37
# ╠═1028224e-38fd-4dde-861e-98d5eb2c85c9
# ╠═e1db2528-f2a9-4b61-80ad-43554ec46961
# ╠═74a64a34-5021-41a2-aa40-b7d7b058f10e
# ╠═858b526b-6441-474a-bbed-f92d487b596e

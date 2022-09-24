### A Pluto.jl notebook ###
# v0.19.9

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
Pkg.activate(@__DIR__)

# ╔═╡ 6c53f862-ace6-4a18-ad8f-4cfd516d8027
using NeuroAnalyzer

# ╔═╡ b91df3a9-fb1d-4f3e-a462-2f7820ec6854
using PlutoUI

# ╔═╡ b1583955-b45d-468f-aef5-c11bf76d0384
using Interact

# ╔═╡ e5587167-7e88-4524-bb78-95c4bc11de2d
using Gtk

# ╔═╡ 53358ed5-e464-4557-85bb-a609051b9963
cd(expanduser("~/Documents/Code/NeuroAnalyzer.jl"))		

# ╔═╡ 81ca33a6-5311-48c3-8c91-32482b887760
# ╠═╡ disabled = true
#=╠═╡
Pkg.update()
  ╠═╡ =#

# ╔═╡ ed3bb59b-f1cb-47ce-91fa-5adbd1db1fd5
# ╠═╡ disabled = true
#=╠═╡
Pkg.resolve()
  ╠═╡ =#

# ╔═╡ 52051ec7-ff32-4414-928e-af9572420e94
# ╠═╡ disabled = true
#=╠═╡
Pkg.instantiate()
  ╠═╡ =#

# ╔═╡ 04d52f0f-e04c-4263-8c26-4bcd4a3361c5
na_info()

# ╔═╡ 633c7cc2-dab8-47f6-8eab-0f18fbdbd6a4
file_name = open_dialog("Select file EEG file to open")

# ╔═╡ b0aebf11-9d67-4f74-a008-ef03162cc39c
edf = eeg_import_edf(file_name);

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

# ╔═╡ 32a12ecd-d939-4885-a768-4acfb1e7527a
begin
	epoch_n = eeg_epoch_n(edf)
	md"""
	$(@bind epoch Slider(1:epoch_n, show_value=true))
	$(@bind epoch_del Button("X"))
	"""
end

# ╔═╡ b18893e4-a8aa-44d0-9e28-ba76b0c351d8
eeg_plot_signal(edf, scaled=true, epoch=epoch)

# ╔═╡ 83f55252-55fb-4d88-8016-e53452d114b2
epoch_to_delete_ref = Ref{Int64}()

# ╔═╡ de92968e-54f2-4cfa-a9f8-79b85ca0b188
epoch_to_delete_ref[] = epoch

# ╔═╡ 44e41c41-8141-459c-8372-d01ceb6496fe
let	epoch_del
	println(epoch_to_delete_ref[])
	eeg_delete_epoch!(edf, epoch=epoch_to_delete_ref[])
	epoch_n = eeg_epoch_n(edf)
end

# ╔═╡ Cell order:
# ╠═53358ed5-e464-4557-85bb-a609051b9963
# ╠═d82bea74-724e-4212-bbf3-e42e241572bb
# ╠═81ca33a6-5311-48c3-8c91-32482b887760
# ╠═ede1280a-1de5-4c77-ad6d-2a26603b6de6
# ╠═ed3bb59b-f1cb-47ce-91fa-5adbd1db1fd5
# ╠═52051ec7-ff32-4414-928e-af9572420e94
# ╠═6c53f862-ace6-4a18-ad8f-4cfd516d8027
# ╠═b91df3a9-fb1d-4f3e-a462-2f7820ec6854
# ╠═b1583955-b45d-468f-aef5-c11bf76d0384
# ╠═04d52f0f-e04c-4263-8c26-4bcd4a3361c5
# ╠═e5587167-7e88-4524-bb78-95c4bc11de2d
# ╠═633c7cc2-dab8-47f6-8eab-0f18fbdbd6a4
# ╠═b0aebf11-9d67-4f74-a008-ef03162cc39c
# ╠═eeddafd4-7018-41c6-999a-91d85863c511
# ╠═938f2be0-819a-4c80-806e-e88f691d77d9
# ╠═0a73c09f-d93d-4422-a090-4225f0819c58
# ╠═0ad7086d-d0cf-49bd-af89-059d976222dc
# ╠═80c58c02-37e1-4839-8091-3ddb563ebf37
# ╠═1028224e-38fd-4dde-861e-98d5eb2c85c9
# ╠═32a12ecd-d939-4885-a768-4acfb1e7527a
# ╠═b18893e4-a8aa-44d0-9e28-ba76b0c351d8
# ╠═83f55252-55fb-4d88-8016-e53452d114b2
# ╠═de92968e-54f2-4cfa-a9f8-79b85ca0b188
# ╠═44e41c41-8141-459c-8372-d01ceb6496fe

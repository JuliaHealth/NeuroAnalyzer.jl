### A Pluto.jl notebook ###
# v0.19.14

#> [frontmatter]
#> title = "NeuroAnalyzer"
#> date = "2022-09-25"
#> description = "NeuroAnalyzer: example session"

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
# ╠═╡ show_logs = false
Pkg.activate(@__DIR__)

# ╔═╡ 6c53f862-ace6-4a18-ad8f-4cfd516d8027
using NeuroAnalyzer

# ╔═╡ b91df3a9-fb1d-4f3e-a462-2f7820ec6854
using PlutoUI

# ╔═╡ e5587167-7e88-4524-bb78-95c4bc11de2d
using Gtk

# ╔═╡ 53358ed5-e464-4557-85bb-a609051b9963
cd(expanduser("~/Documents/Code/NeuroAnalyzer.jl"))		

# ╔═╡ 04d52f0f-e04c-4263-8c26-4bcd4a3361c5
na_info()

# ╔═╡ 633c7cc2-dab8-47f6-8eab-0f18fbdbd6a4
file_name = open_dialog("Select file EEG file to open")

# ╔═╡ b0aebf11-9d67-4f74-a008-ef03162cc39c
eeg = eeg_import(file_name);

# ╔═╡ eeddafd4-7018-41c6-999a-91d85863c511
# eeg_delete_channel!(eeg, channel=[17, 18, 22, 23, 24])

# ╔═╡ 938f2be0-819a-4c80-806e-e88f691d77d9
eeg_filter!(eeg, fprototype=:iirnotch, cutoff=50, bw=8)

# ╔═╡ 0a73c09f-d93d-4422-a090-4225f0819c58
eeg_filter!(eeg, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)

# ╔═╡ 0ad7086d-d0cf-49bd-af89-059d976222dc
eeg_filter!(eeg, fprototype=:butterworth, ftype=:lp, cutoff=45, order=8)

# ╔═╡ 80c58c02-37e1-4839-8091-3ddb563ebf37
eeg_reference_car!(eeg)

# ╔═╡ 1028224e-38fd-4dde-861e-98d5eb2c85c9
eeg_epoch!(eeg, epoch_len=10*eeg_sr(eeg))

# ╔═╡ 15ca8e97-7d4f-4993-a13f-bb6fc9e12889
# Check epochs visually, remove bad epochs 

# ╔═╡ 4b1d7251-7944-43ef-86c6-ede1ef3f7526
epoch_n = eeg_epoch_n(eeg)

# ╔═╡ de3bba50-1c87-406e-8415-df9e70e695e6
@bind epoch Slider(1:epoch_n, show_value=true)

# ╔═╡ b18893e4-a8aa-44d0-9e28-ba76b0c351d8
eeg_plot_signal(eeg, channel=1:24, epoch=epoch)

# ╔═╡ 32a12ecd-d939-4885-a768-4acfb1e7527a
@bind options confirm(
    PlutoUI.combine() do Child
	md"""
		Epoch: $(Child("idx", NumberField(1:epoch_n)))

	
		Add to delete list: $(Child("del", CheckBox(default=false)))"""
    end
)

# ╔═╡ 1e72e2e7-bb85-4872-b373-af060b182633
epochs_to_delete = Vector{Int64}()

# ╔═╡ 609aeaae-89b2-4b25-9a9b-60bd63b1f8fb
begin
	options.del && push!(epochs_to_delete, options.idx)
	unique!(sort!(epochs_to_delete));
end

# ╔═╡ 0aee665f-4979-4bd5-9ede-282edbe075c3
md"""
Delete selected epochs: $(@bind del_epochs confirm(CheckBox(default=false)))"""

# ╔═╡ 5b6677ce-f33d-411d-a663-7a4d0f93e9e7
md"""
Clear list of selected epochs: $(@bind clear_list confirm(CheckBox(default=false)))"""

# ╔═╡ 7b903fbf-41b7-4121-814e-91ef00d4b5ca
begin
	if clear_list && length(epochs_to_delete) > 0
		for idx in length(epochs_to_delete):-1:1
			deleteat!(epochs_to_delete, idx)
		end
	println("Epochs: $(eeg_epoch_n(eeg))")
	println("Epochs to delete: $(epochs_to_delete)")
	end
end

# ╔═╡ 0a313dce-8228-4a74-bb8a-08f889afc6a5
begin
	if del_epochs && length(epochs_to_delete) > 0
		eeg_delete_epoch!(eeg, epoch=epochs_to_delete)
		for idx in length(epochs_to_delete):-1:1
			deleteat!(epochs_to_delete, idx)
		end
	println("Epochs: $(eeg_epoch_n(eeg))")
	println("Epochs to delete: $(epochs_to_delete)")
	end
end

# ╔═╡ Cell order:
# ╠═53358ed5-e464-4557-85bb-a609051b9963
# ╠═d82bea74-724e-4212-bbf3-e42e241572bb
# ╠═ede1280a-1de5-4c77-ad6d-2a26603b6de6
# ╠═6c53f862-ace6-4a18-ad8f-4cfd516d8027
# ╠═b91df3a9-fb1d-4f3e-a462-2f7820ec6854
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
# ╠═15ca8e97-7d4f-4993-a13f-bb6fc9e12889
# ╠═4b1d7251-7944-43ef-86c6-ede1ef3f7526
# ╠═de3bba50-1c87-406e-8415-df9e70e695e6
# ╠═b18893e4-a8aa-44d0-9e28-ba76b0c351d8
# ╟─32a12ecd-d939-4885-a768-4acfb1e7527a
# ╠═1e72e2e7-bb85-4872-b373-af060b182633
# ╠═609aeaae-89b2-4b25-9a9b-60bd63b1f8fb
# ╠═0aee665f-4979-4bd5-9ede-282edbe075c3
# ╠═5b6677ce-f33d-411d-a663-7a4d0f93e9e7
# ╠═7b903fbf-41b7-4121-814e-91ef00d4b5ca
# ╠═0a313dce-8228-4a74-bb8a-08f889afc6a5

# NeuroAnalyzer benchmarking

Set --threads according to your system:
```sh
julia -q -O3 -g0 --threads 48 --cpu-target=native Benchmarking.jl
```

## Results

Here are some results:
- [Lenovo T14](t14_benchmark.md)
- [Threadripper 3960X, NVIDIA Quadro P2200, CUDA enabled](3960x_cuda_benchmark.md)
- [Threadripper 3960X, NVIDIA Quadro P2200, CUDA disabled](3960x_nocuda_benchmark.md)

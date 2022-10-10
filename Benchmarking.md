# NeuroAnalyzer benchmarking

Set --threads according to your system:
```sh
julia -O3 -g0 --threads 48 --cpu-target=native benchmarks/na_benchmarking.jl > benchmarks/results.md
```

## Results

Here are some results:
- [Lenovo T14](benchmarks/t14.md) (12 threads)
- [Threadripper 3960X, NVIDIA Quadro P2200, CUDA enabled](benchmarks/3960x_cuda.md) (48 threads)
- [Threadripper 3960X, CUDA disabled](benchmarks/3960x_nocuda.md) (48 threads)
- [Threadripper 3960X, CUDA disabled](benchmarks/3960x_24threads_nocuda.md) (24 threads)
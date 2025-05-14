## DEPENDÊNCIAS
import Pkg;

Pkg.add("Plots")
using Plots

Pkg.add("MAT")
using MAT

Pkg.add("DSP")
using DSP

Pkg.add("FFTW")
using FFTW

Pkg.add("LaTeXStrings")
using LaTeXStrings

include("flattop.jl")

## IMPORTANDO SINAL
sinal = matread(raw"./sinal.mat")["x"]
N = length(sinal)
n = 0:N-1

## Amostragem e Plot
fa = 1e4
Ta = 1/fa
t = collect(0:Ta:1.023e-1)
plot(t, sinal,label="sinal",title="Sinal.mat")
plot!(t, sinaljanelado,label="Sinal janelado")


## JANELA DE BLACKMAN
sinaljanela_blackman = sinal.*blackman(N)
it = 15:25
Xb = fft(sinaljanela_blackman)
Xbdb = amp2db.(abs.(Xb[it]))
yt = 0:20:120
plot(n[it]*fa/N, Xbdb.+80, xlabel=L"$f$ (Hz)", label = "|X[k]|", line = :stem,
 marker = (:circle, 3), title = "Janela de Blackman", ytick = yt,
 yaxis = (formatter = yt -> string(Int(yt-80))))



## JANELA FLATTOP SOBRE O SINAL NO TEMPO
# sinaljanelado = sinal .* flattop(N)



## FFT do sinal
FFT_sinal = fft(sinal)
abs_FFT_sinal = abs.(FFT_sinal)

plot(abs_FFT_sinal[1:50],label="FFT sinal", title="FFTs")
plot!(abs.(fft(sinaljanelado))[1:50], label="FFT sinal janela")
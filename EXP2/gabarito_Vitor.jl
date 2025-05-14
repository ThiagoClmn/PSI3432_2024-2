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

fa = 10000
Ta = 1/fa
f0 = √3*100
(f0, 2*f0, 5*f0)



N = 1024
σv = sqrt(0.0) # Sem ruído
n = 0:N-1
x = 0.8cos.(2π*f0*n*Ta .- π/5) +  0.1cos.(2π*2*f0*n*Ta .+ π/8) + 0.002cos.(2π*5*f0*n*Ta .+π/2) + σv*randn(N)
plot(n*Ta, x, xlabel = L"$t$ (s)", label = L"x(t)")

file = matopen("sinal.mat","w")
write(file, "x", x)
write(file, "Ta", Ta)
close(file)


s = matread("sinal.mat")
x = s["x"]
Ta = s["Ta"]
N = length(x)


X = fft(x)
plot(n[1:100]*fa/N, abs.(X[1:100]), xlabel=L"$f (Hz)", label = "|X[k]|", line =  marker = (:circle, 3), title = "Janela Retangular")
# Observe que a frequência correspondente à raia k da TDF é k*fa/N


plot(
    n[1:100]*fa/N,
    amp2db.(abs.(X[1:100])),
    xlabel=L"$f (Hz)",
    label = "|X[k]|", marker = (:circle, 3),
    title = "Janela Retangular"
    )

Xh = fft(x.*hamming(N))
Xhdb = amp2db.(abs.(Xh[1:100]))
yt = 0:25:100
plot(
    n[1:100]*fa/N,
    amp2db.(abs.(Xh[1:100])).+60,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|", marker = (:circle, 3),
    title = "Janela de Hamming", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-60)))
    )
# Este truque de somar 60 e mudar a escala é apenas para as barrinhas come


Xb = fft(x.*blackman(N))
Xbdb = amp2db.(abs.(Xb[1:100]))
yt = 0:20:120
plot(
    n[1:100]*fa/N, Xbdb.+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|",
    line = :ste,
    marker = (:circle, 3),
    title = "Janela de Blackman", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )


it = 15:25
Xbdb = amp2db.(abs.(Xb[it]))
yt = 0:20:120
plot(
    n[it]*fa/N,
    Xbdb.+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|",
    line = :stem,
    marker = (:circle, 3),
    title = "Janela de Blackman", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )

k̂ = argmax(Xbdb) + it[1] - 2 # -2 para k valer 0 para sinal DC

f̂0 = k̂*fa/N

((k̂-1)*fa/N, k̂*fa/N)

N1 = 8192
Xb1 = fft([x .* blackman(N); zeros(N1-N)])
k = 0:N1-1
Xbdb1 = amp2db.(abs.(Xb1))
yt = 0:20:120
plot(
    k[1:800]*fa/N1,
    Xbdb1[1:800].+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|",
    line = :stem,
    marker = (:circle, 3),
    title = "Janela de Blackman", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )

it = 120:200
plot(
    k[it]*fa/N1,
    Xbdb1[it].+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|",
    line = :stem,
    marker = (:circle, 3),
    title = "Janela de Blackman", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )

k̂ = argmax(Xbdb1) - 1 # -1 para k valer 0 para sinal DC
Xbdb1[k̂-1:k̂+3]
((k̂-1)*fa/N1, k̂*fa/N1)
2*abs(Xb1[k̂+1]) / sum(blackman(N))

Xf1 = fft([x .* flattop(N); zeros(N1-N)])
k = 0:N1-1
Xfdb1 = amp2db.(abs.(Xf1))
yt = 0:20:120
plot(
    k[1:800]*fa/N1,
    Xfdb1[1:800].+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|",
    line=:stem, marker = (:circle, 3),
    title = "Janela Flattop", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )

k̂ = argmax(Xbdb1) - 1 # -1 para k valer 0 para sinal DC
((k̂-1)*fa/N1, k̂*fa/N1)
2*abs(Xf1[k̂+1]) / sum(flattop(N))
ϕ0 = (angle(Xf1[k̂]),angle(Xf1[k̂+1]))
-π/5

it = 240:400
plot(
    k[it]*fa/N1,
    Xfdb1[it].+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|", line = :stem,
    marker = (:circle, 3),
    title = "Janela Flattop", ytick = yt,
    yaxis = (formatter = yt -> string(Int(yt-80)))
    )

k̂1 = argmax(Xfdb1[it]) + it[1] - 2 # -2 para k valer 0 para sinal DC
Xfdb1[k̂1-1:k̂1+3]
((k̂1-1)*fa/N1, (k̂1+1)*fa/N1)
2*abs(Xf1[k̂1+1]) / sum(flattop(N))
2*abs(Xb1[k̂1+1]) / sum(blackman(N))
(angle(Xf1[k̂1+2]), angle(Xf1[k̂1]))
π/8

it = 600:1000
plot(
    k[it]*fa/N1,
    Xfdb1[it].+80,
    xlabel=L"$f$ (Hz)",
    label = "|X[k]|", line = :stem,
    marker = (:circle, 3), title = "Janela Flattop",
    ytick = yt, yaxis = (formatter = yt -> string(Int(yt-80)))
    )

k̂2 = argmax(Xfdb1[it]) + it[1] - 2 # -2 para k valer 0 para sinal DC
Xfdb1[k̂2-1:k̂2+3]
((k̂2-1)*fa/N1, (k̂2+1)*fa/N1)
2*abs(Xf1[k̂2+1]) / sum(flattop(N))
2*abs(Xb1[k̂2+1]) / sum(blackman(N))
(angle(Xf1[k̂2]), angle(Xf1[k̂2+2]))
π/2




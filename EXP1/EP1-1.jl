#================================================#
# Title : EP1 - Exercício Computacional - PSI 3432																								
# Authors: Thiago da Rocha Calomino Gonçalves, 
# Lucas Gaspar Mendonça															
# Description:	-										
#================================================#

#================================================#
# Packages
#================================================#
import Pkg

Pkg.add("Plots")
using Plots

Pkg.add("MAT")
using MAT

Pkg.add("LaTeXStrings")
using LaTeXStrings

Pkg.add("DSP")
using DSP

Pkg.add("Statistics")
using Statistics


#================================================#
# Parte 2
#================================================#

######### Valores da Parte 1 ###########

# Quantidade de antenas
M = 8

# Velocidade do sinal (m/s)
c = 3e8

######### Valores fornecidos ###########

# Frequência de operação
f = 60e9
Ω0 = 2 * π * f

# Comprimento de onda (λ)
λ = c/f
d = λ/2

# Frequência de amostragem (Hz)
fa = 1e12

# Amplitudes
A0 = 3
A1 = A2 = 1

# Ângulos em graus convertidos para radianos

# θ0 = 20°
θ0 = 20/180 * π

# θ1 = 45°
θ1 = 45/180 * π

# θ2 = -15° 
θ2 = -15/180 * π 


############## Item 1 ##################

# wn  (frequência de corte normalizada)
wn = 2 * π * f/fa

flt = digitalfilter(Lowpass(wn/π), Butterworth(6))

ω = range(0, π, length=500)
H = freqresp(flt, ω)

plot(ω, abs.(H))

# filtragem do sinal 
cd("C:/Users/lucas/Documents/USP/Disciplinas/8_semestre/PES2/EP1")

sinais = matread("sinais.mat")

x_til = sinais["sinal_recebido"]

n = 0:(length(x_til[:, 1])-1)

# xim                           
port_xim = (2) .* cos.(Ω0 * (n/fa))
xim_in = x[:, 1] .* port_xim
xim = filt(flt, xim_in)

plot(n, xim, xlabel =  L"n", label = L"x_{im}")

# xqm
port_xqm = (-2) .* sin.(Ω0 * (n / fa))
xqm_in = x_til[:, 1].* port_xqm
xqm = filt(flt, xqm_in)

plot(n, xqm, xlabel =  L"n", label = L"x_{qm}")

############## Item 2 ##################

# yre
port_yre = (2) .* cos.(Ω0 * (n / fa))
yre_in = x_til .* port_yre
yre = filt(flt, yre_in)

# yim
port_yim = (-2) .* sin.(Ω0 * (n / fa))
yim_im = x_til .* port_yim
yim = filt(flt, yim_in)

z = yre .+ yim .* im

m = 0:(M-1)
v_θ0 = exp.((m.*(Ω0 * d * sin(θ0)/c))im)
w = (1/M)*transpose(conj(v_θ0))

y = transpose(w * transpose(z))

plot(n, real(y), xlabel =  L"n", label = L"\mathbb{R}\{y\}")
plot(n, imag(y), xlabel =  L"n", label = L"\mathbb{Im}\{y\}")

############## Item 3 ##################
Ta = 1e-9
N = round(Int, fa * Ta)

# Número de símbolos transmitidos
qntd_simbol = Int(length(y) / N)

# Calculando o valor médio de cada símbolo
simbols = []
for i in 1:qntd_simbol
    typeof(i)
    symbol_segment = y[(i-1)*N+1:i*N]
    push!(simbols, mean(symbol_segment))
end

# Resultados
println("Número de símbolos transmitidos: ", num_symbols)
println("Valores médios de cada símbolo: ", simbols)

############## Item 4 ##################
symbols_16_qam = [
    -3 + 3im,
    -3 + 1im,
    -3 - 3im,
    -3 - 1im,
    -1 + 3im,
    -1 + 1im,
    -1 - 3im,
    -1 - 1im,
     3 + 3im,
     3 + 1im,
     3 - 3im,
     3 - 1im,
     1 + 3im,
     1 + 1im,
     1 - 3im,
     1 - 1im
]

# Cálculo das distâncias
distances = [abs.(simbol .- symbols_16_qam) for simbol in simbols]

# Seleção das distâncias mínimas
index_symbols_received_message = [argmin(d)-1 for d in distances]

# Indexação dos símbolos tabelados
recon_message = symbols_16_qam[index_symbols_received_message]

# Messagem com os símbolos mapeados
println(recon_message)

############## Item 5 ##################
mat_message = matread("mensagem.mat")
index_symbols_mat_message =  Int64.(vec(mat_message["mensagem"]))

map_message = symbols_16_qam[index_symbols_mat_message]

println(map_message)

# Gráfico de dispersão para os símbolos originais e reconstruídos
scatter(real(map_message), imag(map_message), label="Mensagem Original", color=:blue, markersize = 6)
scatter!(real(recon_message), imag(recon_message), label="Símbolos Reconstruídos", color=:red, markersize = 5)

title!("Comparação de Símbolos no Plano Complexo")
xlabel!("Parte Real")
ylabel!("Parte Imaginária")

############## Item 6 ##################
ruido_σ10 = 10 * randn(size(x_til))
ruido_σ30 = 30 * randn(size(x_til))
ruido_σ50 = 50 * randn(size(x_til))

println(ruido_σ10[1:15,:])
sinal_com_ruido_σ10 = x .+ ruido_σ10
sinal_com_ruido_σ30 = x .+ ruido_σ30
sinal_com_ruido_σ50 = x .+ ruido_σ50

# yre
port_yre = (2) .* cos.(Ω0 * (n / fa))

yre_σ10_in = sinal_com_ruido_σ10 .* port_yre
yre_σ30_in = sinal_com_ruido_σ30 .* port_yre
yre_σ50_in = sinal_com_ruido_σ50 .* port_yre

yre_σ10 = filt(flt, yre_σ10_in)
yre_σ30 = filt(flt, yre_σ30_in)
yre_σ50 = filt(flt, yre_σ50_in)

# yim
port_yim = (-2) .* sin.(Ω0 * (n / fa))
yim_σ10_in = sinal_com_ruido_σ10 .* port_yim
yim_σ30_in = sinal_com_ruido_σ30 .* port_yim
yim_σ50_in = sinal_com_ruido_σ50 .* port_yim

yim_σ10 = filt(flt, yim_σ10_in)
yim_σ30 = filt(flt, yim_σ30_in)
yim_σ50 = filt(flt, yim_σ50_in)

# z's
z_σ10 = yre_σ10 .+ yim_σ10 .* im
z_σ30 = yre_σ30 .+ yim_σ30 .* im
z_σ50 = yre_σ50 .+ yim_σ50 .* im

y_σ10 = transpose(w * transpose(z_σ10))
y_σ30 = transpose(w * transpose(z_σ30))
y_σ50 = transpose(w * transpose(z_σ50))

plot(n, real(y_σ10), xlabel =  L"n", label = L"\mathbb{R}\{y\}")
plot(n, imag(y_σ10), xlabel =  L"n", label = L"\mathbb{Im}\{y\}")

plot(n, real(y_σ30), xlabel =  L"n", label = L"\mathbb{R}\{y\}")
plot(n, imag(y_σ30), xlabel =  L"n", label = L"\mathbb{Im}\{y\}")

plot(n, real(y_σ50), xlabel =  L"n", label = L"\mathbb{R}\{y\}")
plot(n, imag(y_σ50), xlabel =  L"n", label = L"\mathbb{Im}\{y\}")
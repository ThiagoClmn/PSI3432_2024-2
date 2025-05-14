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
# Parte 1
#================================================#
############## Item 1 ##################

function u(theta, d, c)
    return d .* sin.(theta) / c
end

function B(theta, theta_0, M, d, c, omega)
    a = (u(theta, d, c) .- u(theta_0, d, c)) .* omega
    numerador = (1 .- exp.(1im .* M .* a))
    denominador = (1 .- exp.(1im .* a))
    return abs.(1/M .* (numerador ./ denominador))
end

M = 8 # Quantidade de antenas
c = 3e8 # Velocidade do sinal (m/s)

######### Valores fornecidos ###########

# Frequência de operação
f = 60e9
Ω0 = 2 * π * f

# Comprimento de onda (λ)
λ = c / f
d = λ / 2 # distância entre antenas

# Frequência de amostragem (Hz)
fa = 1e12

# Amplitudes
A0 = 3
A1 = A2 = 1

# Ângulos em graus convertidos para radianos

# θ0 = -25°
θ0 = -25 * π / 180

d_array = [λ/4, 2*λ/4, 3*λ/4, λ]
Ω0 = 2 * π * f

θ = range(-π/2, π/2, length=10000)

B_abs = B(θ, θ0, M, d, c, Ω0)

# Criando o gráfico
θ_degrees = θ .* (180 / π)
plot(θ_degrees, B_abs, title="Reforçando o sinal de -25°",
     xlabel=L"\theta (\circ)", ylabel=L"|B(\theta,\theta_{0})|",
     legend=false)

############## Item 2 ##################

# Função para criar subplots
function plot_subplots(θ, B_abs, d_array)
    # Converter θ para graus
    θ_degrees = θ .* (180 / π)
    
    # Ajustar a altura do gráfico usando o parâmetro size
    p = plot(θ_degrees, B_abs[:, 1], label="d = λ/4", xlabel="θ (°)", ylabel="|B(θ,θ₀)|", layout=(4, 1), title="Fenômeno de 'aliasing'", size=(1200, 1200))
    
    plot!(p, θ_degrees, B_abs[:, 2], label="d = λ/2", xlabel="θ (°)", ylabel="|B(θ,θ₀)|", subplot=2)
    plot!(p, θ_degrees, B_abs[:, 3], label="d = 3λ/4", xlabel="θ (°)", ylabel="|B(θ,θ₀)|", subplot=3)
    plot!(p, θ_degrees, B_abs[:, 4], label="d = λ", xlabel="θ (°)", ylabel="|B(θ,θ₀)|", subplot=4)
    
    return p
end

# Gerar os dados para B_abs para cada valor de d em d_array
B_abs = hcat([B(θ, θ0, M, d, c, Ω0) for d in d_array]...)

# Chamar a função de plotagem
p = plot_subplots(θ, B_abs, d_array)





############## Item 3 ##################
M = 8 # antenas
f = 60e9 # Hz
c = 3e8 # m/s
λ = c/f # comprimento de onda
d_array = [λ/4, 2*λ/4]
Ω0 = 2 * π * f

θ = range(-π/2, π/2, length=10000)
θ0 = 90 * (π / 180)

function plot_subplots(θ, B_abs, d_array)
    # Converter θ para graus
    θ_degrees = θ .* (180 / π)
        
    # Ajustar a altura do gráfico usando o parâmetro size
    p = plot(θ_degrees, B_abs[:, 1], label="d = λ/4", xlabel=L"\theta (\circ)",
    ylabel=L"|B(\theta,\theta_{0})|", layout=(2, 1), title="Reforçando o sinal de 90°", size=(900, 900))
        
    plot!(p, θ_degrees, B_abs[:, 2], label="d = λ/2", xlabel=L"\theta (\circ)",
    ylabel=L"|B(\theta,\theta_{0})|", subplot=2)
        
    return p
end

# Gerar os dados para B_abs para cada valor de d em d_array
B_abs = hcat([B(θ, θ0, M, d, c, Ω0) for d in d_array]...)

# Criar o gráfico
p = plot_subplots(θ, B_abs, d_array)



############## Item 4 ##################
# Parâmetros
M = 8
c = 3e8
f = 60
Ω0 = 2 * π * f
λ = c / f
d = λ / 2
θ0 = -25 * (π / 180)
θ2 = range(-π/2, π/2, length=10000)

A_1 = 1
A_2 = 2

# Função para calcular o valor absoluto do beamforming
function beamforming_absolute_value(A_1, A_2, θ2)
    return abs.(A_1 .+ A_2 .* B(θ2, θ0, M, d, c, Ω0))
end

# Calcular ABS_theta_2
ABS_theta_2 = beamforming_absolute_value(A_1, A_2, θ2)

# Converter θ2 para graus
θ2_degrees = θ2 .* (180 / π)

# Plote a amplitude do sinal estimada em função de θ2
plot(θ2_degrees, ABS_theta_2, title="Interferência vindo de diversos lugares",
     xlabel=L"\theta_2 (\circ)", ylabel=L"|y(t)|", legend=false)
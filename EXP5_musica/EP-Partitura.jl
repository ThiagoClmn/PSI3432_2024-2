#================================================#
# PSI 3432 - Processamento de Áudio e Imagem 
# Title : EP - Partitura																							
# Authors: Lucas Gaspar Mendonça, 
# Thiago da Rocha Calomino Gonçalves												
# Description:									
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

Pkg.add("FFTW")
using FFTW

#================================================#
# Parte 1
#================================================#

########### Funções Base ###############

# sinal recebido
cd("C:/Users/lucas/Documents/USP/Disciplinas/8_semestre/PES2/EC5")

arquivo = matread("sinal.mat")
sinal = arquivo["x"]
Ta = arquivo["Ta"]
fa = 1/Ta
n = 0:(length(sinal[:, 1])-1)

f_min = 196.8886 # sol2
Δf_min = f_min*(2.0^(1/12) - 1) 

# N
N_min = Int(round(fa/Δf_min))
N = nextpow(2, N_min)

# M
M_max = (60/320)*fa

# Tomando-se a Janela de Kaiser de M pontos:
# p/ atenuação em lóbulo lateral de 90 dB
Δω = f_min/fa
M_min = Int(round((24π*(120+12))/(155*2π*Δω)+1))

M = Int(round((M_min+M_max)/2))

# Asl - Atenuação do lóbulo lateral
Asl = (60+120)/2

if Asl < 13.26
    β = 0
elseif Asl < 60
    β = 0.76609(Asl-13.26)^(0.4)+0.09834(Asl-13.26)
elseif Asl < 120
    β = 0.12438(Asl+6.3)
end
M = 3000
N = 6000

wk = kaiser(M, β/π)

m_window_signal = []
peak_freqs = Vector{Float64}(undef, 1 + div(length(sinal), length(wk)))

# Função para calcular a STFT
function stft(signal, window, nfft, fa)
    # Tamanho do sinal e janela
    signal_length = length(signal)
    window_length = length(window)
    
    # Número de janelas
    num_windows = 1 + div(signal_length, window_length)

    # Matriz para armazenar o resultado da STFT
    stft_matrix = Matrix{Float64}(undef, num_windows, div(nfft, 2))
    
    # Iterar sobre o sinal com deslocamentos de hop_size
    for i in 1:num_windows
        # Índices do segmento
        start_idx = (i - 1) * window_length + 1
        end_idx = start_idx + window_length - 1
        
        if end_idx > signal_length
            m_window_signal = vcat(sinal[start_idx:length(sinal)], zeros(end_idx - signal_length))
        else
            m_window_signal = sinal[start_idx:end_idx]
        end
        
        # Extraindo o segmento do sinal
        segment = (m_window_signal .* window)
        padded_segment = vcat(segment, zeros(nfft - window_length))

        # Calculando a TDF do segmento
        spectrum = fft(padded_segment)

        # Adicionando à matriz STFT
        stft_matrix[i, :] = abs.(spectrum[1:div(nfft, 2)])
        magnitude = abs.(spectrum[1:div(N, 2)]) 

        # Encontrando o maior pico no espectro na faixa de interesse
        min_idx = floor(Int, 100 * nfft / fa) + 1
        max_idx = ceil(Int, 1000 * nfft / fa)

        peak_val, peak_idx = findmax(magnitude[min_idx:max_idx])
        peak_idx += min_idx - 1

        peak_freq = (peak_idx-1) * fa / N
        peak_freqs[i] = 2*peak_freq
        print(peak_freq)
    end
    return stft_matrix # Organizar resultados em colunas

end

# Tamanho da FFT para cada segmento
nfft = N                            
# Calculando a STFT
stft_result = stft(sinal, wk, nfft, fa)

#print(sinal)
frequencies = range(0, stop=fa/2, length=size(stft_result, 1))
tempo = 1:size(stft_result, 2) 

# Plot do espectrograma com frequências no eixo y e tempo no eixo x
p = heatmap(frequencies, tempo, abs.(stft_result'), xlabel="Tempo (Janelas)", ylabel="Frequências (Hz)", title="STFT - Magnitude")
display(p)

print(peak_freqs)

# Define as notas musicais com base no padrão G2 = 196.886 Hz
function nota_musical(frequencia)
    # Frequência de referência para o G2 (196.886 Hz)
    G2 = 196.886
    # Calcula o número de semitons em relação ao G2
    semitons = round(Int, 12 * log2(frequencia / G2))
    
    # Lista de notas em uma oitava
    notas = ["G", "G#", "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#"]
    
    # Calcula a oitava e a nota
    nota = notas[(semitons % 12) + 1]
    oitava = 3 + div(semitons, 12)  # Começa da oitava 2 para o G2
    
    # Retorna o nome da nota e a oitava
    return "$nota$oitava"
end

# Converte cada frequência para a nota correspondente
notas = [nota_musical(f) for f in peak_freqs]

# Exibe o resultado
println(notas)



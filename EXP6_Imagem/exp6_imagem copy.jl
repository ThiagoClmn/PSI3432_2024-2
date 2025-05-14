##### Dependências
import Pkg;

#Pkg.add("DSP");
using DSP;

#Pkg.add("MAT");
using MAT;

#Pkg.add("Images");
using Images;

#Pkg.add("ImageIO");
using ImageIO;

#Pkg.add("FileIO");
using FileIO;

#Pkg.add("ImageView");
using ImageView;

#Pkg.add("FFTViews");
using FFTViews;

#Pkg.add("Plots");
using Plots

#Pkg.add("FFTW");
using FFTW;


####== Importar imagem ==####
i_aux = load("./cassini-interference.tif")
img = convert(Array{Float32}, i_aux)

# Deslocar a origem da transformada para o centro
centered_img = (-1) .^ ((1:size(img, 1))' .+ (1:size(img, 2))) .* img

# Calcular a TDF-2D
spectrum = fft(centered_img)

# Normalizar o módulo do espectro para visualização
log_spectrum = Gray.(log.(abs.(fftshift(spectrum)) .+ 1) ./ maximum(log.(abs.(spectrum) .+ 1)))

# Exibir o espectro
plot(log_spectrum, title="Espectro em Escala Log", color=:gray, size=(600, 600))

# Criar uma matriz para o filtro notch
filter = ones(size(spectrum))

# Rejeitar regiões identificadas no espectro
notch_width = 11  # Largura de rejeição
cx, cy = size(spectrum) ./ 2  # Centro do espectro

# Coordenadas dos pontos a serem rejeitados (exemplo genérico)
filter[Int(cx) .- notch_width:Int(cx) .+ notch_width, Int(cy)-50:Int(cy)+50] .= 0
filter[Int(cx) .- notch_width:Int(cx) .+ notch_width, Int(cy)+50:Int(cy)-50] .= 0

# Certificar que o filtro é simétrico
filter .= fftshift(filter)


# Filtrar no domínio da frequência
filtered_spectrum = spectrum .* filter

# Calcular a TDF inversa e desfazer o deslocamento
filtered_img = real(ifft(filtered_spectrum))
filtered_img .= (-1) .^ ((1:size(img, 1))' .+ (1:size(img, 2))) .* filtered_img

# Visualizar a imagem filtrada
plot(filtered_img, title="Imagem Filtrada", color=:gray, size=(600, 600))











####== Funções de ImageShifter / ZeroPadding ==####
function ImageShifter(Image)
    n1, n2 = size(Image)
    s_img = similar(Image)
    for i in 1:n1
        for j in 1:n2
            s_img[i, j] = ((-1)^(i + j)) * Image[i, j]
        end
    end
    return s_img
end

function ZeroPadding(Image, L)
    n1, n2 = size(Image)
    zero_pad = zeros(Float32, n1 + 2*L, n2 + 2*L)
    zero_pad[(L+1):(n1+L), (L+1):(n2+L)] .= Image
    return zero_pad
end

##== Dimensões ==##
M_1, M_2 = size(img)
N_1 = 2*M_1
N_2 = 2*M_2

# Zero padding
padded_img = ZeroPadding(img, 2)

# Multiplicar a imagem por (-1)^(n1+n2)
shifted_img = ImageShifter(padded_img)

# FFT da imagem deslocada
fft_img = fft(shifted_img)





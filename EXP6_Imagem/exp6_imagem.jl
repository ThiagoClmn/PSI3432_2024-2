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
using Plots;

using FFTW;


####== Importar imagem ==####
i_aux = load("./cassini-interference.tif")
img = convert(Array{Float32}, i_aux)

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
# imshow(padded_img)

# Multiplicar a imagem por (-1)^(n1+n2)
shifted_img = ImageShifter(padded_img)
save("output.png", shifted_img)

# FFT da imagem deslocada
fft_img = fft(shifted_img)





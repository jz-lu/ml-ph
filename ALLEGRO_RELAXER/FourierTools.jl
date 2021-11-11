### Implementation of GSFE for graphene homo-bilayers.

# 2D Fourier Decomposition of matrix M
function fourierDecompose(M::Array{Float64,2}, X::Array{Float64,2}, Y::Array{Float64,2}, k_list::Array{Float64,2})

    #FFTW.set_num_threads(Sys.CPU_CORES)

     # check for 2D input
     if ndims(M) != 2
         throw(ArgumentError("Array M input to fourierDecompose is not 2D!"))
     end

     coeffs = complex(zeros(size(k_list,2),1))
     for k_idx = 1:size(k_list,2)
         coeffs[k_idx] = fourierInnerProduct(M,X,Y,k_list[1:2,k_idx])
         #=
         @printf("(%d, %d, %d) %lf + %lf im \n",
                    k_list[3,k_idx],k_list[4,k_idx],k_list[5,k_idx],
                    real(coeffs[k_idx]),imag(coeffs[k_idx]))
        =#
     end
     return coeffs

     #show(STDOUT, "text/plain", coeffs)

     #N = size(M,1)

     # pick the

     # do the 2D FFT
     #M_FFT = (fft(M)/N^2)
     #prec = 10;
     #M_FFT = round.(M_FFT,prec)

     #pcolormesh(real((M_FFT)))

     #writedlm("./temp/test.txt", vcat(real(M_FFT),imag(M_FFT)))
     #=
     open(string("./temp/",filename), "w") do f
          for i in 1:N
               if (i == 1)
                    write(f,"$N $N 0 \n")
               end
               for j in 1:N
                    temp_r = real(M_FFT[i,j])
                    write(f, "$temp_r   ")
               end
               write(f,"\n")
          end
          for i in 1:N
               if (i == 1)
                    write(f,"$N $N 1 \n")
               end
               for j in 1:N
                    temp_i = imag(M_FFT[i,j])
                    write(f, "$temp_i   ")
               end
               write(f,"\n")
          end
     end
     =#

end

function fourierInnerProduct(M::Array{Float64,2}, X::Array{Float64,2}, Y::Array{Float64,2}, k::Array{Float64,1})

    Nx = size(M,1)
    Ny = size(M,2)

    phase = -1im*(k[1]*X + k[2]*Y)

    inner_prod = sum( sum((exp.(phase)).*M) )/(Nx*Ny)
    return inner_prod

end

function fourierGetSampling(A::Array{Float64,2}, max_p::Int64)

    # 90 degree CCW rotation
    R = [0 -1; 1 0]

    B = [0.0 0.0; 0.0 0.0]

    B[:,1] = 2*pi*(R*A[:,2])/(dot(A[:,1],R*A[:,2]))
    B[:,2] = 2*pi*(R*A[:,1])/(dot(A[:,2],R*A[:,1]))
    k1 = B[:,1]
    k2 = B[:,2]
    k3 = -(B[:,1]+B[:,2])

    #println(B[:,1])
    #println(B[:,2])

    max_i = max_p
    #max_rot = 6; # only need for checking sym
    max_rot = 1;
    max_k = max_rot*(max_i+1)*(max_i)/2 + 1


    k_list = zeros(5,max_k)
    k_idx = 1
    for i_idx = 0:max_i
        for j_idx = 0:max((i_idx-1),0)
            k_here = i_idx*k1 + j_idx*k2
            for rot_idx = 1:max_rot
                ang = 2.0*pi*(rot_idx-1)/6
                rot_mat = [cos(ang) -sin(ang); sin(ang) cos(ang)]
                k_list[1:2,k_idx] = rot_mat*k_here
                k_list[3:5,k_idx] = [i_idx j_idx rot_idx]
                k_idx = k_idx+1
            end
        end
    end

    return k_list

end

function fourierGetSamplingFull(A::Array{Float64,2}, max_p::Int64)

    # 90 degree CCW rotation
    R = [0 -1; 1 0]

    B = [0.0 0.0; 0.0 0.0]

    B[:,1] = 2*pi*(R*A[:,2])/(dot(A[:,1],R*A[:,2]))
    B[:,2] = 2*pi*(R*A[:,1])/(dot(A[:,2],R*A[:,1]))
    k1 = B[:,1]
    k2 = B[:,2]
    k3 = -(B[:,1]+B[:,2])

    #println(B[:,1])
    #println(B[:,2])

    max_i = max_p
    max_rot = 6; # only need for checking sym
    #max_rot = 1;
    max_k = convert(Int,max_rot*((max_i+1)*(max_i)/2 + 1))


    k_list = zeros(5,max_k)
    k_idx = 1
    for i_idx = 0:max_i
        for j_idx = 0:max((i_idx-1),0)
            k_here = i_idx*k1 + j_idx*k2
            for rot_idx = 1:max_rot
                ang = 2.0*pi*(rot_idx-1)/6
                rot_mat = [cos(ang) -sin(ang); sin(ang) cos(ang)]
                k_list[1:2,k_idx] = rot_mat*k_here
                k_list[3:5,k_idx] = [i_idx j_idx rot_idx]
                k_idx = k_idx+1
            end
        end
    end

    return k_list

end

function fourierGetSamplingConfig(A::Array{Float64,2}, max_p::Int64)

    # 90 degree CCW rotation
    R = [0 -1; 1 0]

    B = [0.0 0.0; 0.0 0.0]

    B[:,1] = 2*pi*(R*A[:,2])/(dot(A[:,1],R*A[:,2]))
    B[:,2] = 2*pi*(R*A[:,1])/(dot(A[:,2],R*A[:,1]))
    k1 = B[:,1]
    k2 = B[:,2]
    k3 = -(B[:,1]+B[:,2])

    #println(B[:,1])
    #println(B[:,2])

    max_i = max_p
    #max_rot = 6; # only need for checking sym
    max_rot = 1;
    max_k = convert(Int,max_rot*((max_i+1)*(max_i)/2 + 1))


    k_list = zeros(5,max_p^2)
    k_idx = 1
    for i_idx = 1:max_i
        for j_idx = 1:max_i#max((i_idx-1),0)
            k_here = (i_idx-1)*k1 + (j_idx-1)*k2
            for rot_idx = 1:max_rot
                ang = 2.0*pi*(rot_idx-1)/6
                rot_mat = [cos(ang) -sin(ang); sin(ang) cos(ang)]
                k_list[1:2,k_idx] = rot_mat*k_here
                k_list[3:5,k_idx] = [i_idx j_idx rot_idx]
                k_idx = k_idx+1
            end
        end
    end

    return k_list

end
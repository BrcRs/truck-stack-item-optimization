using LinearAlgebra
using JuMP
using Cbc

model = Model(Cbc.Optimizer)

@variable(model, zetaT[1:nbPlannedTrucks] >= 0)
@variable(model, zetaE[1:nbExtraTrucks] >= 0)

@variable(model, SS[1:nbStacks] >= 0)
@variable(model, SP[1:nbStacks] >= 0)
@variable(model, SK[1:nbStacks] >= 0)
@variable(model, SPD[1:nbStacks] >= 0)
@variable(model, SU[1:nbStacks] >= 0)
@variable(model, SO[1:nbStacks] >= 0)
@variable(model, SXe[1:nbStacks] >= 0)
@variable(model, SXo[1:nbStacks] >= 0)
@variable(model, SYe[1:nbStacks] >= 0)
@variable(model, SYo[1:nbStacks] >= 0)
@variable(model, SZe[1:nbStacks] >= 0)
@variable(model, betaM[1:nbStacks] >= 0)
@variable(model, betaP[1:nbStacks] >= 0)
@variable(model, nu[1:nbStacks] >= 0)
@variable(model, tau[1:nbStacks] >= 0)
@variable(model, phi[1:nbStacks] >= 0)
@variable(model, SG[1:nbStacks, 1:nbPlants] >= 0)

@variable(model, TI[1:nbTrucks, 1:nbItems], lower_bound = 0, upper_bound = 1)
@variable(model, R[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1)
@variable(model, Theta[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1)
@variable(model, S[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1)
@variable(model, Z[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1)
@variable(model, Omega[1:nbStacks, 1:nbTrucks, 1:nbItems], lower_bound = 0, upper_bound = 1)
@variable(model, ST[1:nbStacks, 1:nbTrucks], lower_bound = 0, upper_bound = 1)
@variable(model, IOV[1:nbItems], lower_bound = 0, upper_bound = 1)
@variable(model, mu[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, eta[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, xi[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, chi[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, r[1:nbStacks], lower_bound = 0, upper_bound = 1)

@variable(model, sigma1[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, sigma2[1:nbStacks], lower_bound = 0, upper_bound = 1)
@variable(model, sigma3[1:nbStacks], lower_bound = 0, upper_bound = 1)

@variable(model, Psi[1:nbStacks, 1:nbTrucks], lower_bound = 0)
@variable(model, Q[1:nbStacks, 1:nbItems], lower_bound = 0)
@variable(model, H[1:nbStacks, 1:nbItems], lower_bound = 0)
@variable(model, V[1:nbStacks, 1:nbItems], lower_bound = 0)
@variable(model, W[1:nbStacks, 1:nbItems], lower_bound = 0)
@variable(model, Gl[1:nbStacks, 1:nbItems], lower_bound = 0)
@variable(model, Gr[1:nbStacks, 1:nbItems], lower_bound = 0)

@variable(model, lambda[[1:nbStacks * (nbStacks+1)/2]], lower_bound = 0)

# MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

MZ = max(IS...) + 1

MQ = max(IU...) + 1

MH = max(IP...) + 1

MV = max(IK...) + 1

MW = max(IPD...) + 1

Gr = max(SO...) + 1

Gl = max(IO...) + 1

Meta = nbTrucks * Mtau/10

MTE = max(TE...)

MTKE = max(TKE...)

MTGE = max(TGE...)

MTW = max(TW...) + 1

MPsi = nbItems

Mlambda = 2 * TL .+ 1 # TODO check if Mlambda big enough

SZo = 0

MST = 2

MOmega = 2

Mmu = 2

Mtau = 10

Xi1 = falses(convert(Int, nbStacks*(nbStacks+1)/2), nbStacks)
Xi2 = falses(convert(Int, nbStacks*(nbStacks+1)/2), nbStacks)

epsilon = 0.001

function fillXi1!(Xi::BitArray)
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        Xi[m:m+n-1, i] .= 1
        m = m + n
        n = n - 1
        i = i + 1
    end
end

function identityMat(n)
    mat = falses(n, n)
    for i in 1:n
        mat[i,i] = 1
    end

    return mat
end

function fillXi2!(Xi::BitArray)
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        print(m, ":")
        println(n)
        # Xi[m:m+n-1, i] .= 1
        # Xi[m:m+n-1, i:size(Xi)[2]] = diagm([1 for _ in 1:n])
        # Xi[m:m+n-1, i:size(Xi)[2]] = identityMat(n)
        Xi[m:m+n-1, i:size(Xi)[2]] = I(n)
        m = m + n
        n = n - 1
        i = i + 1
        println()
    end
end

function ones(n::Int)
    return ones(Int8, n, 1)
end

fillXi1!(Xi1)
fillXi2!(Xi2)

@constraint(model, cZetaT1, -zetaT >= -ones(size(zetaT)[1]))
@constraint(model, cZetaT2, -zetaT >= -TIT * ones(size(zetaT)[1]))

@constraint(model, cZetaE1, -zetaT >= -ones(size(zetaT)[1]))
@constraint(model, cZetaE2, -zetaT >= -TIE * ones(size(zetaT)[1]))

@constraint(model, cTI_F, TI <= F)

@constraint(model, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[1]) <= ones(size(transpose(TI))[1]))

@constraint(model, cTI_TP_IP, transpose(TI) * TP == IP)

@constraint(model, cR_Theta_MI4, R <= Theta * MI4)

@constraint(model, cR_TI_TU, -MI4*(1-Theta) <= R - (transpose(TI) * TU) <= MI4* (1-Theta))

@constraint(model, cR_1_IU, R * ones(size(R)[1]) == IU)

@constraint(model, cTheta_1_1, Theta * ones(size(Theta)[1]) >= ones(size(Theta)[1]))

@constraint(model, cIDE_TI_TDA_IDL, IDE <= transpose(TI) * TDA <= IDL)

@constraint(model, cZ_S_MZ, Z <= S * MZ)

@constraint(model, cPsi_Omega, Psi = transpose(vcat([transpose(Omega[:,i,:]*ones(size(Omega)[1])) for i in 1:size(Omega)[2]]...)))

for j in 1:size(Omega)[2]
    @constraint(model, cOmega_S_MOmega, Omega[:, j, :] <= S * MOmega)
end
for i in 1:size(Omega)[1]
    @constraint(model, cOmega_TI_S, -(1-S)*MOmega <= Omega[i,:,:] - transpose(TI) <= (1 - S)*MOmega)
end

@constraint(model, cPsi_ST_MPsi, Psi <= St * MPsi)

# @constraint(model, )

print(model)
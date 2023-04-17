module TSIInstanceLoader

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

MZ = max(IS...) + 1.0

MQ = max(IU...) + 1.0

MH = max(IP...) + 1.0

MV = max(IK...) + 1.0

MW = max(IPD...) + 1.0

Gr = max(SO...) + 1.0

Gl = max(IO...) + 1.0

Meta = nbTrucks * Mtau/10

MTE = max(TE...)

MTKE = max(TKE...)

MTGE = max(TGE...)

MTW = max(TW...) + 1.0

MPsi = nbItems

Mlambda = 2.0 * TL .+ 1.0 # TODO check if Mlambda big enough

SZo = 0.0

MST = 2.0

MOmega = 2.0

Mmu = 2.0

Mtau = 10.0

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

@constraint(model, cPsi_Omega, Psi = hcat([Omega[:,i,:]*ones(size(Omega)[1]) for i in 1:size(Omega)[2]]...))

for j in 1:size(Omega)[2]
    @constraint(model, cOmega_S_MOmega, Omega[:, j, :] <= S * MOmega)
end
for i in 1:size(Omega)[1]
    @constraint(model, cOmega_TI_S, -(1-S)*MOmega <= Omega[i,:,:] - transpose(TI) <= (1 - S)*MOmega)
end

@constraint(model, cPsi_ST_MPsi, Psi <= St * MPsi)

@constraint(model, cPsi_S_ST, -(1-ST)*MPsi <= Psi - hcat([S * ones(size(S)[1]) for i in 1:size(Psi)[2]]) <= (1-ST)*MPsi)

@constraint(model, cS_IS_Z, S * IS = Z * ones(size(Z)[1]))

@constraint(model, cZ_SS_S, -MZ * (1-S) <= Z - hcat([SS for i in 1:size(Z)[2]]) <= MZ * (1-S))

@constraint(model, cQ_S, Q <= S * MQ)

@constraint(model, cS_IU_Q, S * IU = Q * ones(size(Q)[1]))

@constraint(model, cQ_SU_S, -MQ * (1-S) <= Q - hcat([SU for i in 1:size(Q)[2]]) <= MQ * (1-S))


@constraint(model, cH_S, H <= S * MH)

@constraint(model, cS_IP_H, S * IP = H * ones(size(H)[1]))

@constraint(model, cQ_SP_S, -MH * (1-S) <= H - hcat([SP for i in 1:size(H)[2]]) <= MH * (1-S))


@constraint(model, cV_S, V <= S * MV)

@constraint(model, cS_IK_V, S * IK = V * ones(size(V)[1]))

@constraint(model, cV_SK_S, -MV * (1-S) <= V - hcat([SK for i in 1:size(V)[2]]) <= MV * (1-S))


@constraint(model, cW_S, W <= S * MW)

@constraint(model, cS_IPD_W, S * IPD = W * ones(size(W)[1]))

@constraint(model, cW_SPD_S, -MW * (1-S) <= W - hcat([SPD for i in 1:size(W)[2]]) <= MW * (1-S))


@constraint(model, cGl_S, Gl <= S * MG)
@constraint(model, cGr_S, Gr <= S * MG)

@constraint(model, Gl * ones(size(Gl)[1]) == Gr * ones(size(Gr)[1]))

@constraint(model, cGr_SO_S, -MG * (1-S) <= Gr - hcat([SO for i in 1:size(Gr)[2]]) <= MG * (1-S))
@constraint(model, cGl_IOV_S, -MG * (1-S) <= Gl - hcat([IOV for i in 1:size(Gl)[2]]) <= MG * (1-S))

@constraint(model, cSXe_SXo_SL_SO, SXe - SXo = SL + SO * MTL)
@constraint(model, cSYe_SYo_SW_SO, SYe - SYo = SW + SO * MTW)

@constraint(model, cSXe_SXo_SW_SO, SXe - SXo = SW + (1 - SO) * MTW)
@constraint(model, cSYe_SYo_SL_SO, SYe - SYo = SL + (1 - SO) * MTL)

@constraint(model, cSZe_S_IH, SZe = S* IH)

@constraint(model, cSXe_ST_TL, SXe <= ST * TL)
@constraint(model, cSYe_ST_TW, SYe <= ST * TW)
@constraint(model, cSZe_ST_TH, SZe <= ST * TH)

@constraint(model, cSXo_SXo, vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo <= SXo)

@constraint(model, cXi2SXo_Xi1SXe_betaM_betaP, Xi2 * SXo - Xi1 * SXe - betaM + betaP == -epsilon)

@constraint(model, cbetaM_lambda, betaM <= lambda * Mlambda)

@constraint(model, betaP <= (1-lambda)*Mlambda)

@constraint(model, cmu_betaM, (1-mu) <= betaM * Mmu)

@constraint(model, cXi2ST_Xi1ST_nu, (Xi2 * ST - Xi1 * ST) * diagm([i for i in 1:size(nu)[2]]) == nu)

@constraint(model, ctau_phi_nu, tau - phi <= (nu - nbTrucks) * Mtau)

@constraint(model, ctau_nu, tau >= (nu-nbTrucks)*Mtau/10)

@constraint(model, ctau_eta, tau <= eta * Meta)

@constraint(model, cphi_eta, phi <= (1 - eta)*Meta)

@constraint(model, cXi1SYe_Xi2SYo, Xi1 * SYe <= Xi2 * SYo + xi * MTW + (tau + phi) * MTW + (1 - mu) * MTW)
@constraint(model, cXi2SYe_Xi1SYo, Xi2 * SYe <= Xi1 * SYo + (1-xi) * MTW + (tau + phi) * MTW + (1 - mu) * MTW)

@constraint(model, cXi1SU_Xi2SU, Xi1 * SU * TE <= Xi2 * SU * TE + (tau + phi) * MTE)

@constraint(model, cXi1SU_Xi2SU_chi, Xi1SU - Xi2SU >= chi * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

@constraint(model, cXi2SU_Xi1SU_chi, Xi2SU - Xi1SU >= (1-chi) * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

@constraint(model, cXi2SK_Xi1SK, Xi2*SK*TKE >= Xi1*SK*TKE - (1 - r) * MTKE - (tau + phi) * MTKE)

@constraint(model, cXi1SK_Xi2SK_chi, Xi1*SK*TKE - Xi2*SK*TKE >= chi*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)
@constraint(model, cXi2SK_Xi1SK_chi, Xi2*SK*TKE - Xi1*SK*TKE >= (1-chi)*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)

@constraint(model, cXi2SG_Xi1SG, Xi2*SG*TGE >= Xi1*SG*TGE - (tau + phi)*MTGE - (1 - sigma3) * MTGE)

@constraint(model, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 >= 1)



print(model)

end